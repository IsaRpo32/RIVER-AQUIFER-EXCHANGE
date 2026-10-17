"""
Parallel MODFLOW Simulation for River-Aquifer Exchange

This script runs parallel MODFLOW-NWT simulations to model river-aquifer
exchange during flood events. It uses the Streamflow Routing (SFR2) package
to simulate dynamic surface water-groundwater interaction.

Usage:
    python parallel_modflow.py

Requirements:
    - MODFLOW-NWT executable in ../executables/ (MODFLOW-NWT.exe on Windows,
      mfnwt on Linux)
    - Input CSV in ../data/input/ (built by build_input_csv.py)

Time discretization, durations and exported metrics were validated with the
scripts in ../validation/ (see validation/README.md).
"""

import os
import sys

# Add parent directory to path for src imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import glob
import json
import logging
import multiprocessing as mp
import shutil
import threading
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple

import flopy
import numpy as np
import pandas as pd

from src.modflow_utils import (
    generate_topo,
    generate_river_stress_period_data,
    hydrograph_cooper_rorabaugh,
    extract_sfr_data,
    analyze_peak_integrals,
)
from src import dimensionless_numbers as dn

ROOT = Path(__file__).resolve().parents[1]
DATA_INPUT = ROOT / "data" / "input"
DATA_OUTPUT = ROOT / "data" / "output"
DEFAULT_EXE = str(ROOT / "executables" / ("MODFLOW-NWT.exe" if os.name == "nt" else "mfnwt"))
INPUTS = DATA_INPUT / "model_inputs.csv"
ENSEMBLES = ("main", "high_connectivity")


def load_inputs(ensemble: str = "main", original_n: bool = False) -> pd.DataFrame:
    """
    Rows of data/input/model_inputs.csv for one ensemble ("main" or "high_connectivity").
    The two ensembles share IDMERGE values, so they are always run separately.
    original_n=True puts back the earlier Manning n (ten times larger), kept in n_original.
    """
    if ensemble not in ENSEMBLES:
        raise ValueError(f"ensemble must be one of {ENSEMBLES}, not {ensemble!r}")
    d = pd.read_csv(INPUTS)
    d = d[d["ensemble"] == ensemble].reset_index(drop=True)
    if original_n:
        d["n"] = d["n_original"]
    return d


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('simulation_parallel.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Thread-safe results consolidation
results_lock = threading.Lock()
progress_lock = threading.Lock()

def create_simulation_directory(simulation_id: int, base_dir: str = str(DATA_OUTPUT), subdir: str = "simple_model") -> str:
    """Create unique directory for each simulation"""
    path_results = os.path.join(base_dir, subdir, f"run_{simulation_id:03d}", "")
    Path(path_results).mkdir(parents=True, exist_ok=True)
    return path_results

def cleanup_simulation_files(sim_dir: str, keep_essential: bool = True):
    """
    Clean up simulation files to save disk space
    
    Parameters:
    -----------
    sim_dir : str
        Simulation directory path
    keep_essential : bool
        If True, keeps only essential files (QHidrom.dat, log files)
        If False, removes entire directory
    """
    if not keep_essential:
        # Remove entire directory
        if os.path.exists(sim_dir):
            shutil.rmtree(sim_dir)
            logger.debug(f"Removed directory: {sim_dir}")
        return
    
    # Keep only essential files
    essential_files = [
        'QHidrom.dat',          # Hydrogram data (essential)
        '*.log',                # Log files
        '*.lst',                # List files (small, useful for debugging)
        '*.png'                 # Generated figures (keep for analysis)
    ]
    
    # Files to remove (large files)
    files_to_remove = [
        # '*.hds',                # Head files (large)
        '*.cbc',                # Cell-by-cell flow files (large)
        '*.dis',                # Discretization files (can be regenerated)
        '*.bas',                # Basic package files (can be regenerated)
        '*.upw',                # UPW package files (can be regenerated)
        '*.nwt',                # NWT solver files (can be regenerated)
        '*.sfr',                # SFR package files (except .out)
        '*.ghb',                # GHB package files (can be regenerated)
        '*.oc'                  # Output control files (can be regenerated)
    ]
    
    try:
        # Remove large files
        files_removed = 0
        space_saved = 0
        
        for pattern in files_to_remove:
            for filepath in glob.glob(os.path.join(sim_dir, pattern)):
                if os.path.exists(filepath):
                    file_size = os.path.getsize(filepath)
                    os.remove(filepath)
                    files_removed += 1
                    space_saved += file_size
                    
        if files_removed > 0:
            space_mb = space_saved / (1024 * 1024)
            logger.debug(f"Cleaned {sim_dir}: removed {files_removed} files, saved {space_mb:.1f} MB")
            
    except Exception as e:
        logger.warning(f"Error cleaning up {sim_dir}: {e}")

DEFAULT_TIME_CFG = {
    'mode': 'adaptive',      # 'adaptive' (event-based steps) or 'courant' (old uniform Courant step)
    'n_hyd': 200,            # steps per flood duration t_d during the event
    'courant_sfr': 1.0,      # max wave Courant per SFR routing sub-step (celerity*dt/(numtim*dx))
    'flood_margin': 1.0,     # fine steps continue this many t_d after the flood ends
    'tsmult': 1.1,           # geometric growth of dt during the recession
    'dt_max_hours': 6.0,     # cap of dt during the recession
    'courant_factor': 1.0,   # only for mode='courant': multiplies cteCourant
    'nstp_uniform': None,    # only for mode='courant': fixed number of uniform transient steps (overrides Courant)
    'ghb_cond_factor': 1.0,  # multiplies the GHB conductance (sensitivity tests only)
    'duration_days': None,   # None -> duration_by_basin(da_km2)
    'save_heads': False,     # .hds every step (~250 MB/run); not used by any exported metric
    'save_cbc': True,        # cell-by-cell budget file (.cbc); deleted by the cleanup anyway
    'print_budget': False,   # write the volumetric budget to the .list at every step (mass-balance checks)
}

def duration_by_basin(da_km2: float) -> float:
    """
    Transient length [days]: 90 (1000 km2), 150 (10000 km2), 400 (100000 km2), interpolated in log10(area).
    Set with test_duration.py so that the net exchange falls below 1 % of its peak before the end
    (slowest tested rows: 19 d, 77 d and 322 d for the three basin sizes).
    """
    return float(np.interp(np.log10(da_km2), [3, 4, 5], [90, 150, 400]))

def build_adaptive_periods(duration: float, t_delay: float, t_d: float, v_max: float,
                           reach_len: float, numtim: int, cfg: Dict) -> Tuple[List, List, List]:
    """
    Transient stress periods (after the steady one):
      SP1: uniform fine dt over delay + flood + margin
      SP2: geometric growth (tsmult) from the fine dt up to dt_max
      SP3: uniform dt_max until the end
    Fine dt = min(t_d/n_hyd, SFR routing limit with wave celerity 5/3*v_max).
    """
    celerity = 5/3 * v_max
    dt_route = cfg['courant_sfr'] * numtim * reach_len / celerity
    dt_fine = min(t_d / cfg['n_hyd'], dt_route)
    dt_max = max(cfg['dt_max_hours'] * 3600, dt_fine)
    r = cfg['tsmult']

    len1 = min(t_delay + (1 + cfg['flood_margin']) * t_d, duration)
    nstp1 = int(np.ceil(len1 / dt_fine))
    perlen, nstp, tsmult = [len1], [nstp1], [1.0]
    rest = duration - len1
    if rest <= 0:
        return perlen, nstp, tsmult

    dt0 = len1 / nstp1
    m = int(np.ceil(np.log(dt_max / dt0) / np.log(r))) if r > 1 and dt_max > dt0 else 0
    len_grow = dt0 * r * (r**m - 1) / (r - 1) if m > 0 else 0.0
    if len_grow >= rest:
        # Growth alone fills the remaining time
        n_g = int(np.ceil(np.log(1 + rest * (r - 1) / dt0) / np.log(r)))
        perlen.append(rest); nstp.append(n_g); tsmult.append(r)
    else:
        if m > 0:
            perlen.append(len_grow); nstp.append(m); tsmult.append(r)
        perlen.append(rest - len_grow); nstp.append(int(np.ceil((rest - len_grow) / dt_max))); tsmult.append(1.0)
    return perlen, nstp, tsmult

def step_end_times(perlen: List, nstp: List, tsmult: List) -> np.ndarray:
    """End time of every time step, measured from the start of the first period given."""
    times, t_start = [], 0.0
    for L_p, n_p, r in zip(perlen, nstp, tsmult):
        if r == 1.0:
            dts = np.full(n_p, L_p / n_p)
        else:
            dt1 = L_p * (r - 1) / (r**n_p - 1)
            dts = dt1 * r**np.arange(n_p)
        times.append(t_start + np.cumsum(dts))
        t_start += L_p
    return np.concatenate(times)

def interp_peak(t: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """Time and value of the maximum, refined with a parabola through the max and its two neighbours."""
    i = int(np.argmax(y))
    if i == 0 or i == len(y) - 1:
        return float(t[i]), float(y[i])
    a, b, c = np.polyfit(t[i-1:i+2] - t[i], y[i-1:i+2], 2)
    if a >= 0:
        return float(t[i]), float(y[i])
    dt_v = float(np.clip(-b / (2 * a), t[i-1] - t[i], t[i+1] - t[i]))
    return float(t[i] + dt_v), float(a * dt_v**2 + b * dt_v + c)

def flood_wave_metrics(t_in: np.ndarray, H_in: np.ndarray, Times: np.ndarray,
                       Q_out_reach: np.ndarray, H_out: np.ndarray) -> Dict:
    """
    Peak delay and stage attenuation as in the manuscript (Eq. 16-17):
      t_p,in, H_p,in  -> prescribed upstream stage hydrograph H_b(t) (as injected, beta floor included)
      t_p,out         -> peak of the simulated outlet stage hydrograph; read from the inflow of the outlet
                         reach, from which SFR (ICALC=1) computes the depth with Manning, so both peak at the
                         same instant, but flow is printed with 5 digits and depth only to 1 mm
      H_p,out         -> maximum simulated outlet depth (printed resolution 1 mm)
    All times measured from the start of the transient period.
    """
    t_p_in, H_p_in = interp_peak(t_in, H_in)
    t_p_out, _ = interp_peak(Times, Q_out_reach)
    H_p_out = float(np.max(H_out))
    return {'t_p_in': t_p_in, 't_p_out': t_p_out, 'H_p_in': H_p_in, 'H_p_out': H_p_out,
            'delta_tp': t_p_out - t_p_in, 'delta_H': H_p_in - H_p_out, 'H_ratio': H_p_out / H_p_in}

def return_time(Times: np.ndarray, q_exc: np.ndarray) -> float:
    """
    Return time t_r (manuscript Eq. 19): first time after the positive exchange peak at which the net
    river-aquifer exchange changes sign (river->aquifer to aquifer->river), linearly interpolated.
    Measured from the start of the transient period. NaN if the sign never changes.
    """
    i_peak = int(np.argmax(q_exc))
    neg = np.flatnonzero(q_exc[i_peak:] < 0)
    if len(neg) == 0:
        return np.nan
    k = i_peak + neg[0]
    return float(Times[k-1] + (Times[k] - Times[k-1]) * q_exc[k-1] / (q_exc[k-1] - q_exc[k]))

def dir_size_mb(path: str) -> float:
    return sum(f.stat().st_size for f in Path(path).rglob('*') if f.is_file()) / 1024**2

def run_single_simulation(simulation_params: Tuple) -> Dict:
    """
    Run a single MODFLOW simulation with isolated directory and error handling
    
    Parameters:
    -----------
    simulation_params : tuple
        (simulation_id, parameter_row) where parameter_row contains all simulation parameters
        
    Returns:
    --------
    dict : Results dictionary with simulation outcomes and metadata
    """
    simulation_id, params, base_results_dir, simulation_subdir, modflow_exe_path = simulation_params[:5]
    time_cfg = {**DEFAULT_TIME_CFG, **(simulation_params[5] if len(simulation_params) > 5 else {})}

    start_time = time.time()
    logger.info(f"Starting simulation {simulation_id}")
    
    try:
        # Create isolated directory for this simulation
        path_results = create_simulation_directory(simulation_id, base_results_dir, simulation_subdir)
        name_model = f'water_sim_bagre_{simulation_id}'
        
        # Extract parameters from the row
        d_bkf     = params["d_bkf"]
        w_bkf     = params["w_bkf"]
        L         = params["L"]
        B         = params["B"]
        D         = params["D"]
        Ka        = params["Ka"]
        Kr        = params["Kr"]
        dr        = params["dr"]
        Sy        = params["Sy"]
        Slope     = params["S"]
        Ss        = 1e-7
        
        n_manning = params["n"]
        H_0       = params["H_0"]
        t_d_days  = params["t_d_days"]
        
        deltay = Slope*L
        Slope_trans = 1e-20
        Kv_Kh = 1
        
        # Hydrograph parameters
        H_0 = H_0
        eta_0 = (d_bkf-H_0)/H_0
        t_d = t_d_days*24*3600
        eta_d = 0.5
        H_p = eta_0 * H_0
        
        nrow = int(params["nrow"])
        ncol = 65
        delc = L/nrow
        delr = B/ncol
        nlay = 4
        
        # Create MODFLOW model with isolated workspace
        model = flopy.modflow.Modflow(
            name_model,
            exe_name=modflow_exe_path,
            version="mfnwt",
            model_ws=path_results,
            verbose=False  # Reduce output for parallel processing
        )
        
        hydrom_times = 600
        t = np.arange(0, t_d, hydrom_times)
        T = Ka * D
        
        # Generate hydrograph
        H_b = hydrograph_cooper_rorabaugh(t, H_0, eta_0, eta_d, t_d)
        Q = (1/n_manning) * (H_b**(5/3)) * (Slope**(1/2)) * (w_bkf)
        idx_Qmax = np.argmax(Q)
        
        velocity = Q / (w_bkf * d_bkf)
        T_arrive = (L/velocity)/(3600*24)
        
        v_mean = Q.mean()/(w_bkf*d_bkf)
        v_max = Q.max()/(w_bkf*d_bkf)
        
        t_mean = L/v_mean 
        t_max = L/v_max
        t_channel = L / v_max
        
        V_hid = (np.trapezoid(Q, t))
        tc = dn.get_tc(Sy=Sy, Ka=Ka, B=B, d_bkf=d_bkf)
        
        cell_size = np.sqrt(delr * delc)
        tc3 = Sy * cell_size**2 / (4 * T)
        t_peak = t[idx_Qmax]
        rates = 1/np.array([t_channel, tc, t_peak, tc3])
        
        TC = len(rates) / np.sum(rates)
        delta_t_TC = 0.25 * TC
        reach_len = L / nrow
        area = params["da_km2"]
        dakm2 = params["da_km2"]

        if dakm2 == 1000:
            cteCourant = .5
        elif dakm2 == 10000:
            cteCourant = .35
        else:
            cteCourant = .3
        
        duration_days = time_cfg['duration_days'] or duration_by_basin(dakm2)
        transient_length = duration_days * 24 * 3600
        t_delay = 1*24*3600
        numtim = 50

        if time_cfg['mode'] == 'courant':
            # Old scheme: one uniform transient step limited by Courant (flow velocity) and TC
            delta_t_courant = time_cfg['courant_factor']*cteCourant*reach_len / v_max
            delta_t = min(delta_t_TC, delta_t_courant)
            trans_perlen = [transient_length]
            trans_nstp = [int(time_cfg['nstp_uniform'] or np.ceil(transient_length / delta_t))]
            trans_tsmult = [1.0]
        else:
            trans_perlen, trans_nstp, trans_tsmult = build_adaptive_periods(
                transient_length, t_delay, t_d, v_max, reach_len, numtim, time_cfg)

        perlen = [1*24*3600] + trans_perlen
        nstp = [1] + trans_nstp
        tsmult = [1.0] + trans_tsmult

        # Generate topography
        topo = generate_topo(nrow, ncol, L=L, delta_h=deltay, slope_t=Slope_trans, min_topo=D+d_bkf+1, plot=False)
        botm = np.linspace(D-2, -2, nlay)

        nper = len(perlen)
        steady = [True] + [False]*(nper - 1)
        
        # Create discretization
        dis = flopy.modflow.ModflowDis(
            model,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            top=topo,
            botm=botm,
            itmuni=1,
            lenuni=2,
            nper=nper,
            perlen=perlen,
            nstp=nstp,
            tsmult=tsmult,
            steady=steady,
        )

        # Hydrograph input
        t0 = perlen[0]
        hydrogram = np.c_[t + t0 + t_delay, np.vstack(Q)]
        # Base flow after the event until past the end of the run: SFR averages the table over each step and
        # takes zero flow beyond the last entry, which produced a one-step drop of the inflow at t_d
        hydrogram = np.vstack([hydrogram, [t0 + t_delay + t_d, Q[0]], [t0 + transient_length + 86400, Q[0]]])
        np.savetxt(os.path.join(path_results, "QHidrom.dat"), hydrogram, fmt="%.4f")
        
        model.add_external(
            fname='QHidrom.dat',
            unit=70,
        )
        
        tabfiles_dict = {
            1: {
                'numval': len(hydrogram),
                'inuit': 70
            }
        }
        
        # River definition
        river_params_basic = {
            'A': 0,
            'B': 0.5,
            'C': 0,
            'start_cell': int(ncol/2),
            'orientation': "v",
            'width_riv': 1,
            'plot': False
        }
        
        stage_per_time = [.5]
        orientation = river_params_basic["orientation"]
        
        all_river = generate_river_stress_period_data(
            nrow=nrow, ncol=ncol, topo=topo, nper=nper, stage_per_time=stage_per_time,
            del_rbot=5, river_params=river_params_basic, type_gen='basic'
        )
        
        i_riv = all_river["i"]
        j_riv = all_river["j"]
        river = all_river["river"]
        nc_riv = all_river["nc_river"]
        
        # UPW package
        upw = flopy.modflow.ModflowUpw(
            model,
            hk=nlay*[Ka],
            vka=nlay*[Kv_Kh*Ka],
            sy=nlay*[Sy],
            ss=nlay*[Ss],
            laytyp=nlay*[1],
            hdry=-1e30,
            ipakcb=53 if time_cfg['save_cbc'] else None,
            iphdry=0
        )
        
        # NWT solver
        nwt = flopy.modflow.ModflowNwt(model, iprnwt=1, maxiterout=20000, fluxtol=1000, headtol=1e-5)
        
        # SFR reach data
        reach_data = flopy.modflow.ModflowSfr2.get_empty_reach_data(nreaches=len(i_riv))
        reach_data["k"] = 0
        reach_data["i"] = i_riv
        reach_data["j"] = j_riv
        reach_data["iseg"] = 1
        reach_data["rchlen"] = delc
        reach_data["strtop"] = topo[i_riv, j_riv]-d_bkf
        reach_data["ireach"] = list(range(1, nc_riv+1))
        reach_data["slope"] = 0
        reach_data["strthick"] = 0
        reach_data["strhc1"] = 0
        
        # SFR segment data
        segment_data_steady = flopy.modflow.ModflowSfr2.get_empty_segment_data(1)
        segment_data_steady[0]["nseg"] = 1
        segment_data_steady[0]["icalc"] = 0
        segment_data_steady[0]["outseg"] = 0
        segment_data_steady[0]["flow"] = 0
        segment_data_steady[0]["elevup"] = topo[i_riv, j_riv][0] - d_bkf
        segment_data_steady[0]["elevdn"] = topo[i_riv, j_riv][-1] - d_bkf
        segment_data_steady[0]["hcond1"] = Kr
        segment_data_steady[0]["hcond2"] = Kr
        segment_data_steady[0]["thickm1"] = dr
        segment_data_steady[0]["thickm2"] = dr
        segment_data_steady[0]["roughch"] = n_manning
        segment_data_steady[0]["roughbk"] = 0
        segment_data_steady[0]["width1"] = w_bkf
        segment_data_steady[0]["width2"] = w_bkf
        segment_data_steady[0]["depth1"] = H_0
        segment_data_steady[0]["depth2"] = H_0
        
        segment_data_transient = flopy.modflow.ModflowSfr2.get_empty_segment_data(1)
        segment_data_transient[0]["nseg"] = 1
        segment_data_transient[0]["icalc"] = 1
        segment_data_transient[0]["outseg"] = 0
        segment_data_transient[0]["flow"] = Q.min()
        segment_data_transient[0]["elevup"] = topo[i_riv, j_riv][0] - d_bkf
        segment_data_transient[0]["elevdn"] = topo[i_riv, j_riv][-1] - d_bkf
        segment_data_transient[0]["hcond1"] = Kr
        segment_data_transient[0]["hcond2"] = Kr
        segment_data_transient[0]["thickm1"] = dr
        segment_data_transient[0]["thickm2"] = dr
        segment_data_transient[0]["roughch"] = n_manning
        segment_data_transient[0]["roughbk"] = 0
        segment_data_transient[0]["width1"] = w_bkf
        segment_data_transient[0]["width2"] = w_bkf
        
        # SFR package
        sfr = flopy.modflow.ModflowSfr2(
            model,
            nstrm=-reach_data["ireach"].size,
            nsfrpar=0,
            nss=segment_data_steady["nseg"].size,
            const=1.0,
            dleak=1e-5,
            isfropt=0,
            ipakcb=53 if time_cfg['save_cbc'] else None,
            istcb2=54,
            irtflg=3,
            numtim=numtim,
            weight=.75,
            flwtol=1e-5,
            reach_data=reach_data,
            segment_data={0: segment_data_steady, 1: segment_data_transient},
            transroute=True,
            tabfiles=True,
            tabfiles_dict=tabfiles_dict,
        )

        
        # Boundary conditions
        ibound = np.ones((nlay, nrow, ncol), dtype=int)
        ibound[:, 0, :] = 1
        ibound[:, -1, :] = 1
        ibound[:, :, 0] = 0
        ibound[:, :, -1] = 0
        ibound[-1] = 0

        strt = np.ones((nlay, nrow, ncol), dtype=float)*(topo - d_bkf + H_0)
        ghb_data = []

        connection_thickness = D

        for k in range(nlay):
            for j in range(ncol):
                i_top = 0
                topo_top = topo[i_top, j]
                h_top = topo_top - d_bkf + H_0
                cond_top = time_cfg['ghb_cond_factor'] * Ka * w_bkf * delc / connection_thickness
                ghb_data.append([k, i_top, j, h_top, cond_top])
        
                i_bot = nrow - 1
                topo_bot = topo[i_bot, j]
                h_bot = topo_top - d_bkf + H_0 - Slope * L
                cond_bot = time_cfg['ghb_cond_factor'] * Ka * w_bkf * delc / connection_thickness
                ghb_data.append([k, i_bot, j, h_bot, cond_bot])
        
        ghb = flopy.modflow.ModflowGhb(
            model,
            stress_period_data={0: ghb_data, 1: ghb_data},
            ipakcb=53 if time_cfg['save_cbc'] else None
        )
        
        bas = flopy.modflow.ModflowBas(model, ibound=ibound, strt=strt)
        
        # Output control - Essential for SFR output generation
        stress_period_data = {}
        for per in range(len(perlen)):
            for step in range(nstp[per]):
                # 'save budget' is required for the .sfr.out; heads (.hds) only if requested (large files)
                stress_period_data[(per, step)] = (['save head'] if time_cfg['save_heads'] else []) + ['save budget'] \
                    + (['print budget'] if time_cfg['print_budget'] else [])
        
        oc = flopy.modflow.ModflowOc(
            model,
            stress_period_data=stress_period_data,
            compact=True
        )

        # Write and run model
        model.write_input()
        success, buff = model.run_model(silent=True, report=False)
        
        if not success:
            raise RuntimeError(f"MODFLOW simulation {simulation_id} failed")
        
        # Extract results
        # End time of each transient step, from the start of the transient (non-uniform with 'adaptive')
        Times = step_end_times(perlen[1:], nstp[1:], tsmult[1:])

        sfr_data = extract_sfr_data(model, [-1], ["Qaquifer"])
        aggregate_data = np.vstack([sfr_data[n]["Qaquifer"] for n in range(nrow)]).sum(axis=0)[1:]
        results_Qexc = analyze_peak_integrals(Times, aggregate_data)
        t_r = return_time(Times, aggregate_data)

        # Flood-wave metrics (manuscript Eq. 16-17): prescribed inflow hydrograph vs simulated outlet reach
        sfr_reach = extract_sfr_data(model, [1, nrow-2], ["depth", "Qin"])
        H_in = sfr_reach[1]["depth"][1:]
        H_out = sfr_reach[nrow-2]["depth"][1:]
        wave = flood_wave_metrics(t + t_delay, H_b, Times, sfr_reach[nrow-2]["Qin"][1:], H_out)
        # Previous definition (list_processing/visualization_30082025.ipynb): argmax of the simulated depth
        # at reach 2 and reach nrow-1. Kept only to compare with the Excel master.
        idx_0, idx_f = np.argmax(H_in), np.argmax(H_out)

        # Cumulative volume budget at the end of the transient period (as in list_processing/list_file.ipynb)
        mfl = flopy.utils.MfListBudget(os.path.join(path_results, f"{name_model}.list"))
        budget = mfl.get_dataframes()[1].iloc[-1].to_dict()

        # Return results dictionary
        results = {
            'simulation_id': simulation_id,
            'success': True,
            'T_in': results_Qexc["t_max"],
            'T_out': results_Qexc["t_min"],
            'Qexc_in': results_Qexc["int_pos"],
            'Qexc_out': -results_Qexc["int_neg"],
            'V_hid': V_hid,
            **wave,
            't_r': t_r,
            'delta_tp_reach': Times[idx_f] - Times[idx_0],
            'H_ratio_reach': H_out[idx_f] / H_in[idx_0],
            'delta_H_reach': abs(H_out[idx_f] - H_in[idx_0]),
            **budget,
            # Time discretization diagnostics
            'time_mode': time_cfg['mode'],
            'duration_days': duration_days,
            'nstp_total': int(sum(nstp[1:])),
            'dt_min_s': float(np.diff(np.r_[0, Times]).min()),
            'return_closed': bool(results_Qexc["bounds_neg"][1] < len(aggregate_data) - 1),
            'Qexc_end_ratio': float(abs(aggregate_data[-1]) / np.max(aggregate_data)),
            'disk_mb': dir_size_mb(path_results),
            'execution_time': time.time() - start_time,
            'path_results': path_results
        }
        
        # Clean up large files after successful simulation
        cleanup_simulation_files(path_results, keep_essential=True)
        
        logger.info(f"Simulation {simulation_id} completed successfully in {results['execution_time']:.2f}s")
        return results
        
    except Exception as e:
        error_msg = f"Simulation {simulation_id} failed: {str(e)}"
        logger.error(error_msg)
        
        # Clean up failed simulation directory to save space
        if 'path_results' in locals():
            cleanup_simulation_files(path_results, keep_essential=False)
            
        return {
            'simulation_id': simulation_id,
            'success': False,
            'error': str(e),
            'execution_time': time.time() - start_time
        }

def consolidate_results(results_list: List[Dict], original_dt: pd.DataFrame) -> pd.DataFrame:
    """
    Consolidate simulation results back into the original DataFrame
    
    Parameters:
    -----------
    results_list : list
        List of result dictionaries from parallel simulations
    original_dt : pd.DataFrame
        Original parameter DataFrame
        
    Returns:
    --------
    pd.DataFrame : Updated DataFrame with results
    """
    logger.info("Consolidating results...")
    
    # Create a copy of the original DataFrame
    dt_results = original_dt.copy()
    
    # Initialize result columns if they don't exist
    result_columns = ['T_in', 'T_out', 'Qexc_in', 'Qexc_out', 'V_hid', 'success', 'execution_time']
    for col in result_columns:
        if col not in dt_results.columns:
            dt_results[col] = np.nan
    
    # simulation_id is IDMERGE when the input has it, else the row position
    ids = dt_results["IDMERGE"].tolist() if "IDMERGE" in dt_results.columns else list(range(len(dt_results)))
    row_of = {sid: k for k, sid in enumerate(ids)}

    # Update results for each simulation
    for result in results_list:
        sim_id = row_of[result['simulation_id']]
        if result['success']:
            # Copy every exported metric (exchange, peak travel, list budget)
            for key, value in result.items():
                if key not in ('simulation_id', 'path_results'):
                    dt_results.loc[sim_id, key] = value
        else:
            dt_results.loc[sim_id, 'success'] = False
            logger.warning(f"Simulation {sim_id} failed: {result.get('error', 'Unknown error')}")
        
        dt_results.loc[sim_id, 'execution_time'] = result['execution_time']
    
    return dt_results

def run_parallel_simulations(n_simulations: int = None, max_workers: int = None, cleanup_mode: str = 'essential',
                             input_csv_path: str = str(INPUTS),
                             base_results_dir: str = str(DATA_OUTPUT),
                             simulation_subdir: str = "simple_model",
                             modflow_exe_path: str = DEFAULT_EXE,
                             time_cfg: Dict = None):
    """
    Main function to run parallel MODFLOW simulations with space management
    
    Parameters:
    -----------
    n_simulations : int, optional
        Number of simulations to run. If None, uses all rows in sample_points.csv
    max_workers : int, optional
        Maximum number of parallel workers. If None, uses min(cpu_count, 4)
    cleanup_mode : str, optional
        Cleanup strategy: 'essential' (keep minimal files), 'none' (keep all), 'full' (remove all)
    """
    
    # Load parameters
    logger.info("Loading simulation parameters...")
    dt = pd.read_csv(input_csv_path, sep=None, engine="python", encoding="utf-8-sig")  # accepts ',' or ';'
    if "run" in dt.columns:                      # run = 0 skips a row
        dt = dt[dt["run"] == 1]
    if "ensemble" in dt.columns and dt["ensemble"].nunique() > 1:
        raise ValueError("The input mixes ensembles, whose IDMERGE values repeat: run one ensemble at a time "
                         "(load_inputs() or simulations/run_campaign.py).")
    # simulation_id is positional (dt.iloc[i]); keep index aligned so results land on the right row
    dt = dt.reset_index(drop=True)

    # Add result columns
    dt["T_in"] = 0.0
    dt["T_out"] = 0.0
    dt["Qexc_in"] = 0.0
    dt["Qexc_out"] = 0.0
    dt["V_hid"] = 0.0
    
    # Determine simulation range
    if n_simulations is None:
        n_simulations = len(dt)
    else:
        n_simulations = min(n_simulations, len(dt))
    
    # Determine number of workers (reduced default for space management)
    if max_workers is None:
        max_workers = min(mp.cpu_count(), 4)  # Conservative default for space management
    
    logger.info(f"Starting {n_simulations} simulations with {max_workers} workers")
    logger.info(f"Cleanup mode: {cleanup_mode}")
    
    # Estimate space requirements
    estimated_space_per_sim = 30  # MB per simulation (before cleanup)
    max_concurrent_space = max_workers * estimated_space_per_sim
    logger.info(f"Estimated max concurrent space usage: {max_concurrent_space} MB")
    
    # Prepare simulation parameters
    # Run folders are named by IDMERGE (unique, stable when a half-finished campaign is resumed)
    ids = dt["IDMERGE"].astype(int).tolist() if "IDMERGE" in dt.columns else list(range(len(dt)))
    simulation_params = [(ids[i], dt.iloc[i], base_results_dir, simulation_subdir, modflow_exe_path, time_cfg or {})
                         for i in range(n_simulations)]
    
    # Create base results directory and subdirectory
    base_path = Path(base_results_dir)
    base_path.mkdir(parents=True, exist_ok=True)

    full_results_path = base_path / simulation_subdir
    full_results_path.mkdir(parents=True, exist_ok=True)
    
    start_total = time.time()
    results_list = []
    
    # Run parallel simulations
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all simulations
        future_to_sim = {
            executor.submit(run_single_simulation, params): params[0] 
            for params in simulation_params
        }
        
        # Collect results as they complete
        completed = 0
        for future in as_completed(future_to_sim):
            sim_id = future_to_sim[future]
            try:
                result = future.result()
                results_list.append(result)
                completed += 1
                # One JSON line per finished row: progress tracking and nothing lost if the run stops
                with open(full_results_path / "rows_done.jsonl", "a") as f:
                    f.write(json.dumps({"IDMERGE": sim_id, **result}, default=str) + "\n")
                
                with progress_lock:
                    progress = (completed / n_simulations) * 100
                    logger.info(f"Progress: {completed}/{n_simulations} ({progress:.1f}%)")
                    
            except Exception as e:
                logger.error(f"Simulation {sim_id} generated an exception: {e}")
                results_list.append({
                    'simulation_id': sim_id,
                    'success': False,
                    'error': str(e),
                    'execution_time': 0.0
                })
    
    total_time = time.time() - start_total
    
    # Consolidate and save results
    logger.info("Consolidating and saving results...")
    dt_final = consolidate_results(results_list, dt)
    
    # Save results with timestamp
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    output_filename = str(full_results_path / f"complete_sample_points_parallel_{timestamp}.csv")
    dt_final.to_csv(output_filename, index=False)
    
    # Generate summary report
    successful_sims = sum(1 for r in results_list if r['success'])
    failed_sims = n_simulations - successful_sims
    avg_time = np.mean([r['execution_time'] for r in results_list if r['success']])
    
    logger.info("="*60)
    logger.info("PARALLEL SIMULATION SUMMARY")
    logger.info("="*60)
    logger.info(f"Total simulations: {n_simulations}")
    logger.info(f"Successful: {successful_sims}")
    logger.info(f"Failed: {failed_sims}")
    logger.info(f"Success rate: {(successful_sims/n_simulations)*100:.1f}%")
    logger.info(f"Total execution time: {total_time:.2f}s")
    logger.info(f"Average time per simulation: {avg_time:.2f}s")
    logger.info(f"Speedup factor: {avg_time*n_simulations/total_time:.1f}x")
    logger.info(f"Results saved to: {output_filename}")
    logger.info("="*60)
    
    return dt_final, results_list

if __name__ == "__main__":
    # Configuration - SPACE EFFICIENT SETTINGS
    N_SIMULATIONS = None      # Set to None to run all simulations in CSV
    MAX_WORKERS = 10         # Reduced from 16 to 4 for space management
    CLEANUP_MODE = 'none'  # 'essential', 'none', or 'full'

    # Ensemble to run (rows of data/input/model_inputs.csv); simulations/run_campaign.py does the same, resumable
    ENSEMBLE = "main"          # "main" or "high_connectivity"
    ORIGINAL_N = False         # True: earlier Manning n (n_original), for comparison

    # Path configuration
    BASE_RESULTS_DIR = str(DATA_OUTPUT)                                       # Base directory for all results
    SIMULATION_SUBDIR = ENSEMBLE + ("_n_original" if ORIGINAL_N else "")      # Subdirectory within results
    INPUT_CSV_PATH = Path(BASE_RESULTS_DIR) / SIMULATION_SUBDIR / "_inputs.csv"
    INPUT_CSV_PATH.parent.mkdir(parents=True, exist_ok=True)
    load_inputs(ENSEMBLE, ORIGINAL_N).to_csv(INPUT_CSV_PATH, index=False)
    MODFLOW_EXE_PATH = DEFAULT_EXE                              # Path to MODFLOW executable
    # Time-step scheme validated in test_time_steps.py / test_mass_balance.py ('fine'); no .cbc files
    TIME_CFG = {'n_hyd': 400, 'tsmult': 1.05, 'dt_max_hours': 3.0, 'save_cbc': False}
    
    logger.info("="*60)
    logger.info("SPACE-EFFICIENT PARALLEL MODFLOW SIMULATION")
    logger.info("="*60)
    logger.info(f"Configuration:")
    logger.info(f"  - Max workers: {MAX_WORKERS} (reduced for space efficiency)")
    logger.info(f"  - Cleanup mode: {CLEANUP_MODE}")
    logger.info(f"  - Auto-cleanup: Enabled")
    logger.info("="*60)
    
    # Run parallel simulations
    try:
        results_df, results_list = run_parallel_simulations(
            n_simulations=N_SIMULATIONS,
            max_workers=MAX_WORKERS,
            cleanup_mode=CLEANUP_MODE,
            input_csv_path=str(INPUT_CSV_PATH),
            base_results_dir=BASE_RESULTS_DIR,
            simulation_subdir=SIMULATION_SUBDIR,
            modflow_exe_path=MODFLOW_EXE_PATH,
            time_cfg=TIME_CFG
        )
        logger.info("All simulations completed successfully!")
        
        # Final cleanup summary
        results_path = Path(f"{BASE_RESULTS_DIR}/{SIMULATION_SUBDIR}/")
        total_dirs = len([d for d in results_path.iterdir() if d.is_dir()]) if results_path.exists() else 0
        logger.info(f"Total simulation directories: {total_dirs}")
        logger.info("Space optimization: Large files automatically removed after each simulation")
        
    except KeyboardInterrupt:
        logger.warning("Simulation interrupted by user")
    except Exception as e:
        logger.error(f"Critical error: {e}")
        raise