"""
Parallel MODFLOW Simulation for River-Aquifer Exchange

This script runs parallel MODFLOW-NWT simulations to model river-aquifer
exchange during flood events. It uses the Streamflow Routing (SFR2) package
to simulate dynamic surface water-groundwater interaction.

Usage:
    python parallel_modflow.py

Requirements:
    - MODFLOW-NWT executable in ../executables/
    - Input CSV file with simulation parameters
"""

import os
import sys

# Add parent directory to path for src imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import glob
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
from flopy.utils.sfroutputfile import SfrFile

from src.modflow_utils import (
    generate_topo,
    generate_river_stress_period_data,
    hydrograph_cooper_rorabaugh,
    extract_sfr_data,
    analyze_peak_integrals,
)
from src import dimensionless_numbers as dn


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


def create_simulation_directory(simulation_id: int) -> str:
    """Create unique directory for each simulation."""
    path_results = f"../data/output/run_{simulation_id:03d}/"
    Path(path_results).mkdir(parents=True, exist_ok=True)
    return path_results


def cleanup_simulation_files(sim_dir: str, keep_essential: bool = True):
    """
    Clean up simulation files to save disk space.

    Parameters
    ----------
    sim_dir : str
        Simulation directory path.
    keep_essential : bool
        If True, keeps only essential files (QHidrom.dat, log files).
        If False, removes entire directory.
    """
    if not keep_essential:
        # Remove entire directory
        if os.path.exists(sim_dir):
            shutil.rmtree(sim_dir)
            logger.debug(f"Removed directory: {sim_dir}")
        return

    # Keep only essential files
    essential_files = [
        'QHidrom.dat',          # Hydrograph data (essential)
        '*.log',                # Log files
        '*.lst'                 # List files (small, useful for debugging)
    ]

    # Files to remove (large files)
    files_to_remove = [
        '*.hds',                # Head files (large)
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


def run_single_simulation(simulation_params: Tuple[int, pd.Series]) -> Dict:
    """
    Run a single MODFLOW simulation with isolated directory and error handling.

    Parameters
    ----------
    simulation_params : tuple
        (simulation_id, parameter_row) where parameter_row contains all simulation parameters.

    Returns
    -------
    dict
        Results dictionary with simulation outcomes and metadata.
    """
    simulation_id, params = simulation_params

    start_time = time.time()
    logger.info(f"Starting simulation {simulation_id}")

    try:
        # Create isolated directory for this simulation
        path_results = create_simulation_directory(simulation_id)
        name_model = f'water_sim_{simulation_id}'

        # Extract parameters from the row
        d_bkf = params["d_bkf"]
        w_bkf = params["w_bkf"]
        L = params["L"]
        B = params["B"]
        D = params["D"]
        Ka = params["Ka"]
        Kr = params["Kr"]
        dr = params["dr"]
        Sy = params["Sy"]
        Slope = params["S"]
        Ss = 1e-7

        n_manning = params["n"]
        H_0 = params["H_0"]
        t_d_days = params["t_d_days"]

        deltay = Slope * L
        Slope_trans = 1e-20
        Kv_Kh = 1

        # Hydrograph parameters
        eta_0 = (d_bkf - H_0) / H_0
        t_d = t_d_days * 24 * 3600
        eta_d = 0.5  # NOTE: Review this value
        H_p = eta_0 * H_0

        nrow = int(params["nrow"])
        ncol = 65
        delc = L / nrow
        delr = B / ncol
        nlay = 4

        # Create MODFLOW model with isolated workspace
        model = flopy.modflow.Modflow(
            name_model,
            exe_name="../executables/MODFLOW-NWT.exe",
            version="mfnwt",
            model_ws=path_results,
            verbose=False  # Reduce output for parallel processing
        )

        hydrom_times = 600
        t = np.arange(0, t_d, hydrom_times)
        T = Ka * D

        # Generate hydrograph
        H_b = hydrograph_cooper_rorabaugh(t, H_0, eta_0, eta_d, t_d)
        Q = (1 / n_manning) * (H_b ** (5 / 3)) * (Slope ** (1 / 2)) * w_bkf
        idx_Qmax = np.argmax(Q)

        velocity = Q / (w_bkf * d_bkf)
        T_arrive = (L / velocity) / (3600 * 24)

        v_mean = Q.mean() / (w_bkf * d_bkf)
        v_max = Q.max() / (w_bkf * d_bkf)

        t_mean = L / v_mean
        t_max = L / v_max
        t_canal = L / v_max

        V_hid = np.trapezoid(Q, t)
        tc = dn.get_tc(Sy=Sy, Ka=Ka, B=B, d_bkf=d_bkf)

        a_celda = np.sqrt(delr * delc)
        tc3 = Sy * a_celda ** 2 / (4 * T)
        t_pico = t[idx_Qmax]
        valores = 1 / np.array([t_canal, tc, t_pico, tc3])

        TC = len(valores) / np.sum(valores)
        delta_t_TC = 0.25 * TC
        reach_len = L / nrow
        dakm2 = params["da_km2"]

        if dakm2 == 1000:
            cteCourant = 0.55
        elif dakm2 == 10000:
            cteCourant = 0.5
        else:
            cteCourant = 0.5

        delta_t_courant = cteCourant * reach_len / v_max
        delta_t = min(delta_t_TC, delta_t_courant)
        courant_final = v_max * delta_t / reach_len

        if dakm2 == 100000:
            t_trans = 150
        else:
            t_trans = 90
        duracion_transitorio = t_trans * 24 * 3600

        nstp_trans = int(np.ceil(duracion_transitorio / delta_t))
        perlen = [1 * 24 * 3600, duracion_transitorio]
        nstp = [1, nstp_trans]

        # Generate topography
        topo = generate_topo(nrow, ncol, L=L, delta_h=deltay,
                            slope_t=Slope_trans, min_topo=D + d_bkf + 1, plot=False)
        botm = np.linspace(D - 2, -2, nlay)

        nper = len(perlen)
        steady = [True, False]

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
            steady=steady,
        )

        # Hydrograph input
        t0 = perlen[0]
        t_delay = 1 * 24 * 3600
        hydrogram = np.c_[t + t0 + t_delay, np.vstack(Q)]
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
            'start_cell': int(ncol / 2),
            'orientation': "v",
            'width_riv': 1,
            'plot': False
        }

        stage_per_time = [0.5]
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
            hk=nlay * [Ka],
            vka=nlay * [Kv_Kh * Ka],
            sy=nlay * [Sy],
            ss=nlay * [Ss],
            laytyp=nlay * [1],
            hdry=-1e30,
            ipakcb=53,
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
        reach_data["strtop"] = topo[i_riv, j_riv] - d_bkf
        reach_data["ireach"] = list(range(1, nc_riv + 1))
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
            ipakcb=53,
            istcb2=54,
            irtflg=3,
            numtim=50,
            weight=0.75,
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

        strt = np.ones((nlay, nrow, ncol), dtype=float) * (topo - d_bkf + H_0)
        ghb_data = []

        connection_thickness = D

        for k in range(nlay):
            for j in range(ncol):
                i_top = 0
                topo_top = topo[i_top, j]
                h_top = topo_top - d_bkf + H_0
                cond_top = Ka * w_bkf * delc / connection_thickness
                ghb_data.append([k, i_top, j, h_top, cond_top])

                i_bot = nrow - 1
                topo_bot = topo[i_bot, j]
                h_bot = topo_top - d_bkf + H_0 - Slope * L
                cond_bot = Ka * w_bkf * delc / connection_thickness
                ghb_data.append([k, i_bot, j, h_bot, cond_bot])

        ghb = flopy.modflow.ModflowGhb(
            model,
            stress_period_data={0: ghb_data, 1: ghb_data},
            ipakcb=53
        )

        bas = flopy.modflow.ModflowBas(model, ibound=ibound, strt=strt)

        # Output control
        stress_period_data = {}
        for per in range(len(perlen)):
            for step in range(nstp[per]):
                stress_period_data[(per, step)] = [
                    'save head',
                    'save budget'
                ]

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
        num_steps = nstp[1]
        steps = np.linspace(0, num_steps, num_steps)
        delta_t = perlen[1] / num_steps
        Times = steps * delta_t

        sfr_data = extract_sfr_data(model, [-1], ["Qaquifer"])
        aggregate_data = np.vstack([sfr_data[n]["Qaquifer"] for n in range(nrow)]).sum(axis=0)[1:]
        results_Qexc = analyze_peak_integrals(Times, aggregate_data)

        # Return results dictionary
        results = {
            'simulation_id': simulation_id,
            'success': True,
            'T_in': results_Qexc["t_max"],
            'T_out': results_Qexc["t_min"],
            'Qexc_in': results_Qexc["int_pos"],
            'Qexc_out': -results_Qexc["int_neg"],
            'V_hid': V_hid,
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
    Consolidate simulation results back into the original DataFrame.

    Parameters
    ----------
    results_list : list
        List of result dictionaries from parallel simulations.
    original_dt : pd.DataFrame
        Original parameter DataFrame.

    Returns
    -------
    pd.DataFrame
        Updated DataFrame with results.
    """
    logger.info("Consolidating results...")

    # Create a copy of the original DataFrame
    dt_results = original_dt.copy()

    # Initialize result columns if they don't exist
    result_columns = ['T_in', 'T_out', 'Qexc_in', 'Qexc_out', 'V_hid', 'success', 'execution_time']
    for col in result_columns:
        if col not in dt_results.columns:
            dt_results[col] = np.nan

    # Update results for each simulation
    for result in results_list:
        sim_id = result['simulation_id']
        if result['success']:
            dt_results.loc[sim_id, 'T_in'] = result['T_in']
            dt_results.loc[sim_id, 'T_out'] = result['T_out']
            dt_results.loc[sim_id, 'Qexc_in'] = result['Qexc_in']
            dt_results.loc[sim_id, 'Qexc_out'] = result['Qexc_out']
            dt_results.loc[sim_id, 'V_hid'] = result['V_hid']
            dt_results.loc[sim_id, 'success'] = True
        else:
            dt_results.loc[sim_id, 'success'] = False
            logger.warning(f"Simulation {sim_id} failed: {result.get('error', 'Unknown error')}")

        dt_results.loc[sim_id, 'execution_time'] = result['execution_time']

    return dt_results


def run_parallel_simulations(n_simulations: int = None, max_workers: int = None,
                             cleanup_mode: str = 'essential'):
    """
    Main function to run parallel MODFLOW simulations with space management.

    Parameters
    ----------
    n_simulations : int, optional
        Number of simulations to run. If None, uses all rows in input CSV.
    max_workers : int, optional
        Maximum number of parallel workers. If None, uses min(cpu_count, 4).
    cleanup_mode : str, optional
        Cleanup strategy: 'essential' (keep minimal files), 'none' (keep all), 'full' (remove all).
    """
    # Load parameters
    logger.info("Loading simulation parameters...")
    dt = pd.read_csv('../data/input/sample_points.csv')

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
    simulation_params = [(i, dt.iloc[i]) for i in range(n_simulations)]

    # Create results directory
    Path("../data/output/").mkdir(parents=True, exist_ok=True)

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
    output_filename = f"../data/output/results_{timestamp}.csv"
    dt_final.to_csv(output_filename, index=False)

    # Generate summary report
    successful_sims = sum(1 for r in results_list if r['success'])
    failed_sims = n_simulations - successful_sims
    avg_time = np.mean([r['execution_time'] for r in results_list if r['success']])

    logger.info("=" * 60)
    logger.info("PARALLEL SIMULATION SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total simulations: {n_simulations}")
    logger.info(f"Successful: {successful_sims}")
    logger.info(f"Failed: {failed_sims}")
    logger.info(f"Success rate: {(successful_sims / n_simulations) * 100:.1f}%")
    logger.info(f"Total execution time: {total_time:.2f}s")
    logger.info(f"Average time per simulation: {avg_time:.2f}s")
    logger.info(f"Speedup factor: {avg_time * n_simulations / total_time:.1f}x")
    logger.info(f"Results saved to: {output_filename}")
    logger.info("=" * 60)

    return dt_final, results_list


if __name__ == "__main__":
    # Configuration - SPACE EFFICIENT SETTINGS
    N_SIMULATIONS = None    # Set to None to run all simulations in CSV
    MAX_WORKERS = 5         # Reduced from 16 to 4 for space management
    CLEANUP_MODE = 'essential'  # 'essential', 'none', or 'full'

    logger.info("=" * 60)
    logger.info("SPACE-EFFICIENT PARALLEL MODFLOW SIMULATION")
    logger.info("=" * 60)
    logger.info("Configuration:")
    logger.info(f"  - Max workers: {MAX_WORKERS} (reduced for space efficiency)")
    logger.info(f"  - Cleanup mode: {CLEANUP_MODE}")
    logger.info(f"  - Auto-cleanup: Enabled")
    logger.info("=" * 60)

    # Run parallel simulations
    try:
        results_df, results_list = run_parallel_simulations(
            n_simulations=N_SIMULATIONS,
            max_workers=MAX_WORKERS,
            cleanup_mode=CLEANUP_MODE
        )
        logger.info("All simulations completed successfully!")

        # Final cleanup summary
        output_path = Path("../data/output/")
        total_dirs = len([d for d in output_path.iterdir() if d.is_dir()])
        logger.info(f"Total simulation directories: {total_dirs}")
        logger.info("Space optimization: Large files automatically removed after each simulation")

    except KeyboardInterrupt:
        logger.warning("Simulation interrupted by user")
    except Exception as e:
        logger.error(f"Critical error: {e}")
        raise
