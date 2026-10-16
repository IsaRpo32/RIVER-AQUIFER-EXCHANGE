# MODFLOW-NWT executable

Only the simulations (`simulations/run_campaign.py`, `simulations/parallel_modflow.py`, `validation/`) need
MODFLOW-NWT. The post-processing, statistics and figures run from the results in `data/output/`.

## Expected files

The scripts look for the executable here, with these names:

```
executables/
├── MODFLOW-NWT.exe    Windows
└── mfnwt              Linux / macOS
```

(`DEFAULT_EXE` in `simulations/parallel_modflow.py`.)

## Version

The paper simulations used **MODFLOW-NWT 1.3.0** (USGS): the official Windows 64-bit executable for half B of
the main ensemble and an Intel-compiled Linux build for half A. Both give the same results within 0.00004 %.

Download the official release from the USGS:
https://www.usgs.gov/software/modflow-nwt-newton-formulation-modflow-2005

On Linux and macOS compile it from the source in that release, or install it with FloPy:

```bash
get-modflow :flopy        # downloads the USGS executables, mfnwt included, into FloPy's bin folder
```

and copy `mfnwt` into this folder.

## Note on compilers

Different builds of the same version (e.g. Intel vs gfortran) differ by round-off: heads by up to a few mm
and cumulative exchange by about 0.1 %. Use a single executable for every run of a set that will be compared
or pooled.

## License

MODFLOW is public domain software developed by the U.S. Geological Survey.
