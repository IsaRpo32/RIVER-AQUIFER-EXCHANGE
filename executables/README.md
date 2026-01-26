# MODFLOW Executables

This directory should contain the MODFLOW-NWT executable required to run the simulations.

## Obtaining MODFLOW-NWT

MODFLOW-NWT is distributed by the U.S. Geological Survey (USGS) and can be downloaded from:

**Official USGS Source:**
https://www.usgs.gov/software/modflow-nwt-newton-formulation-modflow-2005

### Download Instructions

1. Visit the USGS MODFLOW-NWT page
2. Download the appropriate version for your operating system:
   - Windows: `mfnwt.exe`
   - Linux: Compile from source or use pre-compiled binaries
   - macOS: Compile from source

3. Place the executable in this directory

### Expected Files

After download, this directory should contain:

```
executables/
├── README.md          # This file
├── mfnwt.exe          # Windows executable (or mfnwt for Linux/macOS)
└── mf2005.exe         # Optional: MODFLOW-2005 executable
```

## Alternative: Using FloPy's Executable Manager

FloPy can automatically download MODFLOW executables:

```python
import flopy

# Get the path to executables
exe_path = flopy.which("mfnwt")

# Or download using flopy's get-modflow utility
# python -m flopy.utils.get_modflow :mfnwt
```

## Verification

To verify the executable works correctly:

```python
import subprocess

# Test execution
result = subprocess.run(
    ["./executables/mfnwt.exe", "-v"],
    capture_output=True,
    text=True
)
print(result.stdout)
```

## Note on Version Compatibility

This project was developed using MODFLOW-NWT version 1.2.0. While newer versions should be compatible, please verify your results if using a different version.

## License

MODFLOW is public domain software developed by the U.S. Geological Survey.
