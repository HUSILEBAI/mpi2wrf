# Climate Data Processing Scripts for MPI-ESM1-2-HR Model

This repository contains scripts for processing climate data from the MPI-ESM1-2-HR model. These scripts handle data interpolation, selection, and conversion to a common grid format suitable for further analysis.

## Table of Contents

- [Scripts](#scripts)
  - [Bash Script](#bash-script)
  - [Python Script](#python-script)
- [Dependencies](#dependencies)
- [Usage](#usage)
  - [Bash Script Usage](#bash-script-usage)
  - [Python Script Usage](#python-script-usage)
- [Author](#author)
- [Acknowledgement](#acknowledgement)

## Scripts

### Bash Script

The Bash script processes the Sea Surface Temperature (SST) data from the MPI-ESM1-2-HR historical dataset by selecting specific time steps and interpolating the data to a target grid.

**File**: `process_sst.sh`

#### Script Details

1. **Input Files**:
    ```bash
    input_files=(
        "tos_3hr_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_198001010300-198501010000.nc"
        "tos_3hr_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_198501010300-199001010000.nc"
        "tos_3hr_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_199001010300-199501010000.nc"
        "tos_3hr_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_199501010300-200001010000.nc"
        "tos_3hr_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_200001010300-200501010000.nc"
        "tos_3hr_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_200501010300-201001010000.nc"
        "tos_3hr_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_201001010300-201501010000.nc"
    )
    ```

2. **Target Grid Description**:
    The target grid is defined in `target_grid.txt`:
    ```bash
    gridtype = lonlat
    gridsize = 73728
    xsize = 384
    ysize = 192
    xfirst = -179.53125
    xinc = 0.9375
    yfirst = -89.4625
    yinc = 0.9375
    ```

3. **Processing Loop**:
    The script loops through each input file, adjusts the start time, selects specific time steps, performs bilinear interpolation, and saves the output file.

### Python Script

The Python script processes climate data from the MPI-ESM1-2-HR model by reading multiple variables, performing boundary checks, and ensuring data integrity by handling missing values. The script then interpolates the data to a common grid and saves the results for further analysis.

**File**: `process_mpi2wrf.py`

#### Script Details

1. **Constants and Directories**:
    - **Start and End Dates**:
        ```python
        START_DATE = DT.datetime(1980, 1, 1, 6, 0, 0)
        END_DATE = DT.datetime(1985, 1, 1, 0, 0, 0)
        ```
    - **Directories**:
        ```python
        INPUT_DIR = '/path/to/input/'
        OUTPUT_DIR = '/path/to/output/'
        ```

2. **Variable File Patterns**:
    The script handles multiple variables such as `tas`, `huss`, `ps`, etc.

3. **Utility Functions**:
    - `format_time()`: Formats datetime objects.
    - `read_and_process_variable()`: Reads and processes 2D variables from NetCDF files.
    - `read_and_process_3d_variable()`: Reads and processes 3D variables from NetCDF files.
    - `map_to_pywinter_name()`: Maps variable names to PyWinter naming conventions.

4. **Processing Functions**:
    - `process_chunks()`: Processes data in chunks.
    - `process_variables()`: Main function to process variables from NetCDF files.

5. **Main Execution Flow**:
    The script defines the main execution flow and processes data in chunks:
    ```python
    if __name__ == "__main__":
        main()
    ```

## Dependencies

Ensure the following dependencies are installed before running the scripts:

- Bash
- [CDO](https://code.mpimet.mpg.de/projects/cdo) (for the Bash script)
- Python 3.x
- numpy
- pywinter
- netCDF4
- scipy

## Usage

### Bash Script Usage

1. **Set Execute Permission**:
    ```bash
    chmod +x process_sst.sh
    ```

2. **Run the Script**:
    ```bash
    ./process_sst.sh
    ```

### Python Script Usage

1. **Set Directories**:
    Ensure the `INPUT_DIR` and `OUTPUT_DIR` paths are correctly set in the script.

2. **Run the Script**:
    ```bash
    python process_mpi2wrf.py
    ```

## Author

**Husile Bai**  
Email: husile.bai@utah.edu

## Acknowledgement

The author extends gratitude to **Prof. Alfonso Fernandez** (alfernandez@udec.cl) for recommending the 'pywinter' package and providing codes for data processing.
