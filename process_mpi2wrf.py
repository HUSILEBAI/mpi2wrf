#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 20:37:20 2024

@author: u1215181
"""

import numpy as np
import pywinter.winter as pyw
import netCDF4 as ncf
import datetime as DT
import os
import glob
from scipy.interpolate import griddata

# Constants
START_DATE = DT.datetime(1980, 1, 1, 6, 0, 0)
END_DATE = DT.datetime(1985, 1, 1, 0, 0, 0)
BASE_DATE = DT.datetime(1850, 1, 1)
MISSING_VALUE = 1.e+19  # This is lower than the missing value, so use ">" given an issue with soil data that does not work with the exact value
FIELDS_3D = [3, 6, 8, 11, 13, 15, 17, 19, 20, 22, 24, 26, 27]
INPUT_DIR = '/uufs/chpc.utah.edu/common/home/strong-group7/husile/gsl/cmip_data/MPI-ESM1-2-HR/historical/'
OUTPUT_DIR = '/uufs/chpc.utah.edu/common/home/strong-group7/husile/gsl/cmip_data/MPI-ESM1-2-HR/output_mpi2int/hist/'

# Variable file patterns
variable_files = {
    'tas': f'tas_6hrPlevPt_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_{START_DATE.strftime("%Y%m%d%H%M")}-{END_DATE.strftime("%Y%m%d%H%M")}.nc',
    'huss': f'huss_6hrPlevPt_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_{START_DATE.strftime("%Y%m%d%H%M")}-{END_DATE.strftime("%Y%m%d%H%M")}.nc',
    'ps': f'ps_6hrLev_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_{START_DATE.strftime("%Y%m%d%H%M")}-{END_DATE.strftime("%Y%m%d%H%M")}.nc',
    'psl': f'psl_6hrPlevPt_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_{START_DATE.strftime("%Y%m%d%H%M")}-{END_DATE.strftime("%Y%m%d%H%M")}.nc',
    'vas': f'vas_6hrPlevPt_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_{START_DATE.strftime("%Y%m%d%H%M")}-{END_DATE.strftime("%Y%m%d%H%M")}.nc',
    'uas': f'uas_6hrPlevPt_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_{START_DATE.strftime("%Y%m%d%H%M")}-{END_DATE.strftime("%Y%m%d%H%M")}.nc',
    'tos': f'tos_6hrPlevPt_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_{START_DATE.strftime("%Y%m%d%H%M")}-{END_DATE.strftime("%Y%m%d%H%M")}.nc',
    'hus': f'hus_6hrPlevPt_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_{START_DATE.strftime("%Y%m%d%H%M")}-{END_DATE.strftime("%Y%m%d%H%M")}.nc',
    'ta': f'ta_6hrPlevPt_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_{START_DATE.strftime("%Y%m%d%H%M")}-{END_DATE.strftime("%Y%m%d%H%M")}.nc',
    'va': f'va_6hrPlevPt_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_{START_DATE.strftime("%Y%m%d%H%M")}-{END_DATE.strftime("%Y%m%d%H%M")}.nc',
    'ua': f'ua_6hrPlevPt_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_{START_DATE.strftime("%Y%m%d%H%M")}-{END_DATE.strftime("%Y%m%d%H%M")}.nc',
    'zg': f'zg_6hrPlevPt_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_{START_DATE.strftime("%Y%m%d%H%M")}-{END_DATE.strftime("%Y%m%d%H%M")}.nc',
    'mrsos': f'mrsos_6hrPlevPt_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_{START_DATE.strftime("%Y%m%d%H%M")}-{END_DATE.strftime("%Y%m%d%H%M")}.nc',
    'tsl': f'tsl_6hrPlevPt_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_{START_DATE.strftime("%Y%m%d%H%M")}-{END_DATE.strftime("%Y%m%d%H%M")}.nc'
}

# Utility Functions
def format_time(date_time):
    return date_time.strftime('%Y%m%d%H%M')

def read_and_process_variable(file_name, var_name, start_timesteps, end_timesteps):
    ds = ncf.Dataset(file_name, 'r')
    var_data = ds.variables[var_name][start_timesteps:end_timesteps, :, :]
    ds.close()
    missing_indices = np.where(var_data > MISSING_VALUE)
    var_data[missing_indices] = np.nan
    return var_data

def read_and_process_3d_variable(file_name, var_name, start_timesteps, end_timesteps):
    ds = ncf.Dataset(file_name, 'r')
    var_data = ds.variables[var_name][start_timesteps:end_timesteps, :, :, :]
    plevs = ds.variables['plev'][:]
    ds.close()
    missing_indices = np.where(var_data > MISSING_VALUE)
    var_data[missing_indices] = np.nan
    return var_data, plevs

def map_to_pywinter_name(var_name):
    var_map = {
        'tas': 'TT',
        'huss': 'SPECHUMD',
        'ps': 'PSFC',
        'psl': 'PMSL',
        'vas': 'VV',
        'uas': 'UU',
        'tos': 'SKINTEMP',
        'hus': 'SPECHUMD',
        'ta': 'TT',
        'va': 'VV',
        'ua': 'UU',
        'zg': 'GHT'
    }
    if var_name in var_map:
        return var_map[var_name]
    else:
        raise ValueError(f"Unknown variable name: {var_name}")

def print_variable_shapes(variables):
    for var_name, var_data_list in variables.items():
        print(f"Shapes for {var_name}:")
        for i, data in enumerate(var_data_list):
            print(f"  - Instance {i}: {data.shape}")

# Processing Functions
def process_chunks(total_steps, chunk_size, operation):
    for start in range(0, total_steps, chunk_size):
        end = min(start + chunk_size, total_steps)
        print()
        print('#' * 50)
        print(f'start = {start} \t end = {end}')
        print('#' * 50)
        print()
        operation(start, end)

def process_variables(start_timesteps, end_timesteps, interpolated_tos=None):
    variables_2d = {}
    variables_3d = {}
    variables_sl = {}
    plevs_3d = {}

    lat, lon, dlat, dlon = None, None, None, None
    date_variable = None

    for var_name, file_pattern in variable_files.items():
        file_list = glob.glob(os.path.join(INPUT_DIR, file_pattern))

        for file_name in file_list:
            print(f"Processing {var_name} from {file_name}")

            if var_name in ['tas', 'huss', 'ps', 'psl', 'vas', 'uas', 'tos']:
                if var_name not in variables_2d:
                    variables_2d[var_name] = []

                var_data = read_and_process_variable(file_name, var_name, start_timesteps, end_timesteps)
                variables_2d[var_name].append(var_data)

                if var_name == 'tas':
                    ds = ncf.Dataset(file_name, 'r')
                    lat = ds.variables['lat'][:]
                    lon = ds.variables['lon'][:]
                    dlat = np.abs(lat[1] - lat[0])
                    dlon = np.abs(lon[1] - lon[0])
                    date_variable = ds.variables['time'][start_timesteps:end_timesteps]
                    date_variable.units = ds.variables['time'].getncattr('units')
                    ds.close()

            elif var_name in ['mrsos', 'tsl']:
                if var_name not in variables_sl:
                    variables_sl[var_name] = []

                var_data = read_and_process_variable(file_name, var_name, start_timesteps, end_timesteps)
                variables_sl[var_name].append(var_data)

            elif var_name in ['hus', 'ta', 'va', 'ua', 'zg']:
                if var_name not in variables_3d:
                    variables_3d[var_name] = []
                    plevs_3d[var_name] = None

                try:
                    var_data, plevs = read_and_process_3d_variable(file_name, var_name, start_timesteps, end_timesteps)
                    variables_3d[var_name].append(var_data)

                    if plevs_3d[var_name] is None:
                        plevs_3d[var_name] = plevs
                except Exception as e:
                    print(f"Error processing {var_name} from {file_name}: {e}")

    print(f"2D Variables: {variables_2d.keys()}")
    print(f"3D Variables: {variables_3d.keys()}")
    print_variable_shapes(variables_2d)
    print_variable_shapes(variables_3d)

    if lat is not None and lon is not None:
        winter_geo = pyw.Geo0(lat[0], lon[0], dlat, dlon)

        for n in range(end_timesteps - start_timesteps):
            total_fields = []

            for var_name, var_data_list in variables_2d.items():
                var_data = var_data_list[0][n, :, :]
                pywinter_var_name = map_to_pywinter_name(var_name)
                total_fields.append(pyw.V2d(pywinter_var_name.upper(), var_data))

            for var_name, var_data_list in variables_3d.items():
                var_data = var_data_list[0][n, :, :, :]
                plevs = plevs_3d[var_name]
                pywinter_var_name = map_to_pywinter_name(var_name)

                try:
                    total_fields.append(pyw.V3dp(pywinter_var_name.upper(), var_data, plevs))
                except Exception as e:
                    print(f"Error adding {var_name} to total_fields: {e}")

            if 'tsl' in variables_sl and 'mrsos' in variables_sl:
                sl_layer = ['000010', '010200']
                tsl_out = np.empty((2, variables_sl['tsl'][0].shape[1], variables_sl['tsl'][0].shape[2]))
                tsl_out[0, :, :] = variables_sl['tsl'][0][n, :, :]
                tsl_out[1, :, :] = variables_sl['tsl'][0][n, :, :]
                winter_soilt_layer = pyw.Vsl('ST', tsl_out, sl_layer)
                total_fields.append(winter_soilt_layer)

                mrsos_out = np.empty((2, variables_sl['mrsos'][0].shape[1], variables_sl['mrsos'][0].shape[2]))
                mrsos_out[0, :, :] = variables_sl['mrsos'][0][n, :, :]
                mrsos_out[1, :, :] = variables_sl['mrsos'][0][n, :, :]
                winter_soilm_layer = pyw.Vsl('SM', mrsos_out, sl_layer)
                total_fields.append(winter_soilm_layer)

            date_value = ncf.num2date(date_variable[n], date_variable.units)
            formatted_date = date_value.strftime('%Y-%m-%d_%H')

            output_file = pyw.cinter('MPI', formatted_date, winter_geo, total_fields, OUTPUT_DIR)
            # print(f"Written file: {output_file}")

# Main Execution Flow
def main():
    total_steps = 7308  # total time steps in the data
    chunk_size = 200
    process_chunks(total_steps, chunk_size, process_variables)

if __name__ == "__main__":
    main()
