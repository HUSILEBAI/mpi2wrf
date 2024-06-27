#!/bin/bash

# List of input files to process
input_files=(
    "tos_3hr_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_198001010300-198501010000.nc"
    "tos_3hr_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_198501010300-199001010000.nc"
    "tos_3hr_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_199001010300-199501010000.nc"
    "tos_3hr_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_199501010300-200001010000.nc"
    "tos_3hr_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_200001010300-200501010000.nc"
    "tos_3hr_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_200501010300-201001010000.nc"
    "tos_3hr_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_201001010300-201501010000.nc"
)

# Template for output file names
output_template="tos_6hrPlevPt_MPI-ESM1-2-HR_historical_r1i1p1f1_gn"

# Define the target grid description based on the tas grid
target_grid="target_grid.txt"

# Create the target grid description file
cat > $target_grid << EOL
gridtype = lonlat
gridsize = 73728
xsize = 384
ysize = 192
xfirst = -179.53125
xinc = 0.9375
yfirst = -89.4625
yinc = 0.9375
EOL

# Loop through each input file
for input_file in "${input_files[@]}"; do
    # Extract the time range part of the input file name
    time_range=$(echo $input_file | grep -oP '\d{12}-\d{12}')

    # Split the time range into start and end times
    start_time=${time_range:0:12}
    end_time=${time_range:13:12}

    # Adjust start time to 0600 hours
    start_date="${start_time:0:8}"
    adjusted_start_time="${start_date}0600"

    # Define intermediate and output file names
    intermediate_file="${input_file%.nc}_s.nc"
    output_file="${output_template}_${adjusted_start_time}-${end_time}.nc"

    # Select timesteps from the input file
    # cdo seltimestep,2/14608/2 selects timesteps at intervals from 2 to 14608 with a step of 2
    cdo seltimestep,2/14608/2 $input_file $intermediate_file

    # Perform interpolation of tos data to match the target grid (tas grid)
    cdo remapbil,$target_grid -selname,tos $intermediate_file $output_file

    # Remove the intermediate file to clean up
    rm $intermediate_file

    echo "Processed $input_file to $output_file"
done

# Remove the target grid file as it is no longer needed
rm $target_grid

