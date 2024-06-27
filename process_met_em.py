from netCDF4 import Dataset
import glob
import numpy as np

# Get a list of all "met_em.*.nc" files
file_list = glob.glob("met_em.*.nc")

for file_name in file_list:
    # Open the NetCDF file in read-write mode
    with Dataset(file_name, 'r+') as nc_file:
        # Set units for variables
        nc_file.variables['PRES'].units           = ''
        nc_file.variables['ST010200'].units       = 'K'
        nc_file.variables['SM010200'].units       = ' m3 m-3'
        nc_file.variables['ST000010'].units       = 'K'
        nc_file.variables['SM000010'].units       = ' m3 m-3'
        nc_file.variables['PMSL'].units           = 'Pa'
        nc_file.variables['SKINTEMP'].units       = 'K'
        nc_file.variables['PSFC'].units           = 'Pa'
        nc_file.variables['GHT'].units            = 'm'
        nc_file.variables['VV'].units             = 'm s-1'
        nc_file.variables['UU'].units             = 'm s-1'
        nc_file.variables['SST'].units            = 'K'
        nc_file.variables['SPECHUMD'].units       = 'kg kg-1'
        nc_file.variables['TT'].units             = 'K'

        nc_file.variables['PRES'].description     = ''
        nc_file.variables['ST010200'].description = ' 10-200 cm soil temp'
        nc_file.variables['SM010200'].description = ' 10-200 cm soil moisture'
        nc_file.variables['ST000010'].description = ' 0-10 cm soil temp'
        nc_file.variables['SM000010'].description = ' 0-10 cm soil moisture'
        nc_file.variables['PMSL'].description     = ' Mean sea-level pressure'
        nc_file.variables['SKINTEMP'].description = ' Skin temperature'
        nc_file.variables['PSFC'].description     = ' Surface pressure'
        nc_file.variables['GHT'].description      = ' 3-d geopotential height'
        nc_file.variables['VV'].description       = ' 3-d wind v-component'
        nc_file.variables['UU'].description       = ' 3-d wind v-component'
        nc_file.variables['SST'].description      = 'Sea surface temperature'
        nc_file.variables['SPECHUMD'].description = '3-d specific humidity'
        nc_file.variables['TT'].description       = '3-d air temperature'

        #Correct met_em SST according to land cover
        # Read the variables from the NetCDF file
        SST = nc_file.variables['SST'][:]  # Assuming 'SS' is the variable name for your field
        LM  = nc_file.variables['LANDMASK'][:]
        SKT = nc_file.variables['SKINTEMP'][:]
        # Convert to Kelvin as tos comes in Celsius
        celsiusE = np.where(SST != 0.0)
        SST[celsiusE] = SST[celsiusE]+273.15
        # Find indices where both SST and LM are zero
        noData = np.where((SST == 0.0) & (LM == 0.0))
        #Fill with values from SKIN Temperature (Ta from original variables)
        SST[noData] = SKT[noData]


        # Write the updated SST array back to the NetCDF file
        nc_file.variables['SST'][:] = SST

