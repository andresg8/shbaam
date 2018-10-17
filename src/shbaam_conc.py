#!/usr/bin/env python
from netCDF4 import Dataset, date2num, num2date
import os
import sys

path = sys.argv[1]
#targets = [(os.path.dirname(path) + f) for f in os.listdir(os.path.dirname(path)) if f.endswith(os.path.basename(path)[1:])]
targets = sys.argv[1:-1]
destination = sys.argv[-1]

# This script assumes all files to be concatenated contain an identical layout of data
# thus any random (in this case the first) file can be used as a template to set the 
# general attributes as a foundation for the data to be hosted. This assumption is made 
# because I can't see the usefuleness of concatenating files hosting varying unrelated 
# data and because, for the specificity of this script's purpose, the assumption should 
# have no negative effect
template = Dataset(targets[0], 'r', format="NETCDF4")

# Pacman represents the aggregate file while the files being concatenated are the
# pellets he's "eating"
pacman = Dataset(destination, 'w', format="NETCDF4")

# Setting Global Attributes based on template file
pacman.setncatts(template.__dict__)

# Setting Dimensions based on template file
for dName, dSize in template.dimensions.items():
	pacman.createDimension(dName, (None if dSize.isunlimited() else len(dSize)))

# Setting Variable Attributes and Setting values for Latitude and Longitude
for vName, vValue in template.variables.items():
	ref = pacman.createVariable(vName, vValue.datatype, vValue.dimensions)
	pacman[vName].setncatts(template[vName].__dict__)
	if vName in ["lat", "lon"]:
		pacman[vName][:] = template[vName][:]

# Concatenating values for Canint, SWE, and time
count = 0
for target in targets:
	pellet = Dataset(target, 'r', format="netCDF4")
	for vName in pellet.variables.keys():
		if vName in ["Canint", "SWE"]:
			pacman[vName][count] = pellet[vName][:]
		if vName == "time":
			pelletTime = num2date(pellet[vName][:], pellet[vName].units, pellet[vName].calendar)[0]
			pacman[vName][count] = date2num(pelletTime, pacman[vName].units, pacman[vName].calendar)
	count += 1

	pellet.close()

# Done
pacman.close()