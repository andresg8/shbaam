#!/usr/bin/env python
import sys
import os.path
import subprocess
import netCDF4
import numpy
import datetime
import fiona
import shapely.geometry
import shapely.prepared
import rtree
import math
import csv


# Parsing command line arguments, point_file can be made a command line arg too if needed
netcdf_file = sys.argv[1]
shape_file = sys.argv[2]
point_file = "../output/SERVIR_STK/GLDAS_VIC.pnt_tst.shp"

output_csv = "../output/SERVIR_STK/timeseries_NorthWestBD_test.csv"
output_nc4 = "../output/SERVIR_STK/map_NorthWestBD_test.nc4"


# Opening files passed by command line args
nc4_data = netCDF4.Dataset(netcdf_file, 'r', format = "NETCDF4")
polygons = fiona.open(shape_file, 'r')


# Establishing required format definers to write a shape file.
# The driver and crs are borrowed from the one provided for consistency.
# The schema is kept simple to accept the co-ordinates of cells of interest
point_file_driver = polygons.driver
point_file_coordrefsys = polygons.crs
point_file_schema = {'geometry': 'Point',
					'properties': {'lon': 'float:4.4',
									'lat': 'float:4.4'}}

 
print("Creating points shape file")

# Writing points shape file which holds a point for every combination 
# (thus 2 for loops) of latitude and longitude that exists in the nc4 file
# provided. These points can later be used to retrieve further data from
# the .nc4 file once the relevant points are filtered.
with fiona.open(point_file, 'w', driver = point_file_driver,
						crs = point_file_coordrefsys,
						schema = point_file_schema) as points:
	for lon_index in range(len(nc4_data.dimensions["lon"])):
		longitude = nc4_data.variables["lon"][lon_index]
		for lat_index in range(len(nc4_data.dimensions["lat"])):
			latitude = nc4_data.variables["lat"][lat_index]
			schema_properties = {"lon": longitude,
									"lat": latitude}
			schema_geometry = shapely.geometry.mapping(
								shapely.geometry.Point((longitude, latitude)))
			points.write({"properties": schema_properties,
									"geometry": schema_geometry,})

print("Success - points created")

print("Creating rtree for points")

# Opening newly created point_file to produce an rtree index
points = fiona.open(point_file, 'r')

# Loading points into an rtree index with the shapely default value 
# for bounds. The index assists in solving for intersections.
index = rtree.index.Index()
for point in points:
	point_id = int(point['id'])
	shape = shapely.geometry.shape(point['geometry'])
	index.insert(point_id, shape.bounds)

print("Tree has been created")

print("Solving for grid cells of interest")

total_interest_cells = 0
interest_longitudes = []
interest_latitudes = []


# Checking every polygon in the shapefile provided for any intersections 
# with the points that were loaded into the rtree index above. If there 
# are intersections we check that the point actually lies within the
# polygon, because the bounds are a default value and thus it is possible
# a point outside of our area of interest intersects with the bounds of
# the polygon. If a point's default bounds intersect with the polygon's and
# the point lies within the polygon, it is added as a locator for the cell 
# of interest
for polygon in polygons:
	poly_geom = shapely.geometry.shape(polygon['geometry'])
	prep_geom = shapely.prepared.prep(poly_geom)
	for point_id in [int(x) for x in list(index.intersection(poly_geom.bounds))]:
		point = points[point_id]
		point_geom = shapely.geometry.shape(point['geometry'])
		if prep_geom.contains(point_geom):
			print(point["geometry"])
			point_lon = point["properties"]['lon']
			point_lat = point["properties"]["lat"]
			interest_longitudes.append(point_lon)
			interest_latitudes.append(point_lat)
			total_interest_cells += 1

print("Done Solving, found the following:")
print("Number of grid cells of interest: " + str(total_interest_cells))
print("Latitudes of interest:", interest_latitudes)
print("Longitude of interest:", interest_longitudes)

print("Solving for long term means for values of SWE across the cells of interest")


# Previously, the actual values for interest latitude and longitude were solved for.
# However, the indices of the values relative to the nc4 data file are required to 
# extract values for SWE. Eg: extracting the SWE value at 179.5 W, 59.5 S is accessed
# by the index [0, 0] because the indices for the example lon, lat are 0. 
# Here, the indices for the values of interest are extracted from the nc4 data file.
int_lat_indices = []
int_lon_indices = []
for index in range(total_interest_cells):
	int_lon_indices.append((nc4_data["lon"][:]).tolist().index(interest_longitudes[index]))
	int_lat_indices.append((nc4_data["lat"][:]).tolist().index(interest_latitudes[index]))


# Using the indices solved for above, the relevant SWE values at every available time
# step are averaged and compiled into a list. 
time_steps = len(nc4_data.dimensions["time"])
long_term_means = [0] * total_interest_cells
for index in range(total_interest_cells):
     lon_index = int_lon_indices[index]
     lat_index = int_lat_indices[index]
     for time_step in range(time_steps):
          long_term_means[index]=long_term_means[index]                        \
                                +nc4_data.variables['SWE']                  \
                                            [time_step,lat_index,lon_index]
long_term_means=[x/time_steps for x in long_term_means]


print(long_term_means)


longitude_step_size = abs(nc4_data["lon"][1]-nc4_data["lon"][0])
latitude_step_size = abs(nc4_data["lat"][1]-nc4_data["lat"][0])


areas= []
for lat_index in int_lat_indices:
     latitude = nc4_data["lat"][lat_index]
     areas.append(6371000*math.radians(latitude_step_size)               \
                           *6371000*math.radians(longitude_step_size)               \
                           *math.cos(math.radians(latitude)))
total_area = sum(areas)


anomalies=[]
for time_step in range(time_steps):
     anomaly_in_time=0
     for index in range(total_interest_cells):
          lon_index=int_lon_indices[index]
          lat_index=int_lat_indices[index]
          area=areas[index]
          long_term_mean=long_term_means[index]
          anomaly_in_area=(nc4_data.variables['SWE']                            \
                                  [time_step,lat_index,lon_index]          \
                      -long_term_mean                                    )/100     \
                    *  area
          anomaly_in_time+=anomaly_in_area
     anomalies.append(100*anomaly_in_time / total_area)

print('Check some computations')

print('- Average of time series: '+str(numpy.average(anomalies)))
print('- Maximum of time series: '+str(numpy.max(anomalies)))
print('- Minimum of time series: '+str(numpy.min(anomalies)))

