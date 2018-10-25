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


# Opening files passed by command line args
nc4_data = netCDF4.Dataset(netcdf_file, 'r', format = "NETCDF4")
polygons = fiona.open(shape_file, 'r')


# Establishing required format definers to write a shape file.
# The driver and crs are borrowed from the one provided for consistency.
# The schema is kept simple to accept the co-ordinates of cells of interest
point_file_driver = polygons.driver
point_file_coordrefsys = polygons.crs
point_file_schema = {'geometry': 'Point',
					'properties': {'lon': 'int:4',
									'lat': 'int:4'}}

 
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
			point_lon = point["properties"]['lon']
			point_lat = point["properties"]["lat"]
			interest_longitudes.append(point_lon)
			interest_latitudes.append(point_lat)
			total_interest_cells += 1

print("Done Solving, found the following:")
print("Number of grid cells of interest: " + str(total_interest_cells))
print("Latitudes of interest:", interest_latitudes)
print("Longitude of interest:", interest_longitudes)



