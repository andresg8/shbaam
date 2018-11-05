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
point_file = sys.argv[3]
output_csv = sys.argv[4]
output_nc4 = sys.argv[5]


# Opening files passed by command line args
nc4_data = netCDF4.Dataset(netcdf_file, 'r', format = "NETCDF4")
polygons = fiona.open(shape_file, 'r')

# Assigning easier variable names to commonly accessed nc4_data
lon_values = nc4_data.variables["lon"]
lat_values = nc4_data.variables["lat"]
time_values = nc4_data.variables["time"]
swe_values = nc4_data.variables['SWE']
lon_dimension = len(nc4_data.dimensions["lon"])
lat_dimension = len(nc4_data.dimensions["lat"])
time_dimension = len(nc4_data.dimensions["time"])


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
  for lon_index in range(lon_dimension):
    longitude = lon_values[lon_index]
    for lat_index in range(lat_dimension):
      latitude = lat_values[lat_index]
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



print("Solving for long term means for values of SWE across the cells of interest")
# Previously, the actual values for interest latitude and longitude were solved for.
# However, the indices of the values relative to the nc4 data file are required to 
# extract values for SWE. Eg: extracting the SWE value at 179.5 W, 59.5 S is accessed
# by the index [0, 0] because the indices for the example lon, lat are 0. 
# Here, the indices for the values of interest are extracted from the nc4 data file.
int_lat_indices = []
int_lon_indices = []
for index in range(total_interest_cells):
  int_lon_indices.append((lon_values[:]).tolist().index(interest_longitudes[index]))
  int_lat_indices.append((lat_values[:]).tolist().index(interest_latitudes[index]))


# Using the indices solved for above, the relevant SWE values at every available time
# step are averaged and compiled into a list. 
long_term_means = [0] * total_interest_cells
for index in range(total_interest_cells):
     lon_index = int_lon_indices[index]
     lat_index = int_lat_indices[index]
     for time_step in range(time_dimension):
          long_term_means[index] = (long_term_means[index] +
                                  swe_values[time_step,lat_index,lon_index])
long_term_means = [x / time_dimension for x in long_term_means]


# The areas of the cells of interest can be evaluated independently based on the
# span of every longitudinal and latitudinal value. A smaller step size implies
# smaller available cells. The total_area variable represents the entire area across
# all cells and is used to determine the numeric impact of a change over the whole 
# land mass of interest. 
longitude_step_size = abs(lon_values[1] - lon_values[0])
latitude_step_size = abs(lat_values[1] - lat_values[0])
areas= []
for lat_index in int_lat_indices:
     latitude = lat_values[lat_index]
     areas.append((6371000*math.radians(latitude_step_size) *             
                           6371000*math.radians(longitude_step_size) *
                           math.cos(math.radians(latitude))))
total_area = sum(areas)


# Now that we have the area of each grid cell and the average value for the 
# snow water equivalent present, we can calculate the change of swe in each 
# cell over the time span available and relative to the total area we're 
# interested in. 
anomalies=[]
for time_step in range(time_dimension):
     anomaly_in_time = 0
     for index in range(total_interest_cells):
          lon_index = int_lon_indices[index]
          lat_index = int_lat_indices[index]
          area = areas[index]
          long_term_mean = long_term_means[index]
          anomaly_in_area = ((swe_values[time_step,lat_index,lon_index] - 
                            long_term_mean)/100 * area)
          anomaly_in_time += anomaly_in_area
     anomalies.append(100 * anomaly_in_time / total_area)



print('Writing CSV file')
# Using a predetermined datetime as a zero, time stamps are created by 
# adding an existing time to the zero. The time stamps are used to write
# their corresponding data into one row on the excel sheet. 
time_stamp_zero = datetime.datetime.strptime('2002-04-01T00:00:00',                
                                         '%Y-%m-%dT%H:%M:%S')
time_stamps = []
for time_step in range(time_dimension):
     delta_t = datetime.timedelta(hours = time_values[time_step])
     time_stamp = (time_stamp_zero + delta_t).strftime('%m/%d/%Y')
     time_stamps.append(time_stamp)

# Writing excel sheet with timestamped data
with open(output_csv, 'wb') as csvfile:
     csvwriter = csv.writer(csvfile, dialect = 'excel')
     for time_step in range(time_dimension):
          row = [time_stamps[time_step],anomalies[time_step]] 
          csvwriter.writerow(row) 



print("Writing new nc4 file")
SWE_nc4 = netCDF4.Dataset(output_nc4, 'w', format="NETCDF4")

# Setting the global attributes for the new nc4 file. 
# Metadata format and values kept consistent with available example. 
dt=datetime.datetime.utcnow()
dt=dt.replace(microsecond=0)
vsn=subprocess.Popen('bash ../version.sh',                                     
                     stdout=subprocess.PIPE,shell=True).communicate()
vsn=vsn[0]
vsn=vsn.rstrip()
SWE_nc4.Conventions='CF-1.6'
SWE_nc4.title=''
SWE_nc4.institution=''
SWE_nc4.source='SHBAAM: '+vsn+', GLDAS_VIC: '+os.path.basename(netcdf_file)   
SWE_nc4.history='date created: '+dt.isoformat()+'+00:00'
SWE_nc4.references='https://github.com/c-h-david/shbaam/'
SWE_nc4.comment=''
SWE_nc4.featureType='timeSeries'


# Setting Dimensions based on nc4 file 
for dName, dSize in nc4_data.dimensions.items():
  SWE_nc4.createDimension(dName, (None if dSize.isunlimited() else len(dSize)))

# Setting Variable Attributes and Setting values for Latitude and Longitude
for vName, vValue in nc4_data.variables.items():
  ref = SWE_nc4.createVariable(vName, vValue.datatype, vValue.dimensions)
  SWE_nc4[vName].setncatts(nc4_data[vName].__dict__)
  if vName in ["lat", "lon", "time"]:
    SWE_nc4[vName][:] = nc4_data[vName][:]

''' COMMENTED OUT TEMPORARILY
print('- Modify CRS variable attributes')
swe.grid_mapping='crs'
crs.grid_mapping_name='latitude_longitude'
crs.semi_major_axis='6378137'
crs.inverse_flattening='298.257223563' 
'''
# Writting values into new nc4 file using only the values extracted from
# our grid cells of interest.
for index in range(total_interest_cells):
     lon_index=int_lon_indices[index]
     lat_index=int_lat_indices[index]
     long_term_mean=long_term_means[index]
     for time_step in range(time_dimension):
          SWE_nc4["SWE"][time_step,lat_index,lon_index] = (
            swe_values[time_step,lat_index,lon_index] - long_term_mean)



# Closing nc4 files
nc4_data.close()
SWE_nc4.close()


print('Check some computations')

print('- Average of time series: '+str(numpy.average(anomalies)))
print('- Maximum of time series: '+str(numpy.max(anomalies)))
print('- Minimum of time series: '+str(numpy.min(anomalies)))

