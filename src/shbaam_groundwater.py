#!/usr/bin/env python
import sys
from netCDF4 import Dataset, num2date, date2num
import fiona
import shapely.geometry
import shapely.prepared
import rtree
import math
import numpy
import csv
import datetime

def createPointShp(input_shp, output_shp, num_lon, num_lat, lons, lats):
    '''
    Create a point with all GLDAS grid cells
    '''
    point_driver = input_shp.driver
    point_crs = input_shp.crs

    point_schema = {'geometry': 'Point',
                 'properties': {'lon' : 'float',
                                'lat' : 'float',
                                'lon_index': 'int:4',
                                'lat_index': 'int:4',
                                }}

    with fiona.open( output_shp, 'w', driver=point_driver, crs=point_crs, schema=point_schema) as pf:
       for lon_index in range(num_lon):
            longitude = lons[lon_index]

            for lat_index in range(num_lat):
                latitude = lats[lat_index]

                pf_property = { 'lon': longitude, 'lat': latitude,
                                'lon_index': lon_index, 'lat_index': lat_index,}
                pf_geometry = shapely.geometry.mapping(
                            shapely.geometry.Point((longitude, latitude)))
                pf.write({
                           'properties': pf_property,
                           'geometry': pf_geometry,
                         })

    print('Success -- created a new shapefile')

def createRtreeIndex(pf):
    idx = rtree.index.Index()
    for point in pf:
        point_id = int(point['id'])
        point_bounds = shapely.geometry.shape(point['geometry']).bounds
        idx.insert(point_id, point_bounds)
    print('Success -- created rtree')
    return idx

def findInterest(pol, pf, idx):
    total = 0
    interest_lon = []
    interest_lat = []

    for polygon in pol:
        polygon_geo = shapely.geometry.shape(polygon['geometry'])
        prepared_polygon = shapely.prepared.prep( polygon_geo )
        for pf_id in [int(i) for i in list(idx.intersection(polygon_geo.bounds))]:
            point = pf[pf_id]
            point_bounds = shapely.geometry.shape(point['geometry'])
            if prepared_polygon.contains(point_bounds):
                #print('cell:', pf_id)
                interest_lon.append( (point['properties']['lon_index'], point['properties']['lon']) )
                interest_lat.append( (point['properties']['lat_index'], point['properties']['lat']) )
                total += 1

    print('The total cells found: ' + str(total) )
    print('The longitude of the cells: ' + str([i[1] for i in interest_lon]) )
    print('The lattitude of the cells: '+ str([i[1] for i in interest_lat]) )
    return total, interest_lon, interest_lat

def findAvg(var, total, interest_lon, interest_lat, time_steps):
    #var: the nc4.variables that you want to get for average such as swe..
    print('--Find long-term mean for each intersecting GRACE grid cell--')
    avg = [0.0] * total
    for cell_index in range(total):
        lon_index = interest_lon[cell_index][0]
        lat_index = interest_lat[cell_index][0]
        for time in range(time_steps):
            avg[cell_index] = avg[cell_index] + var[time, lat_index, lon_index]

    avg = [i / time_steps for i in avg]
    print('--long-term mean:' + str(avg)+'--')
    return avg

def calculateSurfaceArea(total, num_lat, interest_lat, lon_step, lat_step):
    areas = [0] * total
    for cell_index in range(total):
        areas[ cell_index ] = 6371000*math.radians(lat_step) \
                                    *6371000*math.radians(lon_step)\
                                    *math.cos(math.radians(interest_lat[cell_index][1]))
    print('area for each interest cell: ' + str(areas))
    return areas

def anomalyTimeseries(vals, var_factor, avg, time_steps, total, interest_lon, interest_lat, areas):
    print('Compute storage anomaly timeseries')
    anomalies = []
    total_area = sum(areas)
    for time in range(time_steps):
        anomaly_in_time = 0
        for cell_index in range(total):
            lon_index = interest_lon[cell_index][0]
            lat_index = interest_lat[cell_index][0]
            area = areas[cell_index]
            long_term_mean = avg[cell_index]

            anomaly_in_area = ( vals[time, lat_index, lon_index] \
                                - long_term_mean) / var_factor * area
            anomaly_in_time += anomaly_in_area
        anomalies.append(var_factor * anomaly_in_time / total_area)

    avgT = numpy.average(anomalies)
    maxT = numpy.max(anomalies)
    minT = numpy.min(anomalies)
    print('- Average of time series: '+str(avgT))
    print('- Maximum of time series: '+str(maxT))
    print('- Minimum of time series: '+str(minT))

    return anomalies

def createTimes(times):
    dates = num2date(times[:], units=times.units, calendar=times.calendar)
    time_stamps = []
    for i in dates:
        d = datetime.datetime.strptime(str(i), '%Y-%m-%d %H:%M:%S')
        time_stamps.append(d.strftime('%m/%d/%Y'))
    return time_stamps

def outputCSV(output_file, fieldnames, times, anomalies_dict):
    with open(output_file, mode='w') as csvfile:
        writer = csv.DictWriter(csvfile, dialect='excel', fieldnames=fieldnames)
        dates = createTimes(times)
        writer.writeheader()
        for date_idx in range(len(dates)):
            d = dict()
            d[fieldnames[0]] = dates[date_idx]
            for field_idx in range(1, len(fieldnames)):
                d[fieldnames[field_idx]] = anomalies_dict[fieldnames[field_idx]][date_idx]
            writer.writerow(d)
    csvfile.close()
    print('Success -- creating csv file')

def outputNC(output_file, nc4_file, total, interest_lon, interest_lat, time_steps, var_list, avg_dict):
    print('Writing new nc4 file')
    nc4_out = Dataset(output_file, 'w', format='NETCDF4')

    #copying dimensions
    for (name, dim) in nc4_file.dimensions.items():
        nc4_out.createDimension(name, len(dim) if not dim.isunlimited() else None)

    #copying variables with their attributes
    for (name, value) in nc4_file.variables.items():
        nc4_vars = nc4_out.createVariable(name, value.datatype, value.dimensions)
        nc4_vars.setncatts(nc4_file[name].__dict__)
        if name in {"lat", "lon", "time"}:
            nc4_out[name][:] = nc4_file[name][:]

    dt = datetime.datetime.utcnow()
    dt = dt.replace(microsecond=0)

    nc4_out.Conventions='CF-1.6'
    nc4_out.title=''
    nc4_out.institution=''
    nc4_out.history='date created: '+ dt.isoformat() +'+00:00'
    nc4_out.references='https://github.com/c-h-david/shbaam/'
    nc4_out.comment=''
    nc4_out.featureType='timeSeries'

    for var in var_list:
        for i in range(total):
            lon_index, lat_index = interest_lon[i][0],interest_lat[i][0]
            long_term_mean = avg_dict[var][i]
            for time in range(time_steps):
                nc4_out[var][time,lat_index,lon_index] = nc4_file.variables[var][time,lat_index,lon_index] - long_term_mean
    nc4_out.close()
    print("Success -- creating nc4 file")

if __name__ == "__main__":
    files = sys.argv[1:]
    """
    files[0]: concatenating netCDF4 files
    files[1]: given shapefile ('../input/SERVIR_STK/Nepal.shp')
    files[2]: output shapefile with all grid cells from nc4
    files[3]: output csv file
    files[4]: output nc file
    files[5]: specify the variables
    """
    #open netCDF4 file
    nc4_file = Dataset(files[0], 'r', format="NETCDF4") #concatenating files

    num_lon = len(nc4_file.dimensions['lon'])
    num_lat = len(nc4_file.dimensions['lat'])
    time_steps = len(nc4_file.dimensions['time'])

    lons = nc4_file.variables['lon']
    lats = nc4_file.variables['lat']
    times = nc4_file.variables['time']

    lon_step = abs(lons[1] - lons[0])
    lat_step = abs(lats[1] - lats[0])

    all_vars = []
    dimensional = ["lon", "lat", "time"]
    for var in nc4_file.variables.keys():
        if var.lower() not in dimensional:
            all_vars.append(var)

    done = False
    var_list = []
    while not done:
        for i in range(len(all_vars)):
            print(i, all_vars[i])
        print(len(all_vars), "All")
        print(len(all_vars)+1, "Done")
        var = raw_input("What variables should be included: ")
        if var.lower() == "done":
            done = True 
        elif var.lower() == "all":
            done = True
            var_list = all_vars
        elif var in all_vars:
            var_list.append(var)
            print(var + " has been added, current list =", var_list)
        else:
            print("Invalid input, please try again")

    unit_factors = {"mm": 1000, "kg/m2": 1000, "cm": 100, "dm": 10, "m": 1, "km": .001}


    #read polygon shapefile
    polygon = fiona.open(files[1],'r')
    point_file = files[2]

    createPointShp(polygon, point_file, num_lon, num_lat, lons, lats)

    #open newly created point_file
    pf = fiona.open(point_file, 'r')

    #create spatial index for the bounds of each point
    idx = createRtreeIndex(pf)

    total_interest, interest_lon, interest_lat = findInterest(polygon, pf, idx)

    #calculate surface area for each interest cell
    areas = calculateSurfaceArea(total_interest, num_lat, interest_lat, lon_step, lat_step)

    #find long-term mean for swe
    #var_list= files[5:]
    avg_dict = dict()
    anomalies_dict = dict()
    for var in var_list:
        vals = nc4_file.variables[var] #all values in nc4_file
        var_factor = unit_factors[nc4_file[var].__dict__["units"]]
        avg = findAvg(vals, total_interest, interest_lon, interest_lat, time_steps)
        anomalies = anomalyTimeseries(vals, var_factor, avg, time_steps, total_interest, interest_lon, interest_lat, areas)
        avg_dict[var] = avg
        anomalies_dict[var] = anomalies

    output_csv = files[3]
    fieldname = ['date'] + var_list
    #print(fieldname)
    outputCSV(output_csv, fieldname, times, anomalies_dict)

    output_nc = files[4]
    outputNC(output_nc, nc4_file, total_interest, interest_lon, interest_lat, time_steps, var_list, avg_dict)

    nc4_file.close()
    polygon.close()
pf.close()