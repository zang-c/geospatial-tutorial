# -*- coding: utf-8 -*-
"""
Created on 2025.07.07

Purpose: 
    This code is for a Geospatial Analysis tutorial for group meeting.
    It goes throught a basic plotting example that shows some useful packages and how to use them.

    Then we take some of these skills to work with larger xarray.Dataset objects containing sea ice concentration data.
        This will demonstate the building and applicaiton of shapefiles along with some bulk analysis on Dataset objects

@author: Cort L. Zang
"""
#%%
import numpy as np
import os
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import matplotlib.pyplot as plt
import cartopy as cart
from cartopy.geodesic import Geodesic
import cartopy.io.img_tiles as cimgt
import geopandas
from shapely.geometry import Polygon, mapping, Point
from datetime import datetime, timedelta

#%% Pull data from https://www.naturalearthdata.com/ 
'''
Use cartopy.feature.NaturalEarthFeature function to pull some data to add to our map that we'll create later.

define category as 'cultural' or 'physical' (that would probably be better named 'human-related' and 'natural' respectively)
find the file names on https://github.com/nvkelso/natural-earth-vector and omit preamble string 'ne_**m_'

set the scale (resolution) of the data you want to pull. 
    These are generally smaller files so I usually default to the highest (10m) but some have 50m and 110m available

It can be nice to define plotting characterisics when you pull the data as well but this can also be done during plotting.
'''
resol = '10m'
provinc_bodr = cart.feature.NaturalEarthFeature(category='cultural', 
                                                name='admin_1_states_provinces_lines', 
                                                scale=resol, 
                                                facecolor='none', 
                                                edgecolor='k')
urban_area = cart.feature.NaturalEarthFeature(category='cultural', 
                                         name = 'urban_areas_landscan', 
                                         scale = resol, 
                                         facecolor = "#FCC87B",
                                         edgecolor = 'none')
rivers = cart.feature.NaturalEarthFeature(category='physical',
                                          name = 'rivers_north_america',
                                          scale = resol,
                                          facecolor = "none",
                                          edgecolor = "#93b7d4",
                                          ls = ':')
roads = cart.feature.NaturalEarthFeature(category='cultural',
                                          name = 'roads',
                                          scale = '10m',
                                          facecolor = 'none',
                                          edgecolor = '#000000')

imagery = cimgt.OSM() # This defines a layered set of tile images using OpenStreetMap data (https://www.openstreetmap.org/about)

'''
Create a function that plots the different data you selected onto a given fig.axes
'''
def add_map_features(ax):
    ax.coastlines(lw = 0.5)
    ax.add_feature(cfeature.BORDERS, lw = 0.5)
    ax.add_feature(provinc_bodr, lw = 0.5)
    ax.add_feature(rivers, lw = 0.25)
    ax.add_feature(cfeature.RIVERS, lw = 0.25)
    ax.add_feature(cfeature.LAKES, alpha = 0.95)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.LAND, color = "#fff9da")
    ax.add_feature(urban_area)
    ax.add_feature(roads, lw = 0.05)


#%% Define significant locations in your life!
'''
Next we will define some coordinates that are important to us. You can add as many points as you want; I've added three to include
1. Where I grew up, 2. Where I went to college, and 3. Where I'm going to graduate school. For each of these coordinate sets define a label.
'''
lon_1 = -77.6997
lat_1 = 39.4862
label_1 = 'Keedysville, MD'

lon_2 = -75.5994
lat_2 = 38.3607
label_2 = 'Salisbury, MD'

lon_3 = -105.0844
lat_3 = 40.5853
label_3 = 'Fort Collins, CO'

'''
Now, we make lists of each set of Longitudes, Latitudes and labels, as well as define a list of colors to plot with.
    Each of these lists should be the same length.
'''
lons = [lon_1, lon_2, lon_3]
lats = [lat_1, lat_2, lat_3]
labs = [label_1, label_2, label_3]
colors = ["#83BE77", "#BB87D3", "#607BF1"]

# Here I am calculating the longitude (and latitude if needed) that will be the center of our map.
lon_min, lon_max = np.min(lons), np.max(lons)
lat_min, lat_max = np.min(lats), np.max(lats)
lon_middle = (((lon_max + 360) + (lon_min + 360)) / 2) - 360
lat_middle = (((lat_max + 180) + (lat_min + 180)) / 2) - 180
#%% Plot your data
'''
Let's create a figure and then make a subplot within that figure.

Since we are plotting with geospatial data we want to define a projection for the maps.
    Different projections can be found at https://scitools.org.uk/cartopy/docs/v0.15/crs/projections.html

I've started with a Mercator projection which is probably the most 'standard' mapping projection. We'll play around with some others later.
'''
fig = plt.figure(dpi = 300)
ax = fig.add_subplot(projection = ccrs.Mercator(central_longitude=lon_middle))

ax.plot(lons, lats,
            color='crimson', linewidth=1,
            transform=ccrs.PlateCarree())

for lon, lat, lab, color in zip(lons, lats, labs, colors):
    ax.scatter(lon, lat, marker = 'o', s = 10, 
           color = color, label = lab,
           transform = ccrs.PlateCarree(), 
           zorder = 10, edgecolor = 'white', lw = 0.5)

# Choose either our custom map features or OSM mapping! Play around with both, but for the OSM stay at a low zoom and work your way up.
map_features = 'custom'
zoom = 1 
if map_features == 'osm':
    ax.add_image(imagery, zoom)
elif map_features == 'custom':
    add_map_features(ax)
elif map_features == '':
    pass
else:
    print('No recognized mapping style given. Passing')
    pass

# Now we set the extents of our map and create a legend
ax.set_extent([lon_min - 5, lon_max + 5, lat_min - 5, lat_max + 5])
ax.legend(markerscale = 2)

'''
That should get you familar with some of the basics of using matplotlib to create maps.
Things to try:
    Play around with different projections and see how they distort data and such. 
    Try to pull different data from naturalearthdata.com and include it in your maps!
'''
#%% Create a timeseries xr.Dataset of seaice concentrations to save as a netcdf
def concat_seaice():
  '''
  This is the code that I used to make the seaice.nc file that is in our Github repository.
  I pulled the raw data from University of Bremen (https://data.seaice.uni-bremen.de/amsr2/asi_daygrid_swath/n3125/netcdf/2023/)
    Here I am just merging them into one xr.Dataset object and downsampling it to 25km resolution using rioxarray (.rio)
    I then create a timeseries of dates associated with each layer of the xarray using the filename and save the Dataset as a netcdf (.nc) file
  '''
  filepaths = []
  path = r"E:\artofmelt\20230500_ARTofMELT\Data\aom_backtraj_products\data\sea_ice_20250401\\" # Path to folder containing raw UBremen seaice data
  nc = r"E:\misc-projects\group-meetings\geospatial_tutorial\seaice.nc" # Path where we want to save our seaice timeseries
  for root, dirs, files in os.walk(os.path.abspath(path)):
      for file in files:
          filepath = os.path.join(root, file)
          filepaths.append(filepath)
      print(filepaths)
  seaice = xr.open_mfdataset(filepaths, concat_dim= 'time', combine = 'nested', parallel= True, decode_cf = False, engine = 'netcdf4')
  seaice = seaice.rio.write_crs("EPSG:3413")
  seaice = seaice.z.rio.reproject('EPSG:3413', resolution = 25000)
  datetimes = []
  strdatetimes = []
  for timestamp in filepaths:
      timestamp = timestamp[-16:-8]
      strdatetimes.append(timestamp)
  seaice['time'] = strdatetimes
  seaice.to_netcdf(nc, mode = 'w')
  return seaice
# %% Load our netcdf seaice concentration timeseries and make sure it is formatted correctly
'''
Let's open a timeseries seaice concentration file!
    We need to get our timeseries as a series of datetime objects
    then we need to inform the computer what the projection (CRS) of this data is (this sea ice data is in North Polar Stereographic a.k.a. 'EPSG:3413'
'''
seaice = xr.open_dataset(r"E:\misc-projects\group-meetings\geospatial_tutorial\seaice.nc")
seaice['time'] = pd.to_datetime(seaice.time, format = '%Y%m%d')
seaice = seaice.z.rio.write_crs("EPSG:3413")

#%% Define polygons to convert and save as a shapefiles for sea ice concentration analysis
'''
Here we will define a polygon to use as our North East Water (NEW) Polynya
'''
lon1, lat1 = -15, 81.5 # north west corner
lon2, lat2 = -10, 81.5 # north east corner
lon3, lat3 = -7.5, 79 # south east corner
lon4, lat4 = -10, 79 # south west corner

new_coords = [(lon1, lat1), 
              (lon2, lat2), 
              (lon3, lat3), 
              (lon4, lat4), 
              (lon1, lat1)] 
new = Polygon(new_coords)
gdf = geopandas.GeoDataFrame(geometry=[new], crs="EPSG:4326")
gdf = gdf.to_crs('EPSG:3413')
new_outfile_path = r'E:\misc-projects\group-meetings\geospatial_tutorial\new-polynya.shp'
gdf.to_file(new_outfile_path)

'''
Now we will define the Fram Stait region and then subtract the NEW Polynya polygon from the Fram Strait polygon
'''
lon1, lat1 = -22, 82 # north west corner
lon2, lat2 = -5, 82 # north east corner
lon3, lat3 = -5, 78 # south east corner
lon4, lat4 = -22, 78 # south west corner

fram_coords = [(lon1, lat1), 
               (lon2, lat2), 
               (lon3, lat3), 
               (lon4, lat4), 
               (lon1, lat1)] 
fram = Polygon(fram_coords)
fram = fram - new
gdf = geopandas.GeoDataFrame(geometry=[fram], crs="EPSG:4326")
gdf = gdf.to_crs('EPSG:3413')
fram_outfile_path = r'E:\misc-projects\group-meetings\geospatial_tutorial\fram-strait.shp'
gdf.to_file(fram_outfile_path)

'''
Now we will define the Arctic region, here I use a high Arctic threshold (1.5 million meters from the North Pole)
    This code was repurposed from a stackoverflow answer by users rwalsh3750 and MPA (https://stackoverflow.com/a/58735566)
'''
circle_points = Geodesic().circle(lon=0, 
                                  lat=90, 
                                  radius=1.5E6, 
                                  n_samples=64, 
                                  endpoint=False)
arctic = Polygon(circle_points)
gdf = geopandas.GeoDataFrame(geometry=[arctic], 
                             crs="EPSG:4326")
gdf = gdf.to_crs('EPSG:3413')
arctic_outfile_path = r'E:\misc-projects\group-meetings\geospatial_tutorial\arctic.shp'
gdf.to_file(arctic_outfile_path)
#%% Reload the the shapefiles and pull seaice concentration data associated with each shapes region
'''
Here we reload the shapes we defined above and 'clip' the seaice dataset into their regional definitions.
    This is essentially just making multiple datasets of the seaice concentration data focused on different regions
'''
new_shp = geopandas.read_file(new_outfile_path, crs = 'EPSG:3413')
fram_shp = geopandas.read_file(fram_outfile_path, crs = 'EPSG:3413')
arctic_shp = geopandas.read_file(arctic_outfile_path, crs = 'EPSG:3413')

polynya_ice = seaice.rio.clip(new_shp.geometry.apply(mapping), new_shp.crs, drop = False)
seaice_fram = seaice.rio.clip(fram_shp.geometry.apply(mapping), fram_shp.crs, drop = False)
arctic_ice = seaice.rio.clip(arctic_shp.geometry.apply(mapping), arctic_shp.crs, drop = False)
#%% Analyze the data!
'''
Here I want to understand how the concentration of seaice in the different regions changes over time.
'''
fig, ax = plt.subplots()
arctic_ice.mean(dim=['x','y']).plot(marker = '+', 
                                    color = '#000000', 
                                    label = 'Arctic Mean Ice Conc.', 
                                    ax = ax)
polynya_ice.mean(dim=['x', 'y']).plot(marker = 'o', 
                                      color = "#4FA3BD", 
                                      label = 'Polynya Mean Ice Conc.',
                                      ax = ax)
seaice_fram.mean(dim=['x','y']).plot(marker = 'o', 
                                     color = "#D4624E", 
                                     label = 'Fram Strait Mean Ice Conc.',
                                     ax = ax)
ax.set_ylabel('Ice Fraction')
ax.set_title('')
ax.xaxis.set_minor_locator(matplotlib.dates.DayLocator(interval = 1))
ax.xaxis.set_major_formatter(matplotlib.dates.ConciseDateFormatter(matplotlib.dates.AutoDateLocator()))
ax.legend()
plt.show()
#%%
'''
This plot should show that the NEW experiences a strong reduction in seaice concentration from May 15th, 2023 to June 1st, 2023.

Let's look at these dates with more detail...

First we define a new dataset of the seaice data (over the entire globe) and select the dates we are interested in from the 'time' dimension
    Then for plotting I want to reproject the data.

Next I want to plot the data under each day of interest and see how things change between days and over the whole time period
'''

new_polynya_delta = seaice[20:38, :,:] # This pulls seaice concentration data for the dates from May 15th through June 1st, 2023
new_polynya_delta = new_polynya_delta.rio.reproject('EPSG:4326')


fig = plt.figure(figsize = (30, 15))

for day in range(0, len(new_polynya_delta.time)):
    ax = fig.add_subplot(3, 
                         6, 
                         day + 1, 
                         projection = ccrs.Orthographic(central_longitude=0, 
                                                                       central_latitude=80))
    diff = new_polynya_delta[day, :,:] - new_polynya_delta[day - 1, :,:]
    diff.plot(ax = ax, 
              transform = ccrs.PlateCarree(), 
              cmap = 'coolwarm', 
              add_colorbar=True, 
              cbar_kwargs={"label": "sea ice conc. difference between days",
                           'shrink': 0.6})
    ax.set_extent([-20,20,75,85], 
                  crs = ccrs.PlateCarree())
    ax.coastlines('10m', 
                  color = 'gray', 
                  lw = 0.5)
    ax.set_title(datetime.strftime(pd.to_datetime(new_polynya_delta.time[day].item()), 
                                   format = '%Y%m%d') + ' - ' + datetime.strftime(pd.to_datetime(new_polynya_delta.time[day - 1].item()), 
                                                                                  format = '%Y%m%d'))
plt.show()

fig = plt.figure()
ax = fig.add_subplot(projection = ccrs.Orthographic(central_longitude=0, central_latitude=80))
new_polynya_delta.std(['time']).plot(ax = ax, transform = ccrs.PlateCarree(), cmap = 'afmhot_r', cbar_kwargs={"label": "Standard Dev. from 2023.05.15-2023.06.01"})
ax.coastlines('10m', color = 'gray', lw = 0.5)
ax.set_extent([-20,20,75,85], crs = ccrs.PlateCarree())

ax.set_title('')
# %%
