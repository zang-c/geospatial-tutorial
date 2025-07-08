# -*- coding: utf-8 -*-
"""
Created on 2025.06.25

I want to merge all of my GeoTiff files to create a single layered GeoTiff file where the layers
communicate the time series information.

@author: Willis Lab
"""
#%%
import numpy as np
import rioxarray
from rioxarray.rioxarray import affine_to_coords
import xarray
import rasterio
import rasterio.plot
import sys, os, time
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
import matplotlib.patheffects
import matplotlib.colors as mcolors
from matplotlib import rcParams
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy as cart
import cartopy.io.img_tiles as cimgt
from cartopy.io.img_tiles import GoogleTiles
import geopy.distance
from cmcrameri import cm as cmc
import io
from urllib.request import urlopen, Request
from PIL import Image
import contextily as cx

def add_map_features(ax):
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(provinc_bodr)
    ax.add_feature(cfeature.LAND, color = "#fff9da")

resol = '10m'

provinc_bodr = cart.feature.NaturalEarthFeature(category='cultural', 
                                                name='admin_1_states_provinces_lines', 
                                                scale=resol, 
                                                facecolor='none', 
                                                edgecolor='k')
urban_area = cart.feature.NaturalEarthFeature(category='cultural', 
                                         name = 'urban_areas', 
                                         scale = resol, 
                                         facecolor = "#FCA729",
                                         edgecolor = 'k')
rivers = cart.feature.NaturalEarthFeature(category='physical',
                                          name = 'rivers_north_america',
                                          scale = resol,
                                          facecolor = "none",
                                          edgecolor = "#93b7d4",
                                          ls = ':')
topo = cart.feature.NaturalEarthFeature(category='physical',
                                          name = 'land_ocean_seams',
                                          scale = resol)

# lon_start = input('What are the longitudinal coordinates of the city/town/area where you grew up?')
# lat_start = input('What are the latitudinal coordninates of the city/town/area where you grew up?')
lon_start = -77.6997
lat_start = 39.4862
lon_start = float(lon_start) 
lat_start = float(lat_start)
fig = plt.figure()

imagery = cimgt.OSM()
#tiler = TracestrackTopo()
ax = fig.add_subplot(projection = ccrs.Mercator())
#zoom = 11
#ax.add_image(imagery, zoom)
#ax.add_image(tiler, zoom)
#add_map_features(ax)
ax.background_img(name='ETOPO', resolution='high')
ax.scatter(lon_start, lat_start, transform = ccrs.PlateCarree())
#ax.set_extent([lon_start - 1, lon_start + 1, lat_start - 0.5, lat_start + 0.5])
# %%

#!wget http://data.nodc.noaa.gov/thredds/fileServer/woa/WOA09/NetCDFdata/temperature_annual_1deg.nc