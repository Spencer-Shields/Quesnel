#%% Import packages, set working directory

import os
import json
import requests
import geojsonio
import time
import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib as plt
from shapely.geometry import shape

#set working directory to project folder
wd_path = r"D:\Quesnel\data"
os.chdir(wd_path)

#%% Configure session for downloading data via API

# Helper function to printformatted JSON using the json module
def p(data):
    print(json.dumps(data, indent=2))

# if your Planet API Key is not set as an environment variable, you can paste it below
if os.environ.get('PL_API_KEY', ''):
    API_KEY = os.environ.get('PL_API_KEY', '')
else:
    API_KEY = 'PLAK9e79c41de89d4189b6c15ae280dd2c6b'

# construct auth tuple for use in the requests library
BASIC_AUTH = (API_KEY, '')

# Setup Planet Data API base URL
URL = "https://api.planet.com/data/v1"

# Setup the session
session = requests.Session()

# Authenticate
session.auth = (API_KEY, "")

# Specify the sensors/satellites or "item types" to include in our results
item_types = ["PSScene"]

# # Make a GET request to the Planet Data API

res = session.get(URL)

# Print formatted JSON response
p(res.json())

# Print the value of the item-types key from _links
print(res.json()["_links"]["item-types"])

# Setup the stats URL (provides summary of available data based on a search)
stats_url = "{}/stats".format(URL)

# Print the stats URL
print(stats_url)


#%% Configure filter for downloading imagery

# Specify the sensors/satellites or "item types" to include in our results
item_types = ["PSScene"]

# date filter
start_date = '2021-06-01T00:00:00.000Z'
end_date = '2024-11-01T00:00:00.000Z'


date_filter = {
    'type': 'DateRangeFilter',
    'field_name': 'acquired',
    'config': {
        'lt': end_date,
        'gte': start_date
    }
}

#geometry filter
aoi_path = r'Quesnel_thinning\AOI_fullsite_wgs84.geojson' #load geojson of site
with open(aoi_path, 'r') as f:
    aoi_json = json.load(f)
aoi_geom = aoi_json['features'][0]['geometry']
p(aoi_json) #print to verify that file loaded successfully

geometry_filter = {
    'type' : 'GeometryFilter',
    'field_name': 'geometry',
    'config': aoi_geom
}

#instrument filter
instrument_filter = {
    'type': 'StringInFilter',
    'field_name': 'instrument',
    'config': ['PS2.SD', 'PS2', 'PSB.SD']
    }

#cloud filter
min_cloud_proportion = 0
max_cloud_proportion = 0.05

cloud_filter = {
    'type': 'RangeFilter',
    'field_name': 'cloud_cover',
    'config': {
        'gte': min_cloud_proportion,
        'lte': max_cloud_proportion
        }
    }

#view angle filter

# view_angle_filter = {
#     'type': 'RangeFilter',
#     'field_name': 'view_angle',
#     'config': {
#         'gte': 0,
#         'lte': 10
#         }
#     }

#combine all filters
and_filter = {
    'type' : 'AndFilter',
    'config': [
        date_filter
        , geometry_filter
        # , instrument_filter
        , cloud_filter
        ]
    }

#%% Construct request
scenes_request = {
    "item_types" : item_types,
    "interval" : 'month',
    'filter' : and_filter
    }

#send the POST request to the API stats endpoint
res_scenes = session.post(stats_url, json=scenes_request)
p(res_scenes.json()) #see number of results in interval buckets

# %% Set up quick search
quick_url = "{}/quick-search".format(URL)

request = {
    'item_types' : item_types,
    'filter' : and_filter
}

res = session.post(quick_url, json=request)

gj = res.json()
p(gj)

links = gj['_links']
feats = gj['features']
p(feats)
print(type(feats))

#make simple feature object using PSScenes metadata
df = pd.DataFrame([feat['properties'] for feat in feats]) #make dataframe of metadata
df["geometry"] = [shape(feat["geometry"]) for feat in feats]
gdf = gpd.GeoDataFrame(df, geometry="geometry")


#number of scenes for each individual satellite (no useful trends)
n_sat = df['satellite_id'].value_counts()
print(n_sat)

#extract month and year from 'acquired' column, get count of scenes for each month in timeseries
gdf["acquired"] = pd.to_datetime(gdf["acquired"])
gdf['year_month'] = gdf['acquired'].dt.strftime('%Y-%m')  # Extract YYYY-MM
n_tot = gdf['year_month'].value_counts().sort_index()
print(n_tot)

#get scenes which completely contain the AOI
aoi_shapely = shape(aoi_json['features'][0]['geometry'])#convert aoi to shapely object


gdf["aoi_intersection"] = gdf.geometry.intersection(aoi_shapely) #perform intersection
gdf["aoi_proportion"] = gdf["aoi_intersection"].area / aoi_shapely.area #get intersection area divided by aoi area

aoi_prop_threshold = 0
gdf_fullcoverage = gdf[gdf['aoi_proportion'] > aoi_prop_threshold] #get scenes which completely contain the AOI

#get year-month of scenes in gdf_fullcoverage based on 'acquired' column
gdf_fullcoverage = gdf_fullcoverage.copy()
gdf_fullcoverage["acquired"] = pd.to_datetime(gdf_fullcoverage["acquired"])
gdf_fullcoverage['year_month'] = gdf_fullcoverage['acquired'].dt.strftime('%Y-%m')  # Extract YYYY-MM
gdf_fullcoverage['month'] = gdf_fullcoverage['acquired'].dt.strftime('%m')  # Extract MM

year_month_counts = gdf_fullcoverage['year_month'].value_counts().sort_index()

#limit to only imagery from march to november (i.e. growing season)
exclude_months = ['12', '01', '02', '03', '04']
gdf_gs = gdf_fullcoverage[~gdf_fullcoverage['month'].isin(exclude_months)]

n_gs = gdf_gs['year_month'].value_counts().sort_index()
len(n_gs)

# %%
