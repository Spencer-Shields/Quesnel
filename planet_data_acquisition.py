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
import planet
import pathlib

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
    'config': [
        # 'PS2.SD', 
        # 'PS2', 
        'PSB.SD']
    }

#cloud filter
min_cloud_proportion = 0.0
max_cloud_proportion = 0.0

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
        , instrument_filter
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

#%% Aggregate results
all_features = feats
next_link = gj['_links']['_next']

while next_link:
    res = session.get(next_link)
    gj = res.json()
    feats = gj['features']
    all_features.extend(feats)
    next_link = gj['_links'].get('_next', None)



# %% Extract metadata from features
#make simple feature object using PSScenes metadata
df = pd.DataFrame([feat['properties'] for feat in all_features]) #make dataframe of metadata
df["geometry"] = [shape(feat["geometry"]) for feat in all_features]
gdf = gpd.GeoDataFrame(df, geometry="geometry")
gdf['feature_id'] = [feat['id'] for feat in all_features]


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

aoi_prop_threshold = 0.999
gdf_coverage = gdf[gdf['aoi_proportion'] > aoi_prop_threshold] #get scenes which completely contain the AOI

#get year-month of scenes in gdf_fullcoverage based on 'acquired' column
gdf_coverage = gdf_coverage.copy()
gdf_coverage["acquired"] = pd.to_datetime(gdf_coverage["acquired"])
gdf_coverage['year_month'] = gdf_coverage['acquired'].dt.strftime('%Y-%m')  # Extract YYYY-MM
gdf_coverage['month'] = gdf_coverage['acquired'].dt.strftime('%m')  # Extract MM

year_month_counts = gdf_coverage['year_month'].value_counts().sort_index()

# #limit to only imagery from march to november (i.e. growing season)
# exclude_months = ['12', '01', '02', '03', '04']
exclude_months = []
gdf_gs = gdf_coverage[~gdf_coverage['month'].isin(exclude_months)]

len(gdf_gs)

n_gs = gdf_gs['year_month'].value_counts().sort_index()
n_gs


# %% Define order

#get list of asset ids that I want

desired_feat_ids = gdf_gs['feature_id']

orders_url = 'https://api.planet.com/compute/ops/orders/v2' 
response = requests.get(orders_url, auth=session.auth)
response

headers = {'content-type': 'application/json'}

request = {  
   "name":"Quesnel_8b_scenes_harmonized_clipped_fullcoverage_2021to2024",
   "products":[
      {  
         "item_ids": desired_feat_ids.tolist(),
         "item_type":"PSScene",
         "product_bundle":"analytic_8b_sr_udm2" #get analytic 8b
      }
   ],
    'tools': [
        {
            "reproject": {
                "projection": "WGS84",
                "kernel": "cubic"}
                },
        {
            'clip': {
                'aoi': aoi_geom
            }
        },
        {
        'harmonize': {
            'target_sensor' : 'Sentinel-2'
        }
    }
    ]
}

#%% Place order

#define function for placing an order
def place_order(request, auth):
    response = requests.post(orders_url, data=json.dumps(request), auth=auth, headers=headers)
    print(response)
    order_id = response.json()['id']
    print(order_id)
    order_url = orders_url + '/' + order_id
    return order_url

#place the order (uncomment to run)
orders_url = place_order(request, session.auth)

#check status of order
def poll_for_success(order_url, auth, num_loops=150):
    count = 0
    while(count < num_loops):
        count += 1
        r = requests.get(order_url, auth=session.auth)
        response = r.json()
        state = response['state']
        print(state)
        end_states = ['success', 'failed', 'partial']
        if state in end_states:
            break
        time.sleep(10)
        
poll_for_success(orders_url, session.auth)

#see results
r = requests.get(orders_url, auth=session.auth)
response = r.json()
results = response['_links']['results']
[r['name'] for r in results]

#%% Download imagery

#define path of directory where imagery will be downloaded to
download_dir = os.path.join(wd_path, 'planet_scenes', 'raw')

#define function to loop through files and download them
def download_results(results, download_dir='data', overwrite=False):
    results_urls = [r['location'] for r in results]
    results_names = [r['name'] for r in results]
    print('{} items to download'.format(len(results_urls)))
    
    for url, name in zip(results_urls, results_names):
        path = pathlib.Path(os.path.join(download_dir, name))
        
        if overwrite or not path.exists():
            print('Downloading {} to {}'.format(name, path))
            r = requests.get(url, allow_redirects=True)
            path.parent.mkdir(parents=True, exist_ok=True)
            with open(path, 'wb') as f:
                f.write(r.content)
        else:
            print('{} already exists, skipping {}'.format(path, name))

#download the imagery (uncomment to run)
download_results(results, download_dir)