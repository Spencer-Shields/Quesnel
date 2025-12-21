
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
start_date = '2013-01-01T00:00:00.000Z'
end_date = '2025-01-01T00:00:00.000Z'


date_filter = {
    'type': 'DateRangeFilter',
    'field_name': 'acquired',
    'config': {
        'lt': end_date,
        'gte': start_date
    }
}

#instrument filter
instrument_filter = {
    'type': 'StringInFilter',
    'field_name': 'instrument',
    'config': [
        'PS2.SD', 
        'PS2', 
        'PSB.SD']
    }

# #geometry filter
# aoi_path = r'Quesnel_thinning\AOI_fullsite_wgs84.geojson' #load geojson of site
# with open(aoi_path, 'r') as f:
#     aoi_json = json.load(f)
# aoi_geom = aoi_json['features'][0]['geometry']
# p(aoi_json) #print to verify that file loaded successfully
aoi_geom = {'type': 'Polygon',
 'coordinates': [[[-122.78314758004296, 53.33052904360998],
   [-122.70020261459693, 53.3282536221804],
   [-122.69657565349696, 53.37451376544108],
   [-122.779611532762, 53.37679163403286],
   [-122.78314758004296, 53.33052904360998]]]}

geometry_filter = {
    'type' : 'GeometryFilter',
    'field_name': 'geometry',
    'config': aoi_geom
}

#cloud filter
min_cloud_proportion = 0.0
max_cloud_proportion = 100.0

cloud_filter = {
    'type': 'RangeFilter',
    'field_name': 'cloud_cover',
    'config': {
        'gte': min_cloud_proportion,
        'lte': max_cloud_proportion
        }
    }

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


# %% Extract metadata from features, use to filter desired scenes

#make simple feature object using PSScenes metadata
df = pd.DataFrame([feat['properties'] for feat in all_features]) #make dataframe of metadata
df["geometry"] = [shape(feat["geometry"]) for feat in all_features]
gdf = gpd.GeoDataFrame(df, geometry="geometry")
gdf['feature_id'] = [feat['id'] for feat in all_features]

gdf.set_crs(epsg=4326, inplace=True) #set crs to wgs84


# %% Make request

#get list of asset ids that I want
desired_feat_ids = gdf['feature_id'].tolist()
len(desired_feat_ids)

#remove duplicate ids
desired_feat_ids = list(set(desired_feat_ids))
len(desired_feat_ids)


#split desired scenes into chunks (due to order size limits)
chunk_size = 100
desired_feat_chunks = [desired_feat_ids[i:i+chunk_size] for i in range(0, len(desired_feat_ids), chunk_size)]

for x in desired_feat_chunks:
    print(len(x))

#place order
orders_url = 'https://api.planet.com/compute/ops/orders/v2' 
response = requests.get(orders_url, auth=session.auth)
response

headers = {'content-type': 'application/json'}

#%% Define functions for placing order, polling for success, downloading results

#define function for placing an order
def place_order(request, auth):
    response = requests.post(orders_url, data=json.dumps(request), auth=auth, headers=headers)
    print(response)
    order_id = response.json()['id']
    print(order_id)
    order_url = orders_url + '/' + order_id
    return order_url

# ##place the order (uncomment to run)
# orders_url = place_order(request, session.auth)

#check status of order
def poll_for_success(order_url, auth, num_loops=1000):
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
        time.sleep(20)

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

#define a timer function (to help prevent throwing a 405 error when placing orders too quickly)
def countdown_timer(seconds):
    while seconds:
        mins, secs = divmod(seconds, 60)
        time_format = f'Wait to order next chunk: {mins:02}:{secs:02}'
        print(time_format, end='\r')
        time.sleep(1)
        seconds -= 1
    print("00:00\nTime's up!")

#%% extra functions for order checking

def list_orders(orders_url, auth):
    response = requests.get(orders_url, auth=auth)
    response.raise_for_status()
    return response.json().get('orders', [])

def find_existing_order(orders, order_name):
    for order in orders:
        if order.get('name') == order_name:
            return order
    return None

def get_order_status(order_url, auth):
    response = requests.get(order_url, auth=auth)
    response.raise_for_status()
    return response.json()

#%% Iterate over chunks placing order and downloading results

GET_SCENES = True #set to true to order and download scenes. Set to false or comment to skip

order_name = 'Quesnel_2013to2025'

#define path of directory where imagery will be downloaded to
download_dir = os.path.join(wd_path, 'planet_scenes', 'raw_additional',order_name)

if GET_SCENES:
    print('Placing order for ' + order_name)

    download_dir = os.path.join(wd_path, 'planet_scenes', 'raw_additional', order_name)
    if os.path.exists(download_dir):
        print('Download directory exists, check if order is already placed')
    else:
        os.makedirs(download_dir)

    for i, chunk in enumerate(desired_feat_chunks, start=0):
        print(f'Placing order for chunk {i+1} of {len(desired_feat_chunks)}')

        request = {
            "name": order_name + "_chunk_" + str(i + 1),
            # "source_type": "scenes",
            "products": [
                {
                    "item_ids": chunk,
                    "item_type": "PSScene",
                    "product_bundle": "analytic_udm2" #get 4band imagery
                }
            ],
            # "delivery": {
            # "archive_type": "zip",
            # "single_archive": True,
            # "archive_filename": "{{name}}_{{order_id}}.zip"
            # }
            # 'tools': [
            #     {
            #         "coregister": {
            #             "anchor_item": "20210726_191611_86_2307"
            #         }
            #     }
            # ]
        }

        orders_url = place_order(request, session.auth)
        poll_for_success(orders_url, session.auth)
        r = requests.get(orders_url, auth=session.auth)
        response = r.json()
        print(response)
        results = response['_links']['results']
        [r['name'] for r in results]
        print('Downloading results for chunk' + str(i + 1) + ' of ' + str(len(desired_feat_chunks)))
        download_results(results, download_dir)

        countdown_timer(30) #wait x seconds before placing next order
# %%
