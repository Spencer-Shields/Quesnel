
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
start_date = '2020-01-01T00:00:00.000Z'
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
        # 'PS2.SD', 
        # 'PS2', 
        'PSB.SD']
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



# %% Extract metadata from features, use to filter desired scenes

#make simple feature object using PSScenes metadata
df = pd.DataFrame([feat['properties'] for feat in all_features]) #make dataframe of metadata
df["geometry"] = [shape(feat["geometry"]) for feat in all_features]
gdf = gpd.GeoDataFrame(df, geometry="geometry")
gdf['feature_id'] = [feat['id'] for feat in all_features]

gdf.set_crs(epsg=4326, inplace=True) #set crs to wgs84

#number of scenes for each individual satellite (no useful trends)
n_sat = df['satellite_id'].value_counts()
print(n_sat)


#extract month and year from 'acquired' column, get count of scenes for each month in timeseries
gdf["acquired"] = pd.to_datetime(gdf["acquired"])
gdf['year_month'] = gdf['acquired'].dt.strftime('%Y-%m')  # Extract YYYY-MM
n_tot = gdf['year_month'].value_counts().sort_index()
print(n_tot)


#filter by percent cloud cover
max_cloud = 0.1
gdf = gdf[gdf['cloud_cover'] < max_cloud]

#only return images with ground control
gdf = gdf[gdf['ground_control'] == True]

#get scenes which completely contain the AOI
aoi_shapely = shape(aoi_json['features'][0]['geometry'])#convert aoi to shapely object


gdf["aoi_intersection"] = gdf.geometry.intersection(aoi_shapely) #perform intersection
gdf["aoi_proportion"] = gdf["aoi_intersection"].area / aoi_shapely.area #get intersection area divided by aoi area

# aoi_prop_threshold = 0.5
# gdf_coverage = gdf[gdf['aoi_proportion'] > aoi_prop_threshold] #get scenes which completely contain the AOI


####get scenes which completely cover individual thinning plots

block_path = r'Quesnel_thinning\12l_12n_bdy.geojson' #load thinning plot data as simple feature object
with open(block_path, 'r') as f:
    block_json = json.load(f)

block_gdf = gpd.GeoDataFrame.from_features(block_json['features']) #convert block_json dict object to geodataframe
block_gdf.set_crs(epsg=block_json['crs']['properties']['name'].split('::')[1]
                  , inplace=True) #set crs

if gdf.crs != block_gdf.crs: #check if crs of gdf and block_gdf match, convert crs of block_gdf to match gdf if not
    block_gdf = block_gdf.to_crs(gdf.crs)


from geopandas.tools import sjoin
containing_gdf = gpd.sjoin(block_gdf, gdf, how="inner", predicate='within') # Perform a spatial join to find which features in gdf contain the block geometry
gdf_coverage = containing_gdf.drop_duplicates(subset=['feature_id']) #filter duplicate rows from containing_gdf

#get year-month of scenes in gdf_fullcoverage based on 'acquired' column
gdf_coverage = gdf_coverage.copy()
gdf_coverage["acquired"] = pd.to_datetime(gdf_coverage["acquired"])
gdf_coverage['year_month'] = gdf_coverage['acquired'].dt.strftime('%Y-%m')  # Extract YYYY-MM
gdf_coverage['month'] = gdf_coverage['acquired'].dt.strftime('%m')  # Extract MM

year_month_counts = gdf_coverage['year_month'].value_counts().sort_index()
print(year_month_counts)

# #plot number of scenes per month
# import matplotlib.pyplot as plt

# # Plot the number of scenes per month
# plt.figure(figsize=(10, 6))
# year_month_counts.plot(kind='bar', color='skyblue', edgecolor='black')
# plt.title('Number of Scenes per Month', fontsize=16)
# plt.xlabel('Year-Month', fontsize=14)
# plt.ylabel('Number of Scenes', fontsize=14)
# plt.xticks(rotation=45, fontsize=10)
# plt.grid(axis='y', linestyle='--', alpha=0.7)
# plt.tight_layout()
# plt.show()

# #limit to only imagery from certain months (i.e. growing season)
# exclude_months = ['12', '01', '02', '03', '04']
exclude_months = []
gdf_gs = gdf_coverage[~gdf_coverage['month'].isin(exclude_months)]

len(gdf_gs)

n_gs = gdf_gs['year_month'].value_counts().sort_index()
n_gs
n_gs.sum()


# %% Make request, place order

# #get list of assets ids that I already have
# already_have_path = r"D:\Quesnel\data\planet_scenes\raw\a2c2dc75-8eab-428f-894e-6ff0baaea934\manifest.json"
# with open(already_have_path, 'r') as f:
#     already_have_json = json.load(f)
# already_have_ids = [f['annotations']['planet/item_id'] for f in already_have_json['files']]

#get list of asset ids that I want
desired_feat_ids = gdf_gs['feature_id'].tolist()
# desired_feat_ids = [f for f in desired_feat_ids if f not in already_have_ids] #remove already downloaded scenes from list of desired scenes
len(desired_feat_ids)

#place order
orders_url = 'https://api.planet.com/compute/ops/orders/v2' 
response = requests.get(orders_url, auth=session.auth)
response

headers = {'content-type': 'application/json'}

request = {  
   "name":order_name + "_chunk_2",
   "products":[
      {  
         "item_ids": desired_feat_chunks[1],
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

poll_for_success(orders_url, session.auth)

#define path of directory where imagery will be downloaded to
download_dir = os.path.join(wd_path, 'planet_scenes', 'raw')

#see results
r = requests.get(orders_url, auth=session.auth)
response = r.json()
results = response['_links']['results']
[r['name'] for r in results]

#%% Download imagery

#define path of directory where imagery will be downloaded to
download_dir = os.path.join(wd_path, 'planet_scenes', 'raw',order_name)

download_results(results, download_dir)

#define a time function (to help prevent throwing a 405 error when placing orders too quickly)
def countdown_timer(seconds):
    while seconds:
        mins, secs = divmod(seconds, 60)
        time_format = f'Wait to order next chunk: {mins:02}:{secs:02}'
        print(time_format, end='\r')
        time.sleep(1)
        seconds -= 1
    print("00:00\nTime's up!")

## Example: countdown from 10 seconds
#countdown_timer(10)

#split desired scenes into chunks of 500 (max order size)
desired_feat_chunks = [desired_feat_ids[i:i+500] for i in range(0, len(desired_feat_ids), 500)]

#iterate over chunks placing order and downloading results

GET_SCENES = True #set to true to order and download scenes. Set to false or comment to skip

order_name = 'Quesnel_8b_scenes_harmonized_clipped_coveragebyblock_0.1cloud_2021to2024'

if GET_SCENES:
    print('Placing order for ' + order_name)

    download_dir = os.path.join(wd_path, 'planet_scenes', 'raw', order_name)
    if os.path.exists(download_dir):
        print('Download directory exists, check if order is already placed')
    else:
        os.makedirs(download_dir)

        for i, chunk in enumerate(desired_feat_chunks):
            print(f'Placing order for chunk {i+1} of {len(desired_feat_chunks)}')

            request = {
                "name": order_name + "_chunk_" + str(i + 1),
                "products": [
                    {
                        "item_ids": chunk,
                        "item_type": "PSScene",
                        "product_bundle": "analytic_8b_sr_udm2"
                    }
                ],
                'tools': [
                    {
                        "reproject": {
                            "projection": "WGS84",
                            "kernel": "cubic"
                        }
                    },
                    {
                        'clip': {
                            'aoi': aoi_geom
                        }
                    },
                    {
                        'harmonize': {
                            'target_sensor': 'Sentinel-2'
                        }
                    }
                ]
            }

            orders_url = place_order(request, session.auth)
            poll_for_success(orders_url, session.auth)
            r = requests.get(orders_url, auth=session.auth)
            response = r.json()
            results = response['_links']['results']
            [r['name'] for r in results]
            print('Downloading results for chunk' + str(i + 1) + ' of ' + str(len(desired_feat_chunks)))
            download_results(results, download_dir)

            countdown_timer(300) #wait 5 minutes before placing next order