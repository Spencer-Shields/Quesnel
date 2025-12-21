"""
This script is for mass-acquiring PlanetScope scene imagery.
It works in two phases. In the first phase, Planet's Data API is queried to get the IDs of scenes which meet certain criteria.
In the second phase, theses IDs are passed to Planet's Orders API, they are ordered, and then they are downloaded.
Multiple orders are automatically created and placed if you want more than 500 scenes.

Search for the tag TO_EDIT to find places where you may want to modify the script.
"""
#%% import libraries, set API
import json
import os
import pathlib
import time
import pandas as pd
import requests
import shapely.geometry
from requests.auth import HTTPBasicAuth
from planet import Session, DataClient, OrdersClient

#TO_EDIT: input your API key
# if your Planet API Key is not set as an environment variable, you can paste it below
if os.environ.get('PL_API_KEY', ''):
    API_KEY = os.environ.get('PL_API_KEY', '')
else:
    API_KEY = 'PLAK9e79c41de89d4189b6c15ae280dd2c6b'

#setup session
session = requests.Session()

#Authenticate
session.auth = (API_KEY, '')


#%% use Data API to get list of scene IDs for download
"""
Reference https://docs.planet.com/develop/apis/data/item-search/#andfilter
"""
#URL for querying data API
data_url = 'https://api.planet.com/data/v1'

#URL for returning search results from data API query
quick_search_url = '{}/quick-search'.format(data_url)

#define function to print formatted json
def p(data):
    print(json.dumps(data, indent=2))

#date filter
start_date = '2013-01-01'
end_date = '2024-12-31'

date_filter = {
    'type': 'DateRangeFilter',
    'field_name': 'acquired',
    'config': {
        'gte': f'{start_date}T00:00:00.000Z', #start date
        'lte': f'{end_date}T23:59:59.999Z' #end date
    }
}

#Geometry filter
Quesnel_bbox = {'type': 'Polygon',
 'coordinates': [[[-122.78314758004296, 53.33052904360998],
   [-122.70020261459693, 53.3282536221804],
   [-122.69657565349696, 53.37451376544108],
   [-122.779611532762, 53.37679163403286],
   [-122.78314758004296, 53.33052904360998]]]}

geometry_filter = {
    'type': 'GeometryFilter',
    'field_name': 'geometry',
    'relation': 'intersects',
    'config': Quesnel_bbox
}


# instrument_filter = {
#     'type': 'StringInFilter',
#     'field_name': 'instrument',
#     'config': [
#         'PS2',  #dove classic (4 band)
#         'PS2.SD', #dove-r (4 band)
#         'PSB.SD' #super dove (4 band or 8 band)
#     ]
# }

#cloud filter
max_cloud_cover = 0.99 #max proportion cloud cover (1 = 100%)

cloud_filter = {
    'type': 'RangeFilter',
    'field_name': 'cloud_cover',
    'config': {
        'lte': 0.99 #less than or equal to 99% cloud cover
        ,'gte': 0 #greater than or equal to 0% cloud cover
    }
}

#quality filter
quality_filter = {
    'type': 'StringInFilter',
    'field_name': 'quality_category',
    'config': ['standard'] #only standard quality
}

#asset filter (so that I'm only getting scenes where surface reflectance and udm2 assets are available)
# asset_filter = {
#   "type": "AndFilter",
#   "config": [
#     {
#       "type": "AssetFilter",
#       "config": ["analytic_sr"]
#     },
#     {
#       "type": "AssetFilter",
#       "config": ["udm2"]
#     }
#   ]
# }

# asset_filter = {
#   "type": "AndFilter",
#   "config": [
#     {
#       "type": "AssetFilter",
#       "config": ["ortho_analytic_4b_sr"]
#     },
#     {
#       "type": "AssetFilter",
#       "config": ["udm2"]
#     }
#   ]
# }

asset_filter = {
  "type": "AndFilter",
  "config": [
    {
      "type": "AssetFilter",
      "config": ["analytic_sr"]      # analytic SR asset
    }
    # ,
    # {
    #   "type": "AssetFilter",
    #   "config": ["ortho_udm2"]       # UDM2 asset for PSScene
    # }
  ]
}

#permission filter (only get scenes that can be downloaded)
permission_filter ={
      "type": "PermissionFilter",
      "config": [
        #   "assets:download"
        "assets.ortho_analytic_4b_sr:download",
        "assets.ortho_udm2:download"
          ]
    }
  

#TO_EDIT: double check filters before searching for Planet scenes
#combine filters
combined_filter = {
    'type': 'AndFilter',
    'config': [
        date_filter,
        geometry_filter,
        cloud_filter,
        quality_filter,
        # asset_filter,
        permission_filter
        # instrument_filter
    ]
}

#set up request
item_types = ['PSScene'] #look for planetscope scenes (as opposed to RapidEye, SkySat, etc)

data_request = {
    'item_types': item_types,
    'filter': combined_filter
}

#send request to data API quick search endpoint
data_response = session.post(quick_search_url, json=data_request)

#assign response to a variable
data_response_geojson = data_response.json()

#paginate through response to aggregate all results
links = data_response_geojson['_links']
feats = data_response_geojson['features']

all_features = feats
next_link = data_response_geojson['_links']['_next']

while next_link:
    next_response = session.get(next_link)
    next_response_geojson = next_response.json()
    next_feats = next_response_geojson['features']
    all_features.extend(next_feats)
    next_link = next_response_geojson['_links'].get('_next', None)

#TO_EDIT: do additional filtering using dataframe if desired (just make sure that the final dataframe is called all_features_df)
#put all results in a dataframe
all_features_df = pd.DataFrame([feat['properties'] for feat in all_features]) #make dataframe
all_features_df['id'] = [feat['id'] for feat in all_features] #add column for scene ids

#make dataframe column saying what proportion of the aoi is covered by each scene
def compute_aoi_proportion(feature, aoi_geom):
    feature_geom = feature['geometry']
    feature_shape = shapely.geometry.shape(feature_geom)
    aoi_shape = shapely.geometry.shape(aoi_geom)
    intersection = feature_shape.intersection(aoi_shape)
    if intersection.is_empty:
        return 0.0
    else:
        return intersection.area / aoi_shape.area
all_features_df['proportion_aoi_covered'] = [
    compute_aoi_proportion(feat, Quesnel_bbox) for feat in all_features
    ]

# %% set up to place order

#extract ids of scenes to order from features dataframe
desired_feat_ids = all_features_df['id'].tolist()

#remove duplicate ids
desired_feat_ids = list(set(desired_feat_ids))

#set url for interacting with orders API
orders_url = 'https://api.planet.com/compute/ops/orders/v2'

#reorder the strings so that similar dates are grouped together
desired_feat_ids.sort()

#organize ids into chunks due to order size limits
chunk_size = 500
desired_feat_chunks = [
    desired_feat_ids[i:i+chunk_size] 
    for i in range(0, len(desired_feat_ids), chunk_size)
    ]

#check authentication (want response code 200)
order_response = requests.get(orders_url, auth=session.auth)
order_response

#set content type 
headers = {'content-type': 'application/json'}

#make poll for success function to monitor order status
def poll_for_success(order_url, auth, num_loops=30, wait_sec = 30):
    count = 0
    while(count < num_loops):
        count += 1
        r = requests.get(order_url, auth=auth)
        response = r.json()
        state = response['state']
        print(state)
        end_states = ['success', 'failed', 'partial']
        if state in end_states:
            break
        time.sleep(wait_sec)

#make function for downloading results
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

#TO_EDIT: check product bundle selection (analytic_sr_udm2 is 4band, use analytic_8b_sr_udm2 for 8band)
#define product bundle to order
product_bundle = 'analytic_sr_udm2' #4band surface reflectance with UDM2

#TO_EDIT: make basename for the order (chunk number will be added automatically)
#make base name for ordering all the chunks
order_basename = f'Quesnel_Clipped_Harmonized_Coregistered_Product={product_bundle}_Time={start_date}to{end_date}_MaxCloudProp={max_cloud_cover}_MaxChunkSize={chunk_size}'

#%%place orders for each chunk or get urls if orders have already been placed

#get names and ids of previous orders, place in dataframe
orders_response = requests.get(orders_url, auth=session.auth)
previous_orders = order_response.json()['orders']
previous_order_names = [r['name'] for r in previous_orders]
previous_order_ids = [r['id'] for r in previous_orders]
previous_orders_df = pd.DataFrame({
    'name': previous_order_names,
    'id': previous_order_ids})
previous_orders_df = previous_orders_df.assign(url = lambda x: orders_url + '/' + x['id'])

#initiate list to store order urls
order_urls_list = []

#loop over chunks of ids, place order if not already placed, otherwise get order url
for i in range(len(desired_feat_chunks)):

    #make order name
    order_name = order_basename + f'_chunk_{i+1}_of_{len(desired_feat_chunks)}'

    #proceed if order name is not in previous orders
    if order_name not in previous_orders_df['name'].tolist():
        print('Placing order for chunk' + str(i + 1) + ' of ' + str(len(desired_feat_chunks)))

        #get current chunk of scene ids
        chunk = desired_feat_chunks[i]

        ##select scene to use for coregistration (NOTE: THIS WILL NEED TO BE TWEAKED IF THE AOI IS LARGER THAN A SINGLE SCENE)
        #filter dataframe to get scenes which fully cover aoi, with high visible percent, and with high visible confidence percent
        chunk_features_df = all_features_df[all_features_df['id'].isin(chunk)]
        max_aoi_coverage = chunk_features_df['proportion_aoi_covered'].max()
        candidate_scenes_df = chunk_features_df[(chunk_features_df['proportion_aoi_covered'] == max_aoi_coverage)]
        candidate_scenes_df = candidate_scenes_df[(candidate_scenes_df['visible_percent'] == max(candidate_scenes_df['visible_percent']))]
        candidate_scenes_df = candidate_scenes_df[(candidate_scenes_df['visible_confidence_percent'] == max(candidate_scenes_df['visible_confidence_percent']))]

        #get id of first scene in candidate scenes dataframe to use as coregistration anchor
        coreg_anchor_id = candidate_scenes_df.iloc[0]['id']




        #set up order request
        order_request = {
            "name": order_name,
            'order_type': 'partial', #ignore scene ids that cannot be downloaded
            "products": [
                {
                    "item_ids": chunk,
                    "item_type": "PSScene",
                    "product_bundle": product_bundle #get 4band imagery
                }
            ],
            "delivery": {
                "archive_type": "zip",
                "single_archive": True
            },
            "tools": [ #TO_EDIT: add preprocessing tools
                {
                    "clip": {
                        "aoi": Quesnel_bbox  #clip to polygon or multipolygon geojson
                            }
                },
                {
                    "coregister": { #coregister all images to reference item
                        "anchor_item": coreg_anchor_id
                                }
                          },
              {
                    "harmonize": { #radiometric harmonization with Sentinel-2
                          "target_sensor": "Sentinel-2"
                                }
                          }
            ]
}

        #place order
        order_response = requests.post(
            orders_url,
            auth =session.auth,
            data=json.dumps(order_request),
            headers=headers
        )
        print(order_response.json())
        order_id = order_response.json()['id']
        print(order_id)
        order_url = orders_url + '/' + order_id

        order_urls_list.append(order_url)
        countdown_timer(1) #wait 1 second between orders to avoid rate limit errors

    else:
        print('Order exists for chunk' + str(i + 1) + ' of ' + str(len(desired_feat_chunks)))
        order_url = previous_orders_df.loc[previous_orders_df['name'] == order_name, 'url'].values[0]
        order_urls_list.append(order_url)

#%% Download orders using order urls

#TO_EDIT: set ouptut directory for downloading results
#define output directory to download results to
download_dir = pathlib.Path(rf'D:\Quesnel\data\planet_scenes\raw_additional\{order_basename}')
if not download_dir.exists(): #create directory if it doesn't exist
    download_dir.mkdir(parents=True, exist_ok=True)

#loop over order urls, download to download directory
for i in range(0, len(order_urls_list)):

    order_url = order_urls_list[i]

    print('Processing order: ' + order_url)

    #poll for success
    poll_for_success(order_url, session.auth, num_loops=300, wait_sec=60)

    #get order results
    r = requests.get(order_url, auth=session.auth)
    response = r.json()
    results = response['_links']['results']

    #download results
    download_results(results, download_dir=download_dir)




