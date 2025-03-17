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
from urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter

#set working directory to project folder
wd = os.getcwd()

#%% Configure session
# if your Planet API Key is not set as an environment variable, you can paste it below
if os.environ.get('PL_API_KEY', ''):
    API_KEY = os.environ.get('PL_API_KEY', '')
else:
    API_KEY = 'PLAK9e79c41de89d4189b6c15ae280dd2c6b'


#%% Configure session for downloading data via API

# Helper function for printing
import re
import pprint

def display(content):
    result = pprint.pformat(content)
    result = re.sub(r'api_key=\w+', 'api_key=XXXX', result)
    print(result)

# if your Planet API Key is not set as an environment variable, you can paste it below
if os.environ.get('PL_API_KEY', ''):
    API_KEY = os.environ.get('PL_API_KEY', '')
else:
    API_KEY = 'PLAK9e79c41de89d4189b6c15ae280dd2c6b'

# construct auth tuple for use in the requests library
BASIC_AUTH = (API_KEY, '')

# Setup Planet Data API base URL
basemaps_api_url = 'https://api.planet.com/basemaps/v1'
url = f'{basemaps_api_url}/mosaics'

# Setup the session
session = requests.Session()

# Authenticate
session.auth = (API_KEY, "")

# We also need to set up requests to honor and retry common "slow down" status codes
# that the API may respond with as well as other retryable statuses.
retries = Retry(total=10, backoff_factor=1, status_forcelist=[429, 502, 503])
session.mount('https://', HTTPAdapter(max_retries=retries))

# # Make a GET request to the Planet Data API, display all results

rv = session.get(url)
rv.raise_for_status()
result = rv.json()
display(result)

#define function for compiling results from multiple pages
def paginated_get(session, url, item_key, **kwargs):
    while True:
        rv = session.get(url, **kwargs)
        rv.raise_for_status()
        page = rv.json()
        
        for item in page[item_key]:
            yield item
            
        if '_next' in page['_links']:
            url = page['_links']['_next']
        else:
            break

#define helper functions for accessing image series
def series_metadata(session, series_name):
    rv = session.get(f'{basemaps_api_url}/series', params={'name__is': series_name})
    rv.raise_for_status()
    return rv.json()['series'][0]

def mosaics_in_series(session, series_name, start_date=None, end_date=None):
    info = series_metadata(session, series_name)
    url = info['_links']['mosaics']
    mosaics = paginated_get(
        session, 
        url, 
        'mosaics', 
        params={
            'acquired__lt': end_date, 
            'acquired__gt': start_date,
        },
    )
    return mosaics

# see all series I have access to
all_series_list = [] 
all_series = paginated_get(session, f'{basemaps_api_url}/series', 'series')
for item in all_series:
    all_series_list.append(item['name'])
    print(item['name'])

all_series_l = list(all_series)
series_names = [item['name'] for item in all_series_l]

#load area of interest
aoi_path = os.path.join(wd,'data', "Quesnel_thinning", "AOI_fullsite_wgs84.geojson") #load geojson of site
with open(aoi_path, 'r') as f:
    aoi_json = json.load(f)
aoi_geom = aoi_json['features'][0]['geometry']
#%% Configure filter for downloading imagery

# define time series
start_date = '2021-01-01'
end_date = '2024-12-01'

#get monthly imagery
desired_monthly_series = 'Global Quarterly'

monthly_mosaics = mosaics_in_series(session, desired_monthly_series, start_date=start_date, end_date=end_date)
for mosaic in monthly_mosaics:
    display(mosaic)
    print('\n')

#find monthly quads that intersect the aoi
def mosaic_metadata(session, mosaic_name):
    rv = session.get(f'{basemaps_api_url}/mosaics', params={'name__is': mosaic_name})
    rv.raise_for_status()
    return rv.json()['mosaics'][0]

def paginated_query(session, url, payload, item_key, **kwargs):
    rv = session.post(url, json=payload, **kwargs)
    rv.raise_for_status()
    page = rv.json()
    for item in page[item_key]:
        yield item
    
    if '_next' in page['_links']:
        url = page['_links']['_next']
        for item in paginated_get(session, url, item_key, **kwargs):
            yield item
            
def polygon_quad_search(session, mosaic_name, polygon):
    info = mosaic_metadata(session, mosaic_name)
    url = f'{basemaps_api_url}/mosaics/{info["id"]}/quads/search'
    return paginated_query(session, url, polygon, 'items')

#load basemap_client
with open(os.path.join(wd,"basemaps_client.py")) as f:
    exec(f.read())
client = BasemapsClient(api_key=None) 

# quads that intersect with desired series
series = client.series(name=desired_monthly_series)
for mosaic in series.mosaics(start_date=start_date, end_date=end_date):
    coverage = [quad.coverage for quad in mosaic.quads(region=aoi_geom)]
    avg = sum(coverage) / (len(coverage) or 1)
    print(f'{mosaic.name} has {len(coverage)} quads in the AOI averaging {avg:0.1f}% coverage')


#summarize quds


#%% download quads

#convert aoi_geom to bbox
coords = aoi_geom["coordinates"][0]

min_lon = min(coord[0] for coord in coords)  # Minimum longitude
max_lon = max(coord[0] for coord in coords)  # Maximum longitude
min_lat = min(coord[1] for coord in coords)  # Minimum latitude
max_lat = max(coord[1] for coord in coords)  # Maximum latitude

aoi_bbox = (min_lon, min_lat, max_lon, max_lat)

#create download directory and set wd to it
download_dir = os.path.join(wd, "data", "Quesnel_thinning", "planet_basemaps", 'monthly_normalized_8b_sr')

if not os.path.exists(download_dir):
    os.makedirs(download_dir)

#download quads
for item in series.download_quads(bbox=aoi_bbox, start_date=start_date, end_date=end_date):
    print(item)

# import glob
# import shutil

# for dirname in glob.glob('global_monthly_*'):
#     shutil.rmtree(dirname)
