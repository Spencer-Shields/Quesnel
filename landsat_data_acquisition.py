import earthaccess
import os
import geopandas as gpd
from pprint import pprint

#set working directory to project folder
wd_path = r"D:\Quesnel"
os.chdir(wd_path)

#set data search parameters

start_date = '2000-01-01T00:00:00.000Z'
end_date = '2025-01-01T00:00:00.000Z'
aoi_path = 'data\Quesnel_thinning\AOI_fullsite_wgs84.geojson'
aoi = gpd.read_file(aoi_path)
aoi_bbox = (aoi.total_bounds[0],
            aoi.total_bounds[1],
            aoi.total_bounds[2],
            aoi.total_bounds[3])

#authenticate to access Earthdata API
os.environ["EARTHDATA_USERNAME"] = "spenshi"
os.environ["EARTHDATA_PASSWORD"] = "3@rth D@t@ Plz"

earthaccess.login(strategy='environment', persist=True)

#search for HLSL30 results
results = earthaccess.search_data(
    #short_name = "HLSL30",
    short_name = ['HLSL30'],
    bounding_box = aoi_bbox,
    temporal = (start_date, end_date),
    cloud_cover = (0,75)
)

print(f"Found {len(results)} results")

#define output directory, download files
save_hlsl30_dir = r'data\HLS\raw_hlsl30'
os.makedirs(save_hlsl30_dir, exist_ok=True)

if os.path.exists(save_hlsl30_dir) == True:
    files = earthaccess.download(results, save_hlsl30_dir)
    print(f"Downloaded {len(files)} files")

#search for HLSS30 results
results = earthaccess.search_data(
    #short_name = "HLSL30",
    short_name = ['HLSS30' ],
    bounding_box = aoi_bbox,
    temporal = (start_date, end_date),
    cloud_cover = (0,75)
)

print(f"Found {len(results)} results")

#define output directory, download files
save_hlss30_dir = r'data\HLS\raw_hlsl30'
os.makedirs(save_hlss30_dir, exist_ok=True)

if os.path.exists(save_hlss30_dir) == True:
    files = earthaccess.download(results, save_hlss30_dir)
    print(f"Downloaded {len(files)} files")


#search for landsat results

datasets = earthaccess.search_datasets(
            keyword="landsat 7 level 2",
            temporal = (
                '2000-01-01',
                '2012-12-31'
            ),
            bounding_box = aoi_bbox,
)
len(datasets)

for i, dataset in enumerate(datasets[:26]):  # Show first 5
    print(f"\n--- Dataset {i+1} ---")
    print(f"Short Name: {dataset.get('umm', {}).get('ShortName', 'N/A')}")
    print(f"Provider: {dataset.get('meta', {}).get('provider-id', 'N/A')}")
    print(f"Title: {dataset.get('umm', {}).get('EntryTitle', 'N/A')}")

results = earthaccess.search_data(
    short_name = [
        'LEO7'
    ],
    #concept_id = "C3442503899-USGS_EROS",
    bounding_box = aoi_bbox,
    temporal = ('2000-01-01T00:00:00.0000Z', '2012-12-31T23:59:59.999Z'),
    downloadable = True    # <--- Ensure it's available via standard archive
    #cloud_cover = (0,75)
)

print(f"Found {len(results)} results")


import json
from landsatxplore.api import API






#define output directory, download files
save_landsat_dir = r'data\HLS\raw_landsat'
os.makedirs(save_landsat_dir, exist_ok=True)

if os.path.exists(save_landsat_dir) == True:
    files = earthaccess.download(results, save_landsat_dir)
    print(f"Downloaded {len(files)} files")
