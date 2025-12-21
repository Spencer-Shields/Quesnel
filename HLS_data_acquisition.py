import earthaccess
import os
import geopandas as gpd

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

#search for HLS results
results = earthaccess.search_data(
    short_name = "HLSL30",
    bounding_box = aoi_bbox,
    temporal = (start_date, end_date),
    cloud_cover = (0,75)
)

print(f"Found {len(results)} results")

#define output directory, download files
save_dir = r'data\HLS_2000-2025\raw'
os.makedirs(save_dir, exist_ok=True)

if os.path.exists(save_dir) == True:
    files = earthaccess.download(results, save_dir)
    print(f"Downloaded {len(files)} files")
