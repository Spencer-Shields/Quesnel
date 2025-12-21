library(tools)

#The purpose of this script is to check if processed scene files have been created which don't
#correspond with the original raw data

scenes_dir = 'data/planet_scenes/additional_scenes' #directory of processed data
raw_dir = 'data/planet_scenes/raw_additional' #directory for raw data

#define new function to extract metadata
ps_meta <- function(dir, recursive = TRUE) {
  # Returns a tibble of PlanetScope metadata based on metadata.json files in a directory.
  
  # Use pattern argument instead of filtering afterwards
  jsons <- list.files(dir, pattern = "metadata\\.json$", 
                      recursive = recursive, full.names = TRUE)
  
  if (length(jsons) == 0) {
    message('No metadata files found in directory.')
    return(NULL)
  }
  
  # Pre-allocate list for better performance
  jsons_l <- vector("list", length(jsons))
  
  # Read all JSON files
  for (i in seq_along(jsons)) {
    jsons_l[[i]] <- jsonify::from_json(jsons[i])
  }
  
  # Extract id and properties more efficiently
  ids <- vapply(jsons_l, `[[`, character(1), "id")
  props <- lapply(jsons_l, `[[`, "properties")
  
  # Combine into tibble
  json_df <- tibble::tibble(id = ids)
  
  # Get all unique property names across all files
  all_props <- unique(unlist(lapply(props, names)))
  
  # Extract each property
  for (prop in all_props) {
    json_df[[prop]] <- lapply(props, function(x) x[[prop]] %||% NA)
  }
  
  # Unnest single-value columns
  json_df <- tibble::as_tibble(lapply(json_df, function(col) {
    if (all(vapply(col, function(x) length(x) <= 1, logical(1)))) {
      unlist(col)
    } else {
      col
    }
  }))
  
  return(json_df)
}

#make metadata dataframe
meta_df = ps_meta(raw_dir)

#get ids of scenes
ids = meta_df$id

#delete extra files
if(length(ids ==0)){
  message('No ids')
} else {

fl = list.files(scenes_dir, recursive = T, full.names = T, pattern='\\.tif$')

deleted_list = c()

for(i in 1:length(fl)){
  #get file
  f = fl[i]
  
  cat('\rChecking',i,'/',length(fl))
  
  #get file id
  f_id = basename(file_path_sans_ext(f))
  
  #remove _NULL from file id
  f_id = str_replace(f_id, '_NULL', '')
  
  #check if file id exists in ids extracted from metadata, delete if no
  if(file.exists(f)&!f_id %in% ids){
    file.remove(f)
    #add old filepath to deleted list
    deleted_list = append(deleted_list,f)
    cat('\r---DELETED')
    }
}

cat('\n',length(deleted_list), 'items deleted')
}