library(vegan)
library(raster)
library(vegan)
library(raster)
library(factoextra)
library(ggplot2)

# Loading environmental raster data
env_rasters <- stack(list.files(path = "C:\\Users\\doang\\OneDrive\\Desktop\\Gen2 Modeling\\PIma\\Tif Rasters", pattern = "\\.tif$", full.names = TRUE))

# Loading species presence points and removing the 'Species' column
species_points <- read.csv("C:\\Users\\doang\\OneDrive\\Desktop\\Gen2 Modeling\\PIma\\Presance\\PreancesMain.csv")
species_points$Species <- NULL # Assuming you've correctly removed the column

# Setting coordinates for species_points; ensure Lat and Lon columns exist
coordinates(species_points) <- ~Lon+Lat # Corrected order based on typical naming convention

# Extracting environmental data at species presence locations
env_data <- extract(env_rasters, species_points)

# Converting the extracted data to a dataframe
env_data_df <- as.data.frame(env_data)

new_variable_names <- paste("Bio", 1:19, sep = "")

# Assign the new variable names to the dataframe columns
names(env_data_df_clean) <- new_variable_names
# Identifying columns that are entirely NA
na_columns <- apply(env_data_df, 2, function(x) all(is.na(x)))

# Removing columns that are entirely NA from env_data_df
env_data_df_clean <- na.omit(env_data_df)
env_data_df_clean <- as.data.frame(env_data_df_clean)
# Run PCA
pca_results <- prcomp(env_data_df_clean, scale. = TRUE)

# Plot
fviz_pca_var(pca_results, 
             col.var = "red",      # Color for variables
             gradient.cols = "red",
             title = "Pima ENV Variable PCA",# If you want a gradient of colors
             arrows = TRUE)  





##SA
# Loading environmental raster data
SA_env_rasters <- stack(list.files(path = "C:\\Users\\doang\\OneDrive\\Desktop\\Gen2 Modeling\\SA\\Cliped Rasters", pattern = "\\.tif$", full.names = TRUE))

# Loading species presence points and removing the 'Species' column
SA_species_points <- read.csv("C:\\Users\\doang\\OneDrive\\Desktop\\Gen2 Modeling\\SA\\Presance\\SAPrez.csv")
SA_species_points$species <- NULL # Assuming you've correctly removed the column

# Setting coordinates for species_points; ensure Lat and Lon columns exist
coordinates(SA_species_points) <- ~Lat+Lon

# Extracting environmental data at species presence locations
SA_env_data <- extract(SA_env_rasters, SA_species_points)

# Converting the extracted data to a dataframe
SA_env_data_df <- as.data.frame(SA_env_data)

new_variable_names <- paste("Bio", 1:19, sep = "")

# Assign the new variable names to the dataframe columns
names(SA_env_data_df) <- new_variable_names
# Identifying columns that are entirely NA
SA_na_columns <- apply(SA_env_data_df, 2, function(x) all(is.na(x)))

# Removing columns that are entirely NA from env_data_df
SA_env_data_df_clean <- na.omit(SA_env_data_df)
env_data_df_clean <- as.data.frame(env_data_df_clean)
# Run PCA
SA_pca_results <- prcomp(SA_env_data_df_clean, scale. = TRUE)

# Plot
fviz_pca_var(SA_pca_results, 
             col.var = "orange",      # Color for variables
             gradient.cols = "orange",
             title = "South Africa ENV Variable PCA",
             arrows = TRUE)  
