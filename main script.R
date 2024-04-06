#LIBRARYS ----------------------------------------------------------------------
library(vegan)
library(raster)
library(factoextra)
library(ggplot2)
library(sp)
library(sf)
library(scatterplot3d)
library(car)
library(rgl)
library(ggrepel)
# PIMA Analysis -----------------------------------------------------------

# Load Pima environmental raster data
pima_env_rasters <- stack(list.files(path = "C:/Users/doang/OneDrive/Desktop/Gen2 Modeling/PIma/Tif Rasters", 
                                     pattern = "\\.tif$", full.names = TRUE))

# Load Pima species presence points and remove the 'Species' column
pima_species_points <- read.csv("C:/Users/doang/OneDrive/Desktop/Gen2 Modeling/PIma/Presance/PreancesMain.csv")
pima_species_points$Species <- NULL
# Set coordinates for Pima species_points; ensure Lat and Lon columns exist
coordinates(pima_species_points) <- ~Lon+Lat

# Extract Pima environmental data at species presence locations
pima_env_data <- extract(pima_env_rasters, pima_species_points)

# Convert Pima extracted data to a dataframe and remove NAs
pima_env_data_df <- as.data.frame(pima_env_data)
pima_env_data_df <- na.omit(pima_env_data_df)

# Assign new variable names for Pima
pima_new_variable_names <- paste("Bio", 1:ncol(pima_env_data_df), sep = "")
names(pima_env_data_df) <- pima_new_variable_names

# Run PCA for Pima
pima_pca_results <- prcomp(pima_env_data_df, scale. = TRUE)

# Plot PCA for Pima
fviz_pca_var(pima_pca_results, col.var = "red", gradient.cols = c("white", "red"), title = "Pima ENV Variable PCA")

# SOUTH AFRICA Analysis ---------------------------------------------------

# Load South Africa environmental raster data
sa_env_rasters <- stack(list.files(path = "C:/Users/doang/OneDrive/Desktop/Gen2 Modeling/SA/Cliped Rasters", 
                                   pattern = "\\.tif$", full.names = TRUE))

# Load South Africa species presence points and remove the 'Species' column
sa_species_points_df <- read.csv("C:/Users/doang/OneDrive/Desktop/Gen2 Modeling/SA/Presance/SAPrez.csv")
sa_species_points_sf <- st_as_sf(sa_species_points_df, coords = c("Lon", "Lat"), crs = 4326)
sa_species_points_sp <- as(sa_species_points_sf, "Spatial")
sa_species_points <- as.data.frame(sa_species_points_sp)
sa_species_points$species <- NULL
# Set coordinates for South Africa species points; ensure Lat and Lon columns exist
coordinates(sa_species_points) <- ~coords.x1+coords.x2

# Extract South Africa environmental data at species presence locations
sa_env_data <- extract(sa_env_rasters, sa_species_points)

# Convert South Africa extracted data to a dataframe and remove NAs
sa_env_data_df <- as.data.frame(sa_env_data)
sa_env_data_df <- na.omit(sa_env_data_df)

# Assign new variable names for South Africa
sa_new_variable_names <- paste("Bio", 1:ncol(sa_env_data_df), sep = "")
names(sa_env_data_df) <- sa_new_variable_names

# Remove columns with zero variance
sa_env_data_df_clean <- sa_env_data_df[, sapply(sa_env_data_df, function(x) var(x, na.rm = TRUE) != 0)]

# Check for any remaining NA values, which should not exist after na.omit() but just in case
sa_env_data_df_clean <- na.omit(sa_env_data_df_clean)

# Run PCA on the clean and validated dataset
sa_pca_results <- prcomp(sa_env_data_df_clean, scale. = TRUE)

# Plot PCA for South Africa
fviz_pca_var(sa_pca_results, col.var = "orange", gradient.cols = c("white", "orange"), title = "South Africa ENV Variable PCA")

# Variable Impartaince  ---------------------------------------------------------------------------------------------
pima_loadings <- pima_pca_results$rotation[, 1:2] # Using first two principal components for example
sa_loadings <- sa_pca_results$rotation[, 1:2]

# Rank variables based on the absolute values of loadings
pima_importance <- sort(abs(pima_loadings[, 1]), decreasing = TRUE)
sa_importance <- sort(abs(sa_loadings[, 1]), decreasing = TRUE)

# Create a combined data frame for the plot
importance_df <- data.frame(
  Variable = rep(names(pima_importance), 2),
  Importance = c(pima_importance, sa_importance),
  Region = rep(c("Pima", "South Africa"), each = length(pima_importance))
)

# Plot using ggplot2
library(ggplot2)

ggplot(importance_df, aes(x = Variable, y = Importance, fill = Region)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Importance of Environmental Variables in Pima and South Africa",
       x = "Environmental Variables",
       y = "Importance (Absolute PCA Loadings)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#DISPLAYS SPECIES OCRANCE OVER PCA-ENV -------------------------------------------------------------------------------

pima_scores <- predict(pima_pca_results, newdata = pima_env_data_df)
sa_scores <- predict(sa_pca_results, newdata = sa_env_data_df_clean)

combined_pca_scores <- rbind(
  data.frame(PC1 = pima_scores[, 1], PC2 = pima_scores[, 2], Region = 'Pima'),
  data.frame(PC1 = sa_scores[, 1], PC2 = sa_scores[, 2], Region = 'South Africa')
)

ggplot() +
  geom_point(data = combined_pca_scores, aes(x = PC1, y = PC2, color = Region)) +
  geom_hline(yintercept = 0, linetype = "solid") +
  geom_vline(xintercept = 0, linetype = "solid") +
  labs(title = "PCA Biplot with Species Occurrence Data", 
       x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red"))

#Niche Overlap for full range  --------------------------------------------------------
pima_sd <- pima_pca_results$sdev
sa_sd <- sa_pca_results$sdev

# Proportion of variance explained by each principal component
pima_var_explained <- round(pima_sd^2 / sum(pima_sd^2) * 100, 1)
sa_var_explained <- round(sa_sd^2 / sum(sa_sd^2) * 100, 1)

top_vars_pima <- c("Bio1", "Bio2", "Bio8", "Bio5", "Bio10")
top_vars_sa <- c("Bio18", "Bio1", "Bio2", "Bio8", "Bio5")

# Filter the combined_loadings for top variables
top_loadings <- combined_loadings[combined_loadings$Variable %in% top_vars_pima | 
                                    combined_loadings$Variable %in% top_vars_sa, ]

# Calculate a scale factor for the arrows
max_score_extent <- max(abs(combined_pca_scores$PC1), abs(combined_pca_scores$PC2))
loading_scale_factor <- max_score_extent / 5

# Apply the scale factor to the filtered top loadings
top_loadings$PC1 <- top_loadings$PC1 * loading_scale_factor
top_loadings$PC2 <- top_loadings$PC2 * loading_scale_factor

# Create a PCA biplot
pca_biplot <- ggplot(data = combined_pca_scores, aes(x = PC1, y = PC2, color = Region)) +
  geom_point(alpha = 0.5) +
  geom_segment(data = top_loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2, color = Region),
              arrow = arrow(length = unit(0.2, "inches")), size = 1.5, lineend = 'butt') +
  geom_text_repel(data = top_loadings, aes(label = Variable), size = 3,
                  box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines")) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid", size = 0.5) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
  stat_ellipse(data = combined_pca_scores[combined_pca_scores$Region == "Pima", ],
               aes(fill = Region), geom = "polygon", alpha = 0.2, level = 0.95) +
  stat_ellipse(data = combined_pca_scores[combined_pca_scores$Region == "South Africa", ],
               aes(fill = Region), geom = "polygon", alpha = 0.2, level = 0.95) +
  scale_color_manual(values = c("Pima" = "#1b9e77", "South Africa" = "#d95f02"),
   labels = c("Invaded Range (Pima)", "Native Range (South Africa)")) +
  scale_fill_manual(values = c("Pima" = "#66c2a5", "South Africa" = "#fc8d62"))+
  labs(title = "Comparison of Environmental Niches: Stinknet (Oncosiphon pilulifer) in Invaded vs. Native Ranges",
       x = paste("PC1 (", pima_var_explained[1], "% variance explained)"), 
       y = paste("PC2 (", pima_var_explained[2], "% variance explained)")) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  guides(fill = FALSE)
 



show(pca_biplot)
