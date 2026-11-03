# Install packages if needed
# install.packages("rgbif")
# install.packages("dplyr")

library(rgbif)
library(dplyr)

# Search for species name and get GBIF taxon key
name <- "Acanthochromis polyacanthus"

species_key <- name_backbone(name = name)$speciesKey
print(paste("GBIF taxon key:", species_key))

# Download occurrence data (limit to 2000 records for quick testing)
occ_data <- occ_search(taxonKey = species_key,
                       limit = 2000,
                       hasCoordinate = TRUE)

# Extract data frame of occurrences
occ_df <- occ_data$data

# Filter and select only longitude & latitude
coords <- occ_df %>%
  select(species, decimalLongitude, decimalLatitude, country, year) %>%
  rename(longitude = decimalLongitude,
         latitude = decimalLatitude) %>%
  filter(!is.na(longitude), !is.na(latitude))

# Preview
head(coords)
getwd()
setwd("D:/rSDM_marine_fish/")
dir.create("data")
list.files()

# Save to CSV (optional)
write.csv(coords, "D:/rSDM_marine_fish/data/Acanthochromis_polyacanthus.csv", row.names = FALSE)

#view contents if file now exists
list.files()

