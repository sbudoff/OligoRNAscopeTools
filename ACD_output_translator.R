library(dplyr)

data <- readr::read_csv("/home/sam/scRNAseq/RNAscope/CellProfiler/git/measurements_all_v2.csv") 

data <- data %>%
  mutate(ID = NA)

# Iterate over each row in the data
for (i in 2:nrow(data)) {
  # Check if the current row's "Class" is not NA
  if (!is.na(data$Class[i])) {
    # Check if the prior row's "Class" is NA
    if (is.na(data$Class[i - 1])) {
      # Assign the row number of the row with NA to ParentID
      data$ID[i] <- i - 1
      data$ID[i-1] <- i - 1
    } else {
      # Check if the Class of the prior row is also not NA
      # Assign the same value in ParentID as the prior row
      data$ID[i] <- data$ID[i - 1]
    }
  } 
}

data <- data %>%
  filter(!is.na(ID)) %>%
  mutate(GFP = `Cell: Channel 2 max`,
         Cell_Area = `Cell: Area`,
         X = `Centroid X µm`,
         Y = `Centroid Y µm`,
         n_Egr2 = `Subcellular: Channel 3: Num single spots`,
         intesity_Egr2 = `Subcellular cluster: Channel 3: Mean channel intensity`,
         area_Egr2 = `Subcellular cluster: Channel 3: Area`,
         n_Klk6 = `Subcellular: Channel 4: Num single spots`,
         intesity_Klk6 = `Subcellular cluster: Channel 4: Mean channel intensity`,
         area_Klk6 = `Subcellular cluster: Channel 4: Area`,
         n_Ptgds = `Subcellular: Channel 5: Num single spots`,
         intesity_Ptgds = `Subcellular cluster: Channel 5: Mean channel intensity`,
         area_Ptgds = `Subcellular cluster: Channel 5: Area`,
         ) %>%
  select(ID, Image, GFP, Cell_Area, X, Y, n_Egr2, intesity_Egr2, area_Egr2,
         n_Klk6, intesity_Klk6, area_Klk6, n_Ptgds, intesity_Ptgds, area_Ptgds) %>%
  mutate(X = ifelse(is.na(Cell_Area), NA, X),
         Y = ifelse(is.na(Cell_Area), NA, Y),
         GFP = ifelse(is.na(Cell_Area), NA, GFP)) 

flattened_data <- data %>%
  group_by(ID) %>%
  summarise_all(~ if(all(is.na(.))) NA else na.omit(.)[1])


GFP_thresh = 125

OPC_data <- flattened_data %>%
  filter(GFP > GFP_thresh)

write.csv(flattened_data, file = "/home/sam/scRNAseq/RNAscope/CellProfiler/git/measurements_all_useable.csv")
write.csv(OPC_data, file = "/home/sam/scRNAseq/RNAscope/CellProfiler/git/measurements_OPC_useable.csv")
