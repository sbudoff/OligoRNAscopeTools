library(tidyverse)

ACD_ID_Assigner <- function(df) {
  data <- df %>%
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
  return(data)
}

ACDTranslator <- function(input_file_path, path, extra="",
                          GFP_thresh = 0,
                          # "Based on the average gene expression in scRNAseq dataset5, we used a cutoff of 12, 4, 3 ... of Ptgds+, Klk6+, Egr2+"
                          Egr2_thresh = 3,
                          Klk6_thresh = 4,
                          Ptgds_thresh = 0,
                          ReturnSummary=TRUE,
                          SaveResults=TRUE,
                          ComplexFileName=FALSE) {

  # Read the CSV output by the ACD software
  data <- read_csv(input_file_path) 
  
  # Assign IDs matching parents and children
  data <- ACD_ID_Assigner(data)
  
  print('IDs made')
  
  # Convert columns into more human readable names
  data <- data %>%
    filter(!is.na(ID)) %>%
    mutate(GFP = `Subcellular: Channel 2: Num spots estimated`,
           Cell_Area = `Cell: Area`,
           X = `Centroid X µm`,
           Y = `Centroid Y µm`,
           n_Klk6 = `Subcellular: Channel 3: Num single spots`,
           Klk6_puncta = `Subcellular: Channel 3: Num spots estimated`,
           intesity_Klk6 = `Subcellular cluster: Channel 3: Mean channel intensity`,
           area_Klk6 = `Subcellular cluster: Channel 3: Area`,
           n_Egr2 = `Subcellular: Channel 4: Num single spots`,
           Egr2_puncta = `Subcellular: Channel 4: Num spots estimated`,
           intesity_Egr2 = `Subcellular cluster: Channel 4: Mean channel intensity`,
           area_Egr2 = `Subcellular cluster: Channel 4: Area`,
           n_Ptgds = `Subcellular: Channel 5: Num single spots`,
           Ptgds_puncta = `Subcellular: Channel 5: Num spots estimated`,
           intesity_Ptgds = `Subcellular cluster: Channel 5: Mean channel intensity`,
           area_Ptgds = `Subcellular cluster: Channel 5: Area`,
           Location = Parent
    ) %>%
    select(ID, Image, Location, GFP, Cell_Area, X, Y, n_Egr2, Egr2_puncta, intesity_Egr2, area_Egr2,
           n_Klk6, Klk6_puncta, intesity_Klk6, area_Klk6, n_Ptgds, Ptgds_puncta, intesity_Ptgds, area_Ptgds) %>%
    mutate(X = ifelse(is.na(Cell_Area), NA, X),
           Y = ifelse(is.na(Cell_Area), NA, Y),
           Location = ifelse(Location=="PathCellObject", NA, Location),
           GFP = ifelse(is.na(Cell_Area), NA, GFP)) 
  print('Columns renamed for useability')
  
  # Extract the image locations for later
  Image_IDs <- data %>%
    select(ID, Image) %>%
    drop_na()
  
  # remove character columns
  data <- data %>%
    mutate(Location = toupper(str_extract(Location, "^[^_]+")),
           GM = ifelse(Location == "GM", 1, 0),
           CC = ifelse(Location == "CC", 1, 0) + ifelse(Location == "WM", 1, 0) + ifelse(Location == "DORSAL_WM", 1, 0),
           L13 = ifelse(Location == "L13", 1, 0) + ifelse(Location == "L23", 1, 0),
           L4 = ifelse(Location == "L4", 1, 0),
           L56 = ifelse(Location == "L56", 1, 0)) %>%
    select(-c(Image, Location))
  
  print('Layers Identified')
  
  # Remove NAs
  data[is.na(data)] <- 0
  
  # Compute normalized metrics for possible classification schemes
  flattened_data <- data %>%
    group_by(ID) %>%
    mutate(Egr2 = intesity_Egr2*area_Egr2,
           Klk6 = intesity_Klk6*area_Klk6, 
           Ptgds = intesity_Ptgds*area_Ptgds) %>%
    select(-c(intesity_Egr2,area_Egr2,intesity_Klk6,area_Klk6,intesity_Ptgds,area_Ptgds)) %>%
    summarise_all(sum) %>%
    mutate(MOL1 = Egr2_puncta > Egr2_thresh,
           MOL2 = Klk6_puncta > Klk6_thresh,
           `MOL5/6` = Ptgds_puncta > Ptgds_thresh,
           Egr2=Egr2/max(Egr2),
           Klk6=Klk6/max(Klk6),
           Ptgds=Ptgds/max(Ptgds),
           Egr2_pp = Egr2/Egr2_puncta,
           Klk6_pp = Klk6/Klk6_puncta,
           Ptgds_pp = Ptgds/Ptgds_puncta)
  
  print("Integreited intesnity computed")
  
  flattened_data[is.na(flattened_data)] <- 0
  
  # Put image IDs back in
  flattened_data <- left_join(Image_IDs, flattened_data) %>%
    unique()
  
  # Extract animal ID
  if (ComplexFileName){
    flattened_data <- flattened_data %>%
      mutate(filename = str_extract(Image, "_(.*?)\\."),
             filename = str_remove(filename, ".*_"),
             Animal = str_extract(filename, "\\d+"),
             Animal = factor(Animal))
  } else {
    flattened_data <- flattened_data %>%
      mutate(Animal = str_extract(Image, "\\d+"),
             Animal = factor(Animal))
  }
  
  print("Animal name extracted from file name")
  
  # Simplest classification scheme, which normalized value is greatest in a given cell
  Oligo_data <- flattened_data %>%
    group_by(Image) %>%
    filter(GFP > GFP_thresh) %>%
    ungroup() %>%
    mutate(Classification = case_when(
      Egr2 > Klk6 & Egr2 > Ptgds ~ "Egr2",
      Klk6 > Egr2 & Klk6 > Ptgds ~ "Klk6",
      Ptgds > Egr2 & Ptgds > Klk6 ~ "Ptgds",
      TRUE ~ "None"
      ))
  
  print("Simple classification scheme applied")
  
  # Simple summary table counting Classified Oligos in a given region
  Oligo_perRegion_perAnimal_Classified <- Oligo_data %>%
    group_by(Classification, Animal, GM, CC, L13, L4, L56) %>%
    summarise(N = n() ) %>%
    ungroup() %>%
    group_by(Animal, GM, CC, L13, L4, L56) %>%
    mutate(Proportion = N / sum(N))
  
  Oligo_perRegion_Classified <- Oligo_perRegion_perAnimal_Classified %>%
    group_by(Classification, GM, CC, L13, L4, L56) %>%
    summarise(Mean = mean(Proportion) ,
              SD = sd(Proportion)) %>%
    ungroup() 
  
  print("Classification summaries produced")
  
  # MOL summary table counting Oligos based on published identity
  MOL_perRegion_perAnimal <- Oligo_data %>%
    group_by(MOL1, MOL2, `MOL5/6`, Animal, GM, CC, L13, L4, L56) %>%
    summarise(N = n() ) %>%
    ungroup() %>%
    group_by(Animal, GM, CC, L13, L4, L56) %>%
    mutate(Proportion = N / sum(N))
  
  
  MOL_perRegion <- MOL_perRegion_perAnimal %>%
    group_by(MOL1, MOL2, `MOL5/6`, GM, CC, L13, L4, L56) %>%
    summarise(Mean = mean(Proportion) ,
              SD = sd(Proportion)) %>%
    ungroup() 
  
  print('MOL based summaries produced')
  
  if (SaveResults) {
    # Save outputs
    write.csv(flattened_data, file = paste0(path,extra,"measurements_all.csv"))
    write.csv(Oligo_data, file = paste0(path,extra,"measurements_Oligo.csv"))
    write.csv(Oligo_perRegion_perAnimal_Classified, file = paste0(path,extra,"summary_Oligo_by_Animal.csv"))
    write.csv(Oligo_perRegion_Classified, file = paste0(path,extra,"summary_Oligo.csv"))
    write.csv(MOL_perRegion_perAnimal, file = paste0(path,extra,"summary_MOL_by_Animal.csv"))
    write.csv(MOL_perRegion, file = paste0(path,extra,"summary_MOL.csv"))
    print('Data saved')
  }
  
  if (ReturnSummary) {
    return(Oligo_perRegion_perAnimal)
  } else {
    out <- list()
    out[['Summary']] <- MOL_perRegion_perAnimal
    out[['Summary Simple']] <- MOL_perRegion
    out[['Summary Classified']] <- Oligo_perRegion_perAnimal_Classified
    out[['Summary Simple Classified']] <- Oligo_perRegion_Classified
    out[['All']] <- flattened_data
    out[['Oligos']] <- Oligo_data
    return(out)
  }

}

Oligo_perRegion_list_SC <- ACDTranslator(input_file_path = '/media/sam/STOCKTON/measurements_redo08072023.csv', 
                                            path = "/home/sam/Downloads/",
                                            extra = "SC",
                                            ReturnSummary = F)  
Oligo_perRegion_perAnimal_SC <- Oligo_perRegion_list_SC[['Summary']] 

Oligo_perRegion_perAnimal_SC <- Oligo_perRegion_perAnimal_SC %>%
  mutate(Anatomy="SC")
  
Oligo_perRegion_list <- ACDTranslator(input_file_path = '/home/sam/cellprofiler_sam/CellProfiler/measurements_073123.csv',
                                          path = "/home/sam/Downloads/",
                                          extra = "Cortex",
                                          ReturnSummary = F)
Oligo_perRegion_perAnimal_Cortex <- Oligo_perRegion_list[['Summary']] 

Oligo_perRegion_perAnimal_Cortex <- Oligo_perRegion_perAnimal_Cortex%>%
  mutate(Anatomy="Cortex")

cup_animals <- c(6128,6142)

summary_df_by_animal <- rbind(Oligo_perRegion_perAnimal_Cortex, Oligo_perRegion_perAnimal_SC) %>%
  gather(Layer, Val, -c(MOL1, MOL2, `MOL5/6`, Animal, Anatomy, N, Proportion))%>%
  filter(Val != 0 ) %>%
  select(-Val) %>%
  mutate(Healthy = case_when(Animal %in% cup_animals ~ 0,
                             TRUE ~ 1),
         Animal = factor(Animal),
         Anatomy = factor(Anatomy),
         Layer = factor(Layer))

summary_df <- summary_df_by_animal %>%
  group_by(MOL1, MOL2, `MOL5/6`, Layer, Anatomy, Healthy) %>%
  summarise(Mean = mean(Proportion) ,
            SD = sd(Proportion)) %>%
  ungroup() 

# Save outputs
write.csv(summary_df, file = paste0("/home/sam/Downloads/","Oligo_finalSummary.csv"))


summary_df_by_animal_cort <- summary_df_by_animal %>%
  filter(Anatomy=="Cortex",
         Healthy != 1)

library(patchwork)

# Define the desired order of the categories
layer_order <- c("L13", "L4", "L56", "CC", "GM")

# Convert the Category column to a factor with the desired order
summary_df_by_animal_cort$Layer <- factor(summary_df_by_animal_cort$Layer, levels = rev(layer_order))

MOL1_bp <- summary_df_by_animal_cort %>%
  filter(MOL1==1) %>%
  ggplot(aes(x = Proportion, y = Layer)) +
  geom_boxplot() +
  ggtitle('MOL1') +
  theme_minimal()
MOL2_bp <- summary_df_by_animal_cort %>%
  filter(MOL2==1) %>%
  ggplot(aes(x = Proportion, y = Layer)) +
  geom_boxplot() +
  ggtitle('MOL2') +
  theme_minimal()
MOL56_bp <- summary_df_by_animal_cort %>%
  filter(`MOL5/6`==1) %>%
  ggplot(aes(x = Proportion, y = Layer)) +
  geom_boxplot() +
  ggtitle('MOL5/6') +
  theme_minimal()

MOL1_bp+MOL2_bp+MOL56_bp


# ANOVA
summary_df_by_animal <- summary_df_by_animal %>%
  filter(Layer != 'GM')
MOL1_twoWay <- aov(Proportion ~ MOL1 + Layer + MOL1 * Layer, data = summary_df_by_animal_cort)
MOL2_twoWay <- aov(Proportion ~ MOL2 + Layer + MOL2 * Layer, data = summary_df_by_animal_cort)
MOL56_twoWay <- aov(Proportion ~ `MOL5/6` + Layer + `MOL5/6` * Layer, data = summary_df_by_animal_cort)
summary(MOL1_twoWay)
summary(MOL2_twoWay)
summary(MOL56_twoWay)

MOL_plot_df <- Oligo_perRegion_list[['Oligos']] %>%
  select(X, Y, MOL1, MOL2, `MOL5/6`, Animal, CC, L13, L4, L56, Egr2_puncta, Klk6_puncta, Ptgds_puncta) %>%
  gather(Layer, Val, -c(X, Y, MOL1, MOL2, `MOL5/6`, Animal, Egr2_puncta, Klk6_puncta, Ptgds_puncta))%>%
  filter(Val != 0 ) %>%
  select(-Val) %>%
  mutate(Healthy = case_when(Animal %in% cup_animals ~ 0,
                             TRUE ~ 1),
         Animal = factor(Animal),
         Layer = factor(Layer)) 

MOL_plot_df %>%
  filter(MOL1) %>%
  ggplot(aes(X, Y, shape=Layer, color = Egr2_puncta)) +
  geom_point(alpha=0.5) +
  scale_colour_gradient(name = "Egr2", low = "blue", high = "red") +
  theme_minimal() +
  scale_y_reverse() +
  facet_wrap(~Animal)

MOL_plot_df %>%
  filter(MOL2) %>%
  ggplot(aes(X, Y, shape=Layer, color = Egr2_puncta)) +
  geom_point(alpha=0.5) +
  scale_colour_gradient(name = "Egr2", low = "blue", high = "red") +
  theme_minimal() +
  scale_y_reverse() +
  facet_wrap(~Animal)

MOL_plot_df %>%
  filter(`MOL5/6`) %>%
  ggplot(aes(X, Y, shape=Layer, color = Ptgds_puncta)) +
  geom_point(alpha=0.5) +
  scale_colour_gradient(name = "Egr2", low = "blue", high = "red") +
  theme_minimal() +
  scale_y_reverse() +
  facet_wrap(~Animal)


# summary_df_by_animal_cort %>%
#   ggplot(aes(0, y = Proportion, color = Classification, fill = Classification)) +
#   geom_violin() +
#   ggtitle("Distribution of Oligos Parsed By RNAscope Classification") +
#   theme_minimal() +
#   theme(axis.text.x = element_blank(),
#         axis.title.x = element_blank())+
#   facet_grid(~Layer*Healthy)

# library(plotly)
# 
# plot_ly(x=data$Egr2, y=data$Klk6, z=data$Ptgds, type="scatter3d", mode="markers", color=data$Classification)