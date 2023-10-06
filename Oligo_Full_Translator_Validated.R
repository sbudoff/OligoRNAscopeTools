################################################################################
########################  Hyperparameters  #####################################
################################################################################

ID <- "Cortex"
input_file_path <- '/home/sam/scRNAseq/RNAscope/CellProfiler/git/measurements_allmice_09082023.csv'
save_path <- "/home/sam/scRNAseq/RNAscope/CellProfiler/git/output/"
oligos_df_existing <- "SCmeasurements_Oligo.csv"
cup_4d_animals <- c(6127,6140,6142,6126,6141)
cup_7w_animals <- c(5988,5989,5990,5991,5992,5994,5993)
healthy_p60_animals <- c(6740,6741,6794,6795,6904,6905)
dropGM <- T
Egr2_thresh <- 14
Klk6_thresh <- 12
Ptgds_thresh <-7

################################################################################
#####################  Functions/Libraries  ####################################
################################################################################

if (!require(tidyverse) ) {
  install.packages(c("svglite"),dep=TRUE) 
}
if (!require(svglite) ) {
  install.packages(c("svglite"),dep=TRUE) 
}
if (!require(patchwork) ) {
  install.packages(c("patchwork"),dep=TRUE) 
}

library(tidyverse)
library(patchwork)
library(svglite)

#################################################################################

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
  
  # Save outputs
  oligo_path <- paste0(path,extra,"measurements_Oligo.csv")
  write.csv(flattened_data, file = paste0(path,extra,"measurements_all.csv"))
  write.csv(Oligo_data, file = oligo_path)
  print('Data saved')
  
  return(oligo_path)

}

#################################################################################

read_oligos <- function(path, cup_animals,buffer=0, dropGM=T) {
  out <- read_csv(file = path) %>%
    select(-`...1`) %>%
    mutate(Anatomy="Cortex") %>%
    mutate(Condition = case_when(Animal %in% cup_4d_animals ~ "Cup4d",
                                 Animal %in% cup_7w_animals ~ "Cup7w",
                                 Animal %in% healthy_p60_animals ~ "HealthyP60",
                                 TRUE ~ "HealthyP140"),
           Condition = factor(Condition),
           Animal = factor(Animal),
           Anatomy = factor(Anatomy),
           Classification = factor(Classification),
           ID = ID+buffer)  %>%
    mutate(WM = ifelse(GM==0,1,0))
  
  if (dropGM) {
    layer_df <- out %>%
      filter(L13+L4+L56+CC==1) %>%
      select(L13,L4,L56,CC, ID)
  } else {
    layer_df <- out %>%
      select(GM,WM, ID)
  }
  
  layer_df <- layer_df %>%
    gather("Layer", "Val", -ID)%>%
    filter(Val != 0 ) %>%
    select(-Val) %>%
    mutate(Layer=factor(Layer)) %>%
    distinct()
  
  out <- right_join(layer_df,out, by="ID") %>%
    drop_na()
  
  return(out)
}

#################################################################################
SummaryMaker <- function(oligos_df) {
  # Validated cell proportion computer and layer relabeler
  
  oligos_sum <- oligos_df %>%
    group_by(Animal, Layer) %>%
    mutate(TotalCells = n(),
           N_MOL1 = sum(MOL1),
           N_MOL2 = sum(MOL2),
           N_MOL56 = sum(`MOL5/6`),
           `MOL1 Proportion` = N_MOL1/TotalCells,
           `MOL2 Proportion` = N_MOL2/TotalCells,
           `MOL5/6 Proportion` = N_MOL56/TotalCells) %>%
    ungroup() %>%
    distinct() 
  
  # Define the desired order of the categories
  layer_order <- c("L13", "L4", "L56", "CC", "GM", "WM")
  
  # Convert the Category column to a factor with the desired order
  oligos_sum$Layer <- factor(oligos_sum$Layer, levels = rev(layer_order))
  
  return(oligos_sum)
}

#################################################################################

oligo_violin_plotter <- function(oligos_sum_simple) {
  MOL1_bp <- oligos_sum_simple %>%
    filter(MOL1==1) %>%
    group_by(Animal, Layer, Anatomy, Condition) %>%
    select(`MOL1 Proportion`)%>%
    mutate(Proportion = `MOL1 Proportion`) %>%
    distinct() %>%
    ungroup() %>%
    ggplot(aes(x = Proportion, y = Layer)) +
    geom_violin() +
    geom_jitter(aes(color=Animal), alpha=0.6) +
    ggtitle('MOL1') +
    xlim(0,1) +
    theme_minimal() +
    facet_wrap(~Condition, ncol=1)+
    theme(legend.position = 'none') +
    theme(axis.title.x = element_blank())
  MOL2_bp <- oligos_sum_simple %>%
    filter(MOL2==1)  %>%
    group_by(Animal, Layer, Anatomy, Condition) %>%
    select(`MOL2 Proportion`)%>%
    mutate(Proportion = `MOL2 Proportion`) %>%
    distinct() %>%
    ungroup() %>%
    ggplot(aes(x = Proportion, y = Layer)) +
    geom_violin() +
    geom_jitter(aes(color=Animal), alpha=0.6) +
    ggtitle('MOL2') +
    xlim(0,1) +
    theme_minimal() +
    facet_wrap(~Condition,ncol=1)+
    theme(legend.position = 'none',
          axis.title.y = element_blank(),
          axis.text.y = element_blank())
  MOL56_bp <- oligos_sum_simple %>%
    filter(`MOL5/6`==1)  %>%
    group_by(Animal, Layer, Anatomy, Condition) %>%
    select(`MOL5/6 Proportion`)%>%
    mutate(Proportion = `MOL5/6 Proportion`) %>%
    distinct() %>%
    ungroup() %>%
    ggplot(aes(x = Proportion, y = Layer)) +
    geom_violin() +
    geom_jitter(aes(color=Animal), alpha=0.6) +
    ggtitle('MOL5/6') +
    xlim(0,1) +
    theme_minimal() +
    facet_wrap(~Condition, ncol=1) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank())
  
  MOL1_bp+MOL2_bp+MOL56_bp
}
################################################################################
################### Load/Write Data ############################################
################################################################################

# Create measurement CSVs
oligo_path <- ACDTranslator(input_file_path = input_file_path, 
                           path = save_path,
                           extra = ID,
                           Egr2_thresh = Egr2_thresh,
                           Klk6_thresh = Klk6_thresh,
                           Ptgds_thresh = Ptgds_thresh)  

# Load the existing dfs for SC and Cortex
oligos_df <- read_oligos(oligo_path, cup_animals, dropGM = dropGM)

#oligos_df_existing <- read_oligos(existing_path, cup_animals,
#                            buffer=max(oligos_df$ID), dropGM=F) 

# Combine the cortex and SC dfs
#oligos_df <- rbind(oligos_df,oligos_df_existing)
# remove temporary df for memory management
#rm(oligos_df_existing)


################################################################################
######################### Summary Stats ########################################
################################################################################

# Compute summary
oligos_sum_simple <- oligos_df %>%
  select(Animal, Layer, Anatomy, Condition,
         Egr2_puncta, Klk6_puncta, Ptgds_puncta) %>%
  mutate(MOL1 = Egr2_puncta > Egr2_thresh,
         MOL2 = Klk6_puncta > Klk6_thresh,
         `MOL5/6` = Ptgds_puncta > Ptgds_thresh) %>% 
  SummaryMaker() %>%
  group_by(Condition) %>%
  distinct() %>%
  ungroup()


# Plot Violin
violin_grouped <- oligos_sum_simple %>%
  mutate(Layer = case_when(Layer=="L13" ~ "Cortex",
                           Layer=="L4" ~ "Cortex",
                           Layer=="L56" ~ "Cortex",
                           TRUE ~ Layer)) %>%
  oligo_violin_plotter() 

ggsave(file=paste0(save_path,"Violin_NOLayers.svg"), plot=violin_grouped, width=10, height=8)

violin <- oligo_violin_plotter(oligos_sum_simple)

ggsave(file=paste0(save_path,"Violin_Layers.svg"), plot=violin, width=10, height=8)

# Reformat summary to jmp
jmp_df <- oligos_sum_simple %>%
  mutate(`MOL1 Proportion` = `MOL1 Proportion` * MOL1, 
         `MOL2 Proportion` = `MOL2 Proportion` * MOL2, 
         `MOL5/6 Proportion` = `MOL5/6 Proportion` * `MOL5/6`) %>%
  select(Animal, Condition, Layer, `MOL1 Proportion`, `MOL2 Proportion`, 
         `MOL5/6 Proportion`, TotalCells, N_MOL1, N_MOL2, N_MOL56) %>%
  distinct() %>%
  group_by(Animal, Condition, Layer) %>%
  summarize(`MOL1 Proportion` = max(`MOL1 Proportion`), 
            `MOL2 Proportion` = max(`MOL2 Proportion`), 
            `MOL5/6 Proportion` = max(`MOL5/6 Proportion`),
            TotalCells = TotalCells,
            N_MOL1 = N_MOL1,
            N_MOL2 = N_MOL2,
            N_MOL56 = N_MOL56) %>%
  distinct()

jmp_df_grouped <- jmp_df %>%
  group_by(Animal, Condition) %>%
  mutate(Layer = case_when(Layer=="L13" ~ "Cortex",
                           Layer=="L4" ~ "Cortex",
                           Layer=="L56" ~ "Cortex",
                           TRUE ~ Layer)) %>%
  group_by(Animal, Condition, Layer) %>%
  mutate(TotalCells = sum(TotalCells),
         N_MOL1 = sum(N_MOL1),
         N_MOL2 = sum(N_MOL2),
         N_MOL56 = sum(N_MOL56),
         `MOL1 Proportion` = N_MOL1/TotalCells, 
         `MOL2 Proportion` = N_MOL2/TotalCells, 
         `MOL5/6 Proportion` = N_MOL56/TotalCells,
        ) %>%
  distinct()

# Save results
write_csv(jmp_df_grouped, paste0(save_path,"jmp_format_NOLAYERS.csv"))
write_csv(jmp_df, paste0(save_path,"jmp_format_LAYERS.csv"))
write.csv(oligos_sum_simple, paste0(save_path,'summary_byAnimal.csv'))