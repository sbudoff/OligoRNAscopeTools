
library(tidyverse)

GridClusts <- function(mat,
                       max_k=15,
                       dist.methods = c("euclidean", "maximum", "manhattan", "canberra"),#, "DTW")
                       clust.methods = c("ward.D2", "ward.D", "single", "complete", 
                                         "average", "mcquitty", "median", "centroid")) {
  first = T
  for (dist.method in dist.methods){
    distMat <- dist(mat, method=dist.method)
    for (clust.method in clust.methods){
      hc <- hclust(distMat, method = clust.method)
      for (k in 2:max_k) {
        clusts <- cutree(hc, k)
        
        c.index <- clusterSim::index.C(distMat, clusts)
        d.index <- clValid::dunn(distMat, clusts)
        
        if (first) {
          partition_df <- data.frame(Cluster = clusts) %>%
            mutate(ID = 1:length(clusts),
                   dist.method = dist.method,
                   clust.method = clust.method,
                   k = k,
                   c.index = c.index,
                   d.index = d.index)
          first=F
        } else {
          temp <- data.frame(Cluster = clusts) %>%
            mutate(ID = 1:length(clusts),
                   dist.method = dist.method,
                   clust.method = clust.method,
                   k = k,
                   c.index = c.index, 
                   d.index = d.index)
          partition_df <- rbind(partition_df, temp)
          
        }
      }
    }
  }
  return(partition_df)
}


ACDTranslator <- function(input_file_path, path, extra="",
                          GFP_thresh = 0,
                          ReturnSummary=TRUE,
                          SaveResults=TRUE,
                          ComplexFileName=TRUE) {
  data <- read_csv(input_file_path  ) 
  
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
  
  # Convert columns into more human readable names
  data <- data %>%
    filter(!is.na(ID)) %>%
    mutate(GFP = `Subcellular: Channel 2: Num single spots`,
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
           Location = Parent
           ) %>%
    select(ID, Image, Location, GFP, Cell_Area, X, Y, n_Egr2, intesity_Egr2, area_Egr2,
           n_Klk6, intesity_Klk6, area_Klk6, n_Ptgds, intesity_Ptgds, area_Ptgds) %>%
    mutate(X = ifelse(is.na(Cell_Area), NA, X),
           Y = ifelse(is.na(Cell_Area), NA, Y),
           Location = ifelse(Location=="PathCellObject", NA, Location),
           GFP = ifelse(is.na(Cell_Area), NA, GFP)) 
  
  # Extract the image locations for later
  Image_IDs <- data %>%
    select(ID, Image) %>%
    drop_na()
  
  # remove character columns
  data <- data %>%
    mutate(Location = str_extract(Location, "^[^_]+"),
           GM = ifelse(Location == "GM", 1, 0),
           CC = ifelse(Location == "CC", 1, 0) + ifelse(Location == "WM", 1, 0),
           L13 = ifelse(Location == "L13", 1, 0) + ifelse(Location == "L23", 1, 0),
           L4 = ifelse(Location == "L4", 1, 0),
           L56 = ifelse(Location == "L56", 1, 0)) %>%
    select(-c(Image, Location))
  
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
    mutate(Egr2=Egr2/max(Egr2),
           Klk6=Klk6/max(Klk6),
           Ptgds=Ptgds/max(Ptgds),
           Egr2_bin = Egr2 > 0,
           Kl6_bin = Klk6 > 0,
           Ptgds_bin = Ptgds > 0,
           Egr2_pp = Egr2/n_Egr2,
           Klk6_pp = Klk6/n_Klk6,
           Ptgds_pp = Ptgds/n_Ptgds)
  
  flattened_data[is.na(flattened_data)] <- 0
  
  # Put image IDs back in
  flattened_data <- left_join(Image_IDs, flattened_data) %>%
    unique()
  
  # Empirical threshold for OPCs
  # GFP_thresh = 190#300
  
  # Simplest classification scheme, which normalized value is greatest in a given cell
  if (ComplexFileName){
    OPC_data <- flattened_data %>%
      filter(GFP > GFP_thresh) %>%
      mutate(Classification = case_when(
        Egr2 > Klk6 & Egr2 > Ptgds ~ "Egr2",
        Klk6 > Egr2 & Klk6 > Ptgds ~ "Klk6",
        Ptgds > Egr2 & Ptgds > Klk6 ~ "Ptgds",
        TRUE ~ "None"
      ),
      filename = str_extract(Image, "_(.*?)\\."),
      filename = str_remove(filename, ".*_"),
      Animal = str_extract(filename, "\\d+"),
      Animal = factor(Animal))
  } else {
    OPC_data <- flattened_data %>%
      filter(GFP > GFP_thresh) %>%
      mutate(Classification = case_when(
        Egr2 > Klk6 & Egr2 > Ptgds ~ "Egr2",
        Klk6 > Egr2 & Klk6 > Ptgds ~ "Klk6",
        Ptgds > Egr2 & Ptgds > Klk6 ~ "Ptgds",
        TRUE ~ "None"
      ),
      Animal = str_extract(Image, "\\d+"),
      Animal = factor(Animal))
  }
  
  # Simple summary table counting Classified OPCs in a given region
  OPC_perRegion_perAnimal <- OPC_data %>%
    group_by(Classification, Animal, GM, CC, L13, L4, L56) %>%
    summarise(N = n() ) %>%
    ungroup() %>%
    group_by(Animal, GM, CC, L13, L4, L56) %>%
    mutate(Proportion = N / sum(N))
  
  OPC_perRegion <- OPC_perRegion_perAnimal %>%
    group_by(Classification, GM, CC, L13, L4, L56) %>%
    summarise(Mean = mean(Proportion) ,
              SD = sd(Proportion)) %>%
    ungroup() 
    
  if (SaveResults) {
    # Save outputs
    write.csv(flattened_data, file = paste0(path,extra,"measurements_all.csv"))
    write.csv(OPC_data, file = paste0(path,extra,"measurements_Oligo.csv"))
    write.csv(OPC_perRegion_perAnimal, file = paste0(path,extra,"summary_Oligo_by_Animal.csv"))
    write.csv(OPC_perRegion, file = paste0(path,extra,"summary_Oligo.csv"))
  }
  
  if (ReturnSummary) {
    return(OPC_perRegion_perAnimal)
  } else {
    out <- list()
    out[['Summary']] <- OPC_perRegion_perAnimal
    out[['Summary Simple']] <- OPC_perRegion
    out[['All']] <- flattened_data
    out[['Oligos']] <- OPC_data
    return(out)
  }
  
}


# OPC_perRegion_perAnimal_SC <- ACDTranslator(input_file_path = '/home/sam/Downloads/measurements.csv', 
#                                             path = "/home/sam/Downloads/",
#                                             extra = "SC", GFP_thresh = 300,
#                                             ReturnSummary = T)  %>%
#   mutate(Anatomy="SC")


OPC_perRegion_list <- ACDTranslator(input_file_path = '/media/sam/SamHDD/MT/measurements_Run1.csv', 
                                            path = "/home/sam/Downloads/",
                                            extra = "Cortex", GFP_thresh = 0,
                                            ReturnSummary = F, ComplexFileName = F) 
OPC_perRegion_perAnimal_Cortex <- OPC_perRegion_list[['Summary']] 

OPC_perRegion_perAnimal_Cortex <- OPC_perRegion_perAnimal_Cortex%>%
  mutate(Anatomy="Cortex")

summary_df_by_animal <- OPC_perRegion_perAnimal_Cortex %>%#rbind(OPC_perRegion_perAnimal_Cortex, OPC_perRegion_perAnimal_SC) %>%
  gather(Layer, Val, -c(Classification, Animal, Anatomy, N, Proportion))%>%
  filter(Val !=0 ) %>%
  select(-Val) %>%
  mutate(Classification = factor(Classification),
         Animal = factor(Animal),
         Anatomy = factor(Anatomy),
         Layer = factor(Layer))

summary_df <- summary_df_by_animal %>%
  group_by(Classification, Layer, Anatomy) %>%
  summarise(Mean = mean(Proportion) ,
            SD = sd(Proportion)) %>%
  ungroup() 

# Save outputs
write.csv(summary_df, file = paste0("/home/sam/Downloads/","OPC_finalSummary.csv"))

# menage <- aov(Proportion ~ ., data = summary_df_by_animal)
# summary(menage)

summary_df_by_animal_cort <- summary_df_by_animal %>%
  filter(Anatomy=="Cortex") %>%
  filter(Layer != "GM")
twoWay <- aov(Proportion ~ Classification + Layer + Classification * Layer, data = summary_df_by_animal_cort)
summary(twoWay)

summary_df_by_animal_cort %>%
  ggplot(aes(0, y = Proportion, color = Classification, fill = Classification)) +
  geom_violin() +
  ggtitle("Distribution of OPCs Parsed By RNAscope Classification") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  facet_grid(~Layer)

Oligo_df <- OPC_perRegion_list[['Oligos']]

Oligo_mat <- Oligo_df %>%
  mutate(n_Egr2 = n_Egr2/max(n_Egr2),
         n_Klk6 = n_Klk6/max(n_Klk6),
         n_Ptgds = n_Ptgds/ max(n_Ptgds)) %>%
  select(Egr2, Klk6, Ptgds) %>% #n_Egr2, n_Klk6, n_Ptgds,
  as.matrix()

# clust_grid <- GridClusts(Oligo_mat,
#                          max_k=10,
#                          dist.methods = c("euclidean"),#, "maximum", "manhattan"),
#                          clust.methods = c("ward.D2", "ward.D", "single", "complete", 
#                                             "average", "mcquitty", "median", "centroid"))  
# 
# 
# clust_grid %>%
#   select(-ID,-Cluster) %>%
#   distinct() %>%
#   gather(key='Index', value="Value", -c(k,dist.method,clust.method)) %>%
#   ggplot(aes(x = k, y=Value, color=clust.method)) +
#   geom_line() +
#   facet_grid(~Index*dist.method)

# 
# test <- OPC_perRegion_list[['All']]
# GFP_thresh=0
# 
# OPC_data <- test %>%
#   filter(GFP > GFP_thresh) %>%
#   mutate(Classification = case_when(
#     Egr2 > Klk6 & Egr2 > Ptgds ~ "Egr2",
#     Klk6 > Egr2 & Klk6 > Ptgds ~ "Klk6",
#     Ptgds > Egr2 & Ptgds > Klk6 ~ "Ptgds",
#     TRUE ~ "None"
#   ),
#   Animal = as.integer(str_extract(Image, "\\d+")),)
# 
# 
# 
# OPC_perRegion_perAnimal <- OPC_data %>%
#   group_by(Classification, Animal, GM) %>%
#   summarise(N = n() ) %>%
#   ungroup() %>%
#   group_by(Animal, GM) %>%
#   mutate(Proportion = N / sum(N))
# 
# OPC_perRegion <- OPC_perRegion_perAnimal %>%
#   group_by(Classification, GM) %>%
#   summarise(Mean = mean(Proportion) ,
#             SD = sd(Proportion)) %>%
  ungroup() 

  
 test <- summary_df_by_animal_cort %>%
    group_by(Animal,Layer) %>%
    summarise(sum(Proportion))
  