################################################################################
########################  Hyperparameters  #####################################
################################################################################

input_file_path <- "C:/Users/Sam/Downloads/Cortexmeasurements_all_allmice.csv"
save_path <- "C:/Users/Sam/Downloads/"
extra <- ""
GFP_thresh <- 0
dropGM <- F
layers <- T

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
read_cells <- function(path, GFP_thresh, buffer=0, dropGM=T, layers=F) {
  out <- read_csv(file = path) %>%
    select(-`...1`) %>%
    group_by(Image) %>%
    mutate(Oligodendrocyte = GFP > GFP_thresh,
           Animal = factor(Animal),
           ID = ID+buffer) %>%
    ungroup()
  if (layers) {
    out <- out %>%
      select(Oligodendrocyte, Animal, L13, L4, L56, CC)
    
  } else {
    out <- out %>%
      mutate(`White Matter` = CC,
             `Grey Matter` = (CC == 0) + 0) %>%
      select(Oligodendrocyte, Animal, `Grey Matter`,`White Matter`)
  }
           
  return(out)
}

################################################################################
################### Load/Write Data ############################################
################################################################################

data <- read_cells(input_file_path, 
                   GFP_thresh = GFP_thresh, dropGM = dropGM, layers = layers) 

################################################################################
######################### Summary Stats ########################################
################################################################################

if (layers) {
  Layer_data <- data %>%
    gather("Layer", "val", -c(Oligodendrocyte, Animal)) %>%
    filter(val==1) %>%
    select(-val) %>%
    group_by(Animal, Layer) %>%
    summarize(`Percent Oligodendrocytes` = sum(Oligodendrocyte)/n(),
              Oligodendrocytes = sum(Oligodendrocyte),
              Cells = n()) %>%
    ungroup()
  #Save data
  write_csv(Layer_data, paste0(save_path,extra,"OligoPercentANDCounts_LAYERS.csv"))
  
} else {
  GMWM_data_onehot <- data %>%
    select(Oligodendrocyte, Animal, `Grey Matter`, `White Matter`) %>%
    group_by(Animal, `Grey Matter`, `White Matter`) %>%
    summarize(`Percent Oligodendrocytes` = sum(Oligodendrocyte)/n()) %>%
    ungroup()
  
  
  GMWM_data <- GMWM_data_onehot %>% 
    gather("Anatomy", "val", -c(Animal, `Percent Oligodendrocytes`)) %>%
    filter(val == 1 ) %>%
    select(-val)
  
  GMWM_data_onehot_raw <- data %>%
    select(Oligodendrocyte, Animal, `Grey Matter`, `White Matter`) %>%
    group_by(Animal, `Grey Matter`, `White Matter`) %>%
    summarize(Oligodendrocytes = sum(Oligodendrocyte),
              Cells = n()) %>%
    ungroup()
  
  
  GMWM_data_raw <- GMWM_data_onehot_raw %>% 
    gather("Anatomy", "val", -c(Animal, Oligodendrocytes, Cells)) %>%
    filter(val == 1 ) %>%
    select(-val)
  
  
  GMWM_bar <- GMWM_data %>%
    ggplot(aes(x = Animal, y = `Percent Oligodendrocytes`, fill = Anatomy)) +
    geom_bar(stat='identity', position='dodge') +
    ggtitle("Oligodendrocyte Percent by Area and Location") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual("Anatomy", values = c("Grey Matter" = "grey40", "White Matter" = "gray87"))
  
  # Save results
  ggsave(file=paste0(save_path,extra,"GMWM_bar.svg"), plot=GMWM_bar, width=10, height=8)
  write_csv(GMWM_data, paste0(save_path,extra,"OligoPercent_GMWM_NOLAYERS.csv"))
  write_csv(GMWM_data_onehot, paste0(save_path,extra,"OligoPercent_GMWM_NOLAYERS_1hot.csv"))
  write_csv(GMWM_data_raw, paste0(save_path,extra,"OligoPercent_GMWM_NOLAYERS_raw.csv"))
  write_csv(GMWM_data_onehot_raw, paste0(save_path,extra,"OligoPercent_GMWM_NOLAYERS_1hot_raw.csv"))
  
}

