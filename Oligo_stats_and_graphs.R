library(tidyverse)
library(patchwork)

# Load data and format
# Set the working directory to the script's directory
root <- normalizePath(commandArgs()[1], mustWork = FALSE)
root <- dirname(root)
setwd(dirname(root))

# Print the current working directory (root directory)
setwd('/home/sam/cellprofiler_sam/CellProfiler/')
cat("Root directory:", getwd(), "\n")

cup_animals <- c(6126, 6128,6142)

summary_df <- read_csv(file = "Oligo_finalSummary.csv") %>%
  select(-`...1`)


read_oligos <- function(path, cup_animals,buffer=0, dropGM=T) {
  out <- read_csv(file = path) %>%
    select(-`...1`) %>%
    mutate(Anatomy="Cortex") %>%
    mutate(Healthy = case_when(Animal %in% cup_animals ~ 0,
                               TRUE ~ 1),
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

oligos_df <- read_oligos("Cortexmeasurements_Oligo.csv", cup_animals)

oligos_df_sc <- read_oligos("/media/sam/STOCKTON/073123_CSVs/SCmeasurements_Oligo.csv", cup_animals,
                            buffer=max(oligos_df$ID), dropGM=F) 

oligos_df <- rbind(oligos_df,oligos_df_sc)

rm(oligos_df_sc)


################################################################################
######################### Summary Stats ########################################

SummaryMaker <- function(oligos_df) {

  oligos_sum <- oligos_df %>%
    group_by(Animal, Layer, Anatomy, Healthy) %>%
    mutate(TotalCells = n()) %>%
    ungroup() %>%
    group_by(MOL1, Animal, Layer, Anatomy, Healthy) %>%
    mutate(N_mol1 = n()) %>%
    ungroup() %>%
    group_by(MOL2, Animal, Layer, Anatomy, Healthy) %>%
    mutate(N_mol2 = n()) %>%
    ungroup() %>%
    group_by(`MOL5/6`, Animal, Layer, Anatomy, Healthy) %>%
    mutate(N_mol56 = n()) %>%
    ungroup() %>%
    mutate(Healthy = case_when(Healthy == 0 ~ "Cuprizone",
                     T ~ "Healthy")) %>%
    ungroup() %>%
    distinct() 
  
  # Define the desired order of the categories
  layer_order <- c("L13", "L4", "L56", "CC", "GM", "WM")
  
  # Convert the Category column to a factor with the desired order
  summary_df$Layer <- factor(summary_df$Layer, levels = rev(layer_order))
  oligos_sum$Layer <- factor(oligos_sum$Layer, levels = rev(layer_order))
  
  oligos_sum_simple <- oligos_sum  %>%
    group_by(MOL1, MOL2, `MOL5/6`, Animal, Layer, Anatomy, Healthy) %>%
    mutate(`mean N_1` = mean(N_mol1),
           `mean N_2` = mean(N_mol2),
           `mean N_5/6` = mean(N_mol56)) %>%
    distinct() %>%
    ungroup()
  
  MOL1_summary <- oligos_sum_simple %>%
    select(MOL1, Animal, Layer, Anatomy, Healthy, `mean N_1`) %>%
    group_by(MOL1, Animal, Layer, Anatomy, Healthy) %>%
    summarize(Cells = sum(`mean N_1`)) %>%
    ungroup() %>%
    group_by(Animal, Layer, Anatomy, Healthy) %>%
    mutate(`MOL1 Proportion` = Cells/sum(Cells)) %>%
    select(-Cells) %>%
    ungroup()
  
  MOL2_summary <- oligos_sum_simple %>%
    select(MOL2, Animal, Layer, Anatomy, Healthy, `mean N_2`) %>%
    group_by(MOL2, Animal, Layer, Anatomy, Healthy) %>%
    summarize(Cells = sum(`mean N_2`)) %>%
    ungroup() %>%
    group_by(Animal, Layer, Anatomy, Healthy) %>%
    mutate(`MOL2 Proportion` = Cells/sum(Cells)) %>%
    select(-Cells) %>%
    ungroup()
  
  MOL56_summary <- oligos_sum_simple %>%
    select(`MOL5/6`, Animal, Layer, Anatomy, Healthy, `mean N_5/6`) %>%
    group_by(`MOL5/6`, Animal, Layer, Anatomy, Healthy) %>%
    summarize(Cells = sum(`mean N_5/6`)) %>%
    ungroup() %>%
    group_by(Animal, Layer, Anatomy, Healthy) %>%
    mutate(`MOL5/6 Proportion` = Cells/sum(Cells)) %>%
    select(-Cells) %>%
    ungroup()
    
  oligos_sum_simple <- left_join(oligos_sum_simple, MOL1_summary, by = c("MOL1", "Animal", "Layer", "Anatomy", "Healthy"))
  oligos_sum_simple <- left_join(oligos_sum_simple, MOL2_summary, by = c("MOL2", "Animal", "Layer", "Anatomy", "Healthy"))
  oligos_sum_simple <- left_join(oligos_sum_simple, MOL56_summary, by = c("MOL5/6", "Animal", "Layer", "Anatomy", "Healthy"))

  return(oligos_sum_simple)
}

oligo_violin_plotter <- function(oligos_sum_simple) {
  MOL1_bp <- oligos_sum_simple %>%
    filter(MOL1==1) %>%
    group_by(Animal, Layer, Anatomy, Healthy) %>%
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
    facet_wrap(~Healthy, ncol=1)+
    theme(legend.position = 'none') +
    theme(axis.title.x = element_blank())
  MOL2_bp <- oligos_sum_simple %>%
    filter(MOL2==1)  %>%
    group_by(Animal, Layer, Anatomy, Healthy) %>%
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
    facet_wrap(~Healthy,ncol=1)+
    theme(legend.position = 'none',
          axis.title.y = element_blank(),
          axis.text.y = element_blank())
  MOL56_bp <- oligos_sum_simple %>%
    filter(`MOL5/6`==1)  %>%
  group_by(Animal, Layer, Anatomy, Healthy) %>%
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
    facet_wrap(~Healthy, ncol=1) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank())
  
  MOL1_bp+MOL2_bp+MOL56_bp
}

################################################################################
################################################################################
################################################################################

oligos_sum_simple <- oligos_df %>%
  select(Animal, Layer, Anatomy, Healthy,
         Egr2_puncta, Klk6_puncta, Ptgds_puncta) %>%
  mutate(MOL1 = Egr2_puncta > 3,
         MOL2 = Klk6_puncta > 4,
         `MOL5/6` = Ptgds_puncta > 0) %>% 
  SummaryMaker() %>%
  group_by(Healthy) %>%
  distinct() %>%
  ungroup()

oligos_sum_simple %>%
  mutate(Layer = case_when(Layer=="L13" ~ "Cortex",
                           Layer=="L4" ~ "Cortex",
                           Layer=="L56" ~ "Cortex",
                           TRUE ~ Layer)) %>%
  oligo_violin_plotter() 

jmp_df <- oligos_sum_simple %>%
  mutate(Layer = case_when(Layer=="L13" ~ "Cortex",
                           Layer=="L4" ~ "Cortex",
                           Layer=="L56" ~ "Cortex",
                           TRUE ~ Layer),
         `MOL1 Proportion` = `MOL1 Proportion` * MOL1, 
         `MOL2 Proportion` = `MOL2 Proportion` * MOL2, 
         `MOL5/6 Proportion` = `MOL5/6 Proportion` * `MOL5/6`) %>%
  select(Animal, Healthy, Layer, `MOL1 Proportion`, `MOL2 Proportion`, `MOL5/6 Proportion`) %>%
  distinct() %>%
  group_by(Animal, Healthy, Layer) %>%
  summarize(`MOL1 Proportion` = max(`MOL1 Proportion`), 
            `MOL2 Proportion` = max(`MOL2 Proportion`), 
            `MOL5/6 Proportion` = max(`MOL5/6 Proportion`))

write_csv(jmp_df, '/media/sam/STOCKTON/073123_CSVs/073123_jmp_format_NOLAYERS.csv')
write.csv(oligos_sum_simple, '/media/sam/STOCKTON/073123_CSVs/073123_summary_byAnimal.csv')


thesh_list <- list()
for (test_thresh in 1:150) {
  
  thesh_list[[test_thresh]] <- oligos_df %>%
    select(Animal, Layer, Anatomy, Healthy,
           Egr2_puncta, Klk6_puncta, Ptgds_puncta) %>%
    mutate(MOL1 = Egr2_puncta > test_thresh,
           MOL2 = Klk6_puncta > test_thresh,
           `MOL5/6` = Ptgds_puncta > test_thresh) %>% 
    SummaryMaker() %>%
    group_by(Healthy) %>%
    distinct() %>%
    summarize(`MOL1 Proportion` = mean(`MOL1 Proportion`),
              `MOL2 Proportion` = mean(`MOL2 Proportion`),
              `MOL5/6 Proportion` = mean(`MOL5/6 Proportion`),
              Thresh = test_thresh)
}


mol_1_trace <- do.call(rbind, thesh_list)  %>%
  ggplot(aes(x=Thresh, y = `MOL1 Proportion`)) +
  geom_line() +
  theme_minimal() +
  theme(axis.title.x  = element_blank()) +
  facet_wrap(~Healthy, ncol=2, scales = "free_y")

mol_2_trace <- do.call(rbind, thesh_list)  %>%
  ggplot(aes(x=Thresh, y = `MOL2 Proportion`)) +
  geom_line() +
  theme_minimal() +
  theme(axis.title.x  = element_blank()) +
  facet_wrap(~Healthy, ncol=2, scales = "free_y")

mol_56_trace <- do.call(rbind, thesh_list)  %>%
  ggplot(aes(x=Thresh, y = `MOL5/6 Proportion`)) +
  geom_line() +
  theme_minimal() +
  xlim(0,50) +
  theme(axis.title.x  = element_blank()) +
  facet_wrap(~Healthy, ncol=2, scales = "free_y")

mol_1_trace / mol_2_trace / mol_56_trace

oligos_df %>%
  select(Animal, Layer, Anatomy, Healthy,
         Egr2_puncta, Klk6_puncta, Ptgds_puncta) %>%
  group_by(Animal, Layer, Anatomy, Healthy) %>%
  ggplot(aes( y = Klk6_puncta)) +
  geom_histogram() +
  facet_wrap(~Healthy*Layer)


Floridia_thresh <- oligos_df %>%
  select(Animal, Layer, Anatomy, Healthy,
         Egr2_puncta, Klk6_puncta, Ptgds_puncta) %>%
  mutate(MOL1 = Egr2_puncta > 3,
         MOL2 = Klk6_puncta > 4,
         `MOL5/6` = Ptgds_puncta > 12) %>% 
  SummaryMaker() %>%
  group_by(Healthy) %>%
  distinct() 


oligo_violin_plotter(Floridia_thresh) 

CU_thresh <- oligos_df %>%
  select(Animal, Layer, Anatomy, Healthy,
         Egr2_puncta, Klk6_puncta, Ptgds_puncta) %>%
  mutate(MOL1 = Egr2_puncta > 3,
         MOL2 = Klk6_puncta > 5,
         `MOL5/6` = Ptgds_puncta > 0) %>% 
  SummaryMaker() %>%
  group_by(Healthy) %>%
  distinct() 


oligo_violin_plotter(CU_thresh) 







  



oligos_sum2 <- oligos_sum %>%
  group_by(MOL1, MOL2, `MOL5/6`, Layer, Anatomy, Healthy) %>%
  summarise(Mean = mean(Proportion) ,
            SD = sd(Proportion)) %>%
  ungroup() 

# oligo_sum %>%
#   filter(Anatomy=='Cortex') %>%
#   select(Animal, Layer, Healthy, MOL1, MOL2, `MOL5/6`) %>%
#   gather(ID, Val, )
# 
# oligos

MOL1_bp <- oligos_sum %>%
  ungroup() %>%
  filter(MOL1==1) %>%
  group_by(Animal, Layer, Anatomy, Healthy) %>%
  mutate(test = mean(Proportion))

test <- oligos_df %>%
  filter(Healthy == 0 ) %>%
  select(Animal, Healthy, n_Ptgds, Layer)
