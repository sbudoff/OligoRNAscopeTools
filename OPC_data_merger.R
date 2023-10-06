library(dplyr)
library(ggplot2)

# load data without metadata
data <- read.csv("/media/sam/SamHDD/MT/single channel_roiMasks/Processed/Healthy_Cell.csv") 
data <- rbind(data, read.csv("/media/sam/SamHDD/MT/single channel_roiMasks/Processed/Healthy1_Cell.csv") )


# Subset to most interesting statistics
data <- data %>%
  mutate(ID = Parent_Nuclei,
         Condition = Metadata_Condition,
         Tissue = Metadata_ID,
         X = Location_Center_X,
         Y = Location_Center_Y,
         Klk6_Intensity = Intensity_IntegratedIntensity_Klk6_masked,
         Klk6_Count = Children_Klk6_puncta_Count,
         Egr2_Intensity = Intensity_IntegratedIntensity_Egr2_masked,
         Egr2_Count = Children_Egr2_puncta_Count,
         Ptgds_Intensity = Intensity_IntegratedIntensity_Ptgds_masked,
         Ptgds_Count = Children_Ptgds_puncta_Count,
         CorpusCallosum = Parent_CorpusCallosum,
         LayerI.III = Parent_LayerI_III,
         LayerIV = Parent_LayerIV,
         LayerV.VI = Parent_LayerV_VI,
         GM = LayerI.III + LayerIV + LayerV.VI,
         Location = case_when(CorpusCallosum==1 ~ "Corpus Callosum",
                              LayerI.III==1 ~ "Layer I-III",
                              LayerIV==1 ~ "Layer IV",
                              LayerV.VI==1 ~ "Layer V-VI",
                              TRUE~"Border")
         ) %>%
  select(ID, Condition, Tissue, X, Y, 
         Klk6_Intensity, Klk6_Count, 
         Egr2_Intensity,Egr2_Count,
         Ptgds_Intensity, Ptgds_Count,
         CorpusCallosum, GM, LayerI.III, LayerIV, LayerV.VI, Location) %>%
  group_by(Tissue) %>%
  mutate(X=X/max(X),
         Y=Y/max(Y)) %>%
  ungroup()


data %>%
  ggplot(aes(x=Location, y=Egr2_Intensity, color = Tissue)) +
  geom_violin()

Klk6_GM <- data %>%
  filter(GM==1) %>%
  pull(Klk6_Intensity)
Klk6_CC <- data %>%
  filter(CorpusCallosum==1) %>%
  pull(Klk6_Intensity)
t.test(Klk6_CC,Klk6_GM)
