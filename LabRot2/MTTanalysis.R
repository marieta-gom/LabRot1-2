# Load necessary libraries
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

# Load the data from the CSV file
data <- read_csv("C:\\Users\\Lenovo\\Desktop\\MolBio\\LabRotation2\\all_data_MTT.csv", show_col_types = F)  

# Preliminary calculations: Mean and standard deviation for each sample group
grouped_replicates <- data %>%
  group_by(Sample) %>%
    summarise(
      Mean = mean(`Absorbance(554nm)`, na.rm = TRUE),  
      SD = sd(`Absorbance(554nm)`, na.rm = TRUE),  
      #PerctUncert = SD/Mean*100  
)

# Processing of each replicate: Removal of the blank (no_cell) and normalization (no_gel)
  ## The desired value is (mean_of_the_sample - no_cell_condition)/no_gel_condition, for each plate
  
processed_replicates <- grouped_replicates %>%
  separate(Sample, c("Fluorophore", "Modification", "Amount", "Replicate")) %>%
  transform(Amount = as.numeric(Amount)) %>%
  transform(Replicate = as.numeric(Replicate)) %>%
  group_by(Fluorophore, Replicate) %>%
  mutate(Mean_blank = Mean - Mean[Modification == "NoCells"]) %>%
  mutate(SD_blank = sqrt(SD^2 + SD[Modification == "NoCells"]^2)) %>% # Pitagoras form to calculate the SD of the new calculated variable
  #mutate(PerctUncert_blank = SD_blank/Mean_blank*100) %>%
  mutate(Mean_norm = Mean_blank/Mean_blank[Modification == "NoPolymer"]) %>%
  #mutate(PerctUncert_norm = sqrt(PerctUncert_blank^2 + PerctUncert_blank[Modification == "NoPolymer"]^2)) %>%
  mutate(SD_norm = sqrt(SD_blank^2 + (SD_blank[Modification == "NoPolymer"]/(-Mean_blank[Modification == "NoPolymer"]^-2))^2)) %>%
  ungroup()

## Optionally, write the results to a new CSV file
#write.csv(processed_replicates, "C:\\Users\\Lenovo\\Desktop\\MolBio\\LabRotation2\\output_processed_replicates.csv", row.names = FALSE)


# Final data analysis: Integration of the replicates

final_data <- processed_replicates %>%
  mutate(Condition = paste(processed_replicates$Fluorophore, 
                      processed_replicates$Modification, 
                      processed_replicates$Amount, sep = "_")) %>%
  subset(select = -c(5:10)) %>%
  group_by(Condition) %>%
  mutate(Mean_final = mean(Mean_norm)) %>%
  #mutate(PerctUncert_final = sqrt(PerctUncert_norm[Replicate == "1"]^2 + PerctUncert_norm[Replicate == "2"]^2 + PerctUncert_norm[Replicate == "3"]^2)) %>%
  mutate(SD_final = sqrt(SD_norm[Replicate == "1"]^2 + SD_norm[Replicate == "2"]^2 + SD_norm[Replicate == "3"]^2)) %>%
  distinct(Condition, .keep_all = T) %>%
  subset(select = -c(4:7)) %>%
  ungroup()

final_data_fluor <- subset(final_data, Fluorophore == "Fluorescein") %>%
  subset(Amount != '0') %>% #Eliminate the control conditions (do not give any extra information and make the plot more complex)
  arrange(Amount)
final_data_TAMRA <- subset(final_data, Fluorophore == "TAMRA") %>%
  final_data_TAMRA[final_data_TAMRA$Amount != '0',] %>% #Eliminate the control conditions (do not give any extra information and make the plot more complex)
  arrange(Amount)
final_data_Cy5 <- subset(final_data, Fluorophore == "Cy5") %>%
  final_data_Cy5[final_data_Cy5$Amount != '0',] %>% #Eliminate the control conditions (do not give any extra information and make the plot more complex)
  arrange(Amount)
  

#Plot the results


  group_by(Fluorophore) %>%
  ggplot(aes(x = factor(Amount, level=c('0', '1', '2', '10', '20', '100', '200', '1000')), 
             y = Mean_final, 
             fill = Modification)) +
  geom_point(stat="identity", position = "dodge") +
  geom_line(aes(x = Amount, y = Mean_final) +
  scale_fill_manual(values = c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")) + #colorblind palette
  xlab("Polymer amount [mg/well]") +
  ylab("Mean absorbance [554 nm]") +
  ggtitle("MTT assay") +
  theme(plot.title = element_text(size=19)) +
  theme_classic() #+
  #geom_pointrange(aes(x = factor(Amount, level=c('0', '1', '2', '10', '20', '100', '200', '1000')), y = Mean_final, ymin=Mean_final-SD_final, ymax=Mean_final+SD_final), colour="black")

