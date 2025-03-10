##***************************************##
##   R script to generate main Figure 3  ##
##***************************************##
## Project: project_aging
## Sarah Hoelzl
## Last modification 08.2024
## Creation: 04.2022 by SH

library(gplots)
library(tidyverse)
library(ggplot2)
library(data.table)
library(dplyr)
library(plyr)
library(RColorBrewer)
library(gdata)
library(ggpubr)
library(ggrepel)
library(readxl)
library(gridExtra)
library(ggplotify)

######------ Set environment  ------###### 
input <- "/Users/shoelzl/Desktop/mnt_new2/ge43sik2/WP1_aging_SH/manuscript/EscaperPaper/20241216_nature_aging_revision/revised_supplemetary_tables/"

######------ Read in snRNA-seq analysis table  ------###### 
snRNA_analysis <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "n snRNASeq analysis")

######------ Read in allelic ratios  ------###### 
allelic_ratios <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "o snPseudobulk AllelicRatio")

# set and apply min_total_reads cutoff
min=20
allelic_ratios <- allelic_ratios %>% dplyr::filter(total_reads >= min)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------#


######---- Figure 3c: Cell population ----###### 
#Count cell populations
population <- ddply(snRNA_analysis,c("age","cluster"),nrow)
population <- population %>% dplyr::rename(count =V1)

# Compute percentages 
summary <- population %>%
  dplyr::group_by(age) %>%
  dplyr::summarize(sum = sum(count))
population <- population%>%left_join(summary,by="age")
population <- population %>%  dplyr::mutate(fraction = (count/sum)*100)

# set order
population$cluster <- factor(population$cluster,
                              levels = c("Glial","Tcell","FB2","LEC","Adipo","SMC","Pericyte","Epicard","Macro","atrCM","EC2","EC1","FB1","ventCM"),ordered = TRUE)
population$age <- factor(population$age,
                             levels = c("Aged", "Adult"),ordered = TRUE)

#set colors
cluster_colors <- c("#EC68A3", "#CC6FA9", "#459AD5", "#11B6EB", "#9188C0", "#2EB7BE", "#3EAD58", "#AD7AB3", "#36B28F", "#59B031", "#9AA921", "#C49B05", "#E38903", "#EE766F")

#plot
Figure_3c <- ggplot(population, aes(x = fraction, y = age, fill = cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values= cluster_colors) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Figure 3c",
       x = "Cell composition [%]",
       y = "",
       fill = "")
  
######---- Figure 3d: Fractions of Xa states ----###### 

#Count Xa states
Xa_states <- ddply(snRNA_analysis,c("age","Xa"),nrow)
Xa_states <- Xa_states %>% dplyr::rename(count =V1)
Xa_states <- Xa_states %>% dplyr::filter(Xa != "not_assignable")

# Compute percentages 
summary <- Xa_states %>%
  dplyr::group_by(age) %>%
  dplyr::summarize(sum = sum(count))
Xa_states <- Xa_states%>%left_join(summary,by="age")
Xa_states <- Xa_states %>%  dplyr::mutate(fraction = (count/sum)*100)

# set order
Xa_states$Xa <- factor(Xa_states$Xa,
                             levels = c("CAST","Bl6","bial"),ordered = TRUE)
Xa_states$age <- factor(Xa_states$age,
                         levels = c("Adult", "Aged"),ordered = TRUE)

#set colors
Xa_colors <- c("#CC9B24", "#000000", "#DADADA")

#plot
Figure_3d <- ggplot(Xa_states, aes(x = Xa, y = fraction, fill = Xa)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(~age)+
  scale_fill_manual(values= Xa_colors) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Figure 3d",
       x = "",
       y = "Fractions of Xa states [%]",
       fill = "")

######---- Figure 3f: Age-specific escape ----###### 

#save age-specific escape genes detected in the whole heart
age_spec_esc <- c("Shroom4","Tspan7", "Sh3kbp1", "Med14", "Kctd12b", "2210013O21Rik")

#Filter dataframe for age-specific escape genes and prepare for stripchart
Stripchart <- allelic_ratios %>%  dplyr::filter(name %in% age_spec_esc)

# Filter for genes expressed in both time points
count_table <- ddply(Stripchart,c("name","organ"),nrow)
colnames(count_table)[3] <- "timepoint_count"
Stripchart <- merge(Stripchart,count_table)
Stripchart <- Stripchart[Stripchart$timepoint_count == 2,]

# set order
Stripchart$organ <- factor(Stripchart$organ, levels = c("CM","EC","FB","MP","Pericyte","SMC"))

#plot

Figure_3f <- ggplot(Stripchart, aes(x = organ, y = allelic_ratio))+
  geom_point(aes(col=organ))+
  facet_wrap(~age) +
  scale_y_continuous(limits= c(0.6,1), breaks = seq(0.6, 1, by = 0.1)) +
  geom_text(aes(label = ifelse(allelic_ratio<=0.90 & age == "Aged", name ,"")), size=2, na.rm = TRUE, hjust = -0.3) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "#8B1913") + 
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust=0.1)) +
  labs(title = "Figure 3f", 
       x = "", 
       y = "Allelic ratio") 


######---- Plot Figure 3 ----###### 
grid.arrange(Figure_3c, Figure_3d, Figure_3f,  
             layout_matrix = rbind(c( 1 , 1),
                                   c(2, 3)))
