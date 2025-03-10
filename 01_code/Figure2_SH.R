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
library(ComplexHeatmap)
library(RColorBrewer)
library(gdata)
library(ggpubr)
library(ggrepel)
library(eulerr)
library(alluvial)
library(readxl)
library(gridExtra)
library(ggplotify)

######------ Set environment  ------###### 
input <- "/Users/shoelzl/Desktop/mnt_new2/ge43sik2/WP1_aging_SH/manuscript/EscaperPaper/20241216_nature_aging_revision/revised_supplemetary_tables/"
#output <- "/Users/shoelzl/Desktop/mnt_new2/ge43sik2/WP1_aging_SH/10_script_submission/Figure_03/"

######------ Read in Allelic Ratios of all time points ------###### 
allelic_ratios_Embryonic <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "f AllelicRatio_Embryonic")
allelic_ratios_Young <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "h AllelicRatio_Young")
allelic_ratios_Adult <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "b AllelicRatio_Adult")
allelic_ratios_Aged <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "j AllelicRatio_Aged")

# set and apply min_total_reads cutoff
min=20

#apply cutoff and extract median from all allelic ratio tables
table_names <- c("allelic_ratios_Embryonic", "allelic_ratios_Young", "allelic_ratios_Adult", "allelic_ratios_Aged")

# Loop through the table names
for (name in table_names) {
  table <- get(name)
#apply total_read cutoff
  table <- table %>% dplyr::filter(min_total_reads >= min)
# extract median data
  table$sample <- gsub("_1","",table$sample)
  table$sample <- gsub("_2","",table$sample)
  table$sample <- gsub("_3","",table$sample)
  table <- table %>% dplyr::select(-"allelic_ratio", -"total_reads")
  table <- distinct(table)
#Assign name
  assign(name, table)
}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#

######------ Read in TMPs of all time points ------###### 
TPM_Embryonic <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "g TPM_Embryonic")
TPM_Young <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "i TPM_Young")
TPM_Adult <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "c TPM_Adult")
TPM_Aged <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "k TPM_Aged")

#extract median from all TPM tables
table_names <- c("TPM_Embryonic", "TPM_Young", "TPM_Adult", "TPM_Aged")
# Loop through the table names
for (name in table_names) {
  table <- get(name)
  # extract median data
  table$sample <- gsub("_1","",table$sample)
  table$sample <- gsub("_2","",table$sample)
  table$sample <- gsub("_3","",table$sample)
  table <- table %>% dplyr::select(-"TPM")
  table <- distinct(table)
  #Assign name
  assign(paste0(name,"_median"), table)
}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#

# combine all tables
allelic_ratios <- rbind(allelic_ratios_Embryonic, rbind(allelic_ratios_Young, rbind(allelic_ratios_Adult, allelic_ratios_Aged)))
TPM <- rbind(TPM_Embryonic, rbind(TPM_Young, rbind(TPM_Adult, TPM_Aged)))
TPM_median <- rbind(TPM_Embryonic_median, rbind(TPM_Young_median, rbind(TPM_Adult_median, TPM_Aged_median)))

# determine male escape genes
escapers_XY <- allelic_ratios%>%dplyr::filter(chr == "chrX" & sex == "male" & median_allelic_ratio<=0.90)
drop <- unique(escapers_XY$name)
drop <- c(drop , "Xist", "Tsix", "G530011O06Rik")

######---- Figure 3b: Boxplot of percentage of escape throughout development and aging ----###### 
Boxplot <- allelic_ratios %>% dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop)

Boxplot <- Boxplot %>% dplyr::mutate(escape_status = ifelse(median_allelic_ratio<=0.90, "escape","non-escape"))

#Count X-linked genes
Boxplot <- ddply(Boxplot,c("age","organ","sex","escape_status"),nrow)
colnames(Boxplot)[5] <- "Xlinked_genes"

#Add percentages for labels
count_table_sum <- Boxplot %>% dplyr::group_by(sex, organ, age) %>% dplyr::summarise(sum=sum(Xlinked_genes))
Boxplot <- merge(Boxplot, count_table_sum, by = c("sex", "organ", "age"))
Boxplot <- Boxplot %>%  dplyr::mutate(percentage= (Xlinked_genes/sum)*100)
Boxplot$percentage <- round(Boxplot$percentage, 3)

#Prepare dataframe for plotting
Boxplot <- Boxplot %>%  dplyr::filter(escape_status == "escape") %>% dplyr::select(-"sum")

#set colors and levels
Boxplot$age <- factor(Boxplot$age, levels = c("Embryonic", "Young","Adult", "Aged"))
age_colors <- c("#DCDCDD", "#E8D4AB", "#99AAB3", "#244357")
organ_colors <- c("#DCB465", "#8B1812","#244C51", "#C97A41", "#A77A76", "#555463", "#97A092")


#plot
Figure_2b <- ggplot(Boxplot, aes(x=age, y= percentage))+
  geom_boxplot(aes(fill=age), outliers= F) +
  geom_jitter(aes(color = organ),size=2, width = 0.1) +
  scale_fill_manual(values = age_colors) +
  scale_color_manual(values = organ_colors) +
  ylim(0,10) +
  theme(axis.text.x=element_text(size=9,angle = 90),
        axis.line = element_line(colour = 'black', size = 0.5),
        axis.ticks.y=element_line(colour = "black",size=0.5),
        axis.title.y=element_text(size=9),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="none",
        legend.title=element_text(size=9),
        legend.spacing = unit(0.5, 'cm')) +
  labs(x="",y="Percentage of escape") +
  ggtitle("Figure 2b")


######---- Figure 3c: Venn diagram of escape throughout development and aging ----###### 
Venn <- allelic_ratios %>% dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop & median_allelic_ratio<=0.90)

escapees_Embryonic <- unique(Venn$name[Venn$age == "Embryonic"])
escapees_Young <- unique(Venn$name[Venn$age == "Young"])
escapees_Adult <- unique(Venn$name[Venn$age == "Adult"])
escapees_Aged <- unique(Venn$name[Venn$age == "Aged"])

Venn_list <- list(escapees_Embryonic = escapees_Embryonic, escapees_Young = escapees_Young, escapees_Adult = escapees_Adult, escapees_Aged = escapees_Aged)
Figure_2c <- plot(euler(Venn_list), shape = "circle", quantities = list(type = c("percent", "counts")), input = "disjoint", 
                  legend = FALSE ,fills = paste0(age_colors, 60), edges = list(col=age_colors, lex = 5), main = "Figure 2c")


######---- Figure 3d: Barplot for number of escape genes per organ throughout development and aging ----###### 
Barplot <- allelic_ratios %>% dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop & median_allelic_ratio<=0.90)

#Count X-linked genes
Barplot <- ddply(Barplot,c("age","organ"),nrow)
colnames(Barplot)[3] <- "escape_genes"

#set levels
Barplot$age <- factor(Barplot$age, levels = c("Embryonic", "Young","Adult", "Aged"))
Barplot$organ <- factor(Barplot$organ,
                         levels = c('Ki','Lu', 'Mu','He','Br','Li','Sp'),ordered = TRUE)

#plot
Figure_2d <- ggplot(Barplot, aes( x=organ, y=escape_genes, fill=age)) + 
  geom_bar(color= "black", position="dodge", stat="identity") +
  scale_fill_manual(values = age_colors) +
  labs(x="",y="Number of escape genes") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust=0.1),
        axis.ticks.x=element_blank(),
        axis.line = element_line(colour = 'black', size = 0.2),
        axis.ticks.y=element_line(colour = "black",size=0.2),
        axis.title.y=element_text(size=9),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.spacing = unit(1.0, 'cm')) +
  ggtitle("Figure 2d")


######---- Figure 3e: Float chart of escape throughout development and aging ----###### 

#extract escape genes per organ and filter for them
Float <- allelic_ratios %>% dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop & median_allelic_ratio<=0.90)
Float$name_organ <- paste(Float$name, Float$organ, sep= "_")
escapees <- unique(Float$name_organ)

#
Float <- allelic_ratios %>% dplyr::mutate(name_organ= paste0(name,"_",organ)) %>% dplyr::filter(name_organ %in% escapees & sex == "female")
#add escape status
Float <- Float %>% dplyr::mutate(escape_status = ifelse(median_allelic_ratio<=0.90, "escape","non-escape"))
Float <- Float %>% select(name_organ, escape_status, age)

#prepare for alluvial plotting
Float <- pivot_wider(Float, names_from = age, values_from = escape_status)
Float[is.na(Float)] <- "NI" 
# NI (non-informative) means either not expressed (below selected expression cutoff) or not investigated (for Muscle and Spleen at the embryonic stage)

#to track escape status between each time points, dataframe was separated into 3
Float1 <- Float[,c("Embryonic", "Young")]
Float2 <- Float[,c("Young", "Adult")]
Float3 <- Float[,c("Adult", "Aged")]

# calculate frequencies of combinations in each dataframe
Float1 <- Float1 %>% 
  group_by(Float1[,1:2]) %>%
  mutate(freq = n()) %>%
  ungroup %>% 
  distinct(Float1[,1:2], .keep_all = TRUE) 

Float2 <- Float2 %>% 
  group_by(Float2[,1:2]) %>%
  mutate(freq = n()) %>%
  ungroup %>% 
  distinct(Float2[,1:2], .keep_all = TRUE) 

Float3 <- Float3 %>% 
  group_by(Float3[,1:2]) %>%
  mutate(freq = n()) %>%
  ungroup %>% 
  distinct(Float3[,1:2], .keep_all = TRUE) 

#plot
Float1_plot <- as.ggplot(expression(alluvial(Float1[,1:2], freq=Float1$freq,cex = 0.7)))
Float2_plot <- as.ggplot(expression(alluvial(Float2[,1:2], freq=Float2$freq, cex = 0.7)))
Float3_plot <- as.ggplot(expression(alluvial(Float3[,1:2], freq=Float3$freq, cex = 0.7)))
         
#combine to one plot
Figure_2e <- grid.arrange(Float1_plot, Float2_plot, Float3_plot, ncol=3)


######---- Figure 2f: Escape heatmap ----###### 
### select necessary columns
Heatmap <- allelic_ratios %>% dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop) %>% 
  dplyr::select(name, median_allelic_ratio, age, organ) 

Heatmap$age <- factor(Heatmap$age, levels = c("Embryonic", "Young", "Adult", "Aged")) 

#set colors and breaks
my_colors <- c("#1c642d", "#0F7031", "#357B30", "#4E8330", "#658C2D", "#78962A" ,"#8D9F25", "#A2A71D", "#B3B112", "#C97314" ,"#8b1913")
my_breaks <- seq(0.5,1,by = 0.05)

#LOOP

# organs with Embryonic samples
samples_E <- c("He", "Br", "Ki", "Li", "Lu")
for(sample in samples_E){
  ### get escape genes throughout all ages for each organ
  escape_list <- Heatmap %>% dplyr::filter(median_allelic_ratio<=0.90)
  escape_list <- escape_list[escape_list[,"organ"]==sample,]
  escape_list <- unique(escape_list$name)
  # separate organs
  organ <- Heatmap[Heatmap[,"organ"]==sample,]
  # select needed columns and genes in escaper list
  organ <- organ%>% dplyr::select(name, median_allelic_ratio, age) %>% dplyr::filter(name %in% escape_list)
  # prepare dataframe for heatmap plotting (make timepoints as column name with median_ratio as value)
  organ <- organ %>% pivot_wider(names_from = age, values_from = median_allelic_ratio)
  #reorder by AR
  organ <- organ[order(-organ$`Embryonic`),]
  #prepare for plotting
  organ <- organ %>% remove_rownames %>% column_to_rownames(var="name")
  organ <- as.matrix(organ)
  #plot
  plot <- as.ggplot(function() {heatmap.2(organ, Rowv = NA, Colv = NA, tracecol=NULL, sepcolor="black", keysize=1, dendrogram = "none",
                                             cexCol=1,cexRow = 0.5,na.color = "white", breaks = my_breaks,
                                             sepwidth=c(0.005,0.01), colsep=1:ncol(organ), rowsep=1:nrow(organ),
                                             col=colorRampPalette(my_colors), main = sample)})
  #write output
  assign(paste0("Heatmap_",sample), plot)
}

# organs without Embryonic samples
samples <- c("Sp", "Mu")
for(sample in samples){
  ### get escape genes throughout all ages for each organ
  escape_list <- Heatmap %>% dplyr::filter(median_allelic_ratio<=0.90)
  escape_list <- escape_list[escape_list[,"organ"]==sample,]
  escape_list <- unique(escape_list$name)
  # separate organs
  organ <- Heatmap[Heatmap[,"organ"]==sample,]
  # select needed columns and genes in escaper list
  organ <- organ%>% dplyr::select(name, median_allelic_ratio, age) %>% dplyr::filter(name %in% escape_list)
  # prepare dataframe for heatmap plotting (make timepoints as column name with median_ratio as value)
  organ <- organ %>% pivot_wider(names_from = age, values_from = median_allelic_ratio)
  #reorder by AR
  organ <- organ[order(-organ$`Young`),]
  #prepare for plotting
  organ <- organ %>% remove_rownames %>% column_to_rownames(var="name")
  organ <- as.matrix(organ)
  #plot
  plot <- as.ggplot(function() {heatmap.2(organ, Rowv = NA, Colv = NA, tracecol=NULL, sepcolor="black", keysize=1, dendrogram = "none",
                                          cexCol=1,cexRow = 0.5,na.color = "white", breaks = my_breaks,
                                          sepwidth=c(0.005,0.01), colsep=1:ncol(organ), rowsep=1:nrow(organ),
                                          col=colorRampPalette(my_colors), main = sample)})
  #write output
  assign(paste0("Heatmap_",sample), plot)
}

Figure_2f <- grid.arrange(Heatmap_Ki, Heatmap_He, Heatmap_Lu, Heatmap_Mu, Heatmap_Br, Heatmap_Sp, Heatmap_Li, 
             layout_matrix = rbind(c(1, 2 , 3, 4),
                                   c(NA, 5, 6, 7)))


######---- Figure 2g and h: Allelic ratio drops ----###### 

#Defining early and age-specific escape genes per organ
Early_timepoints <- c("Embryonic", "Young", "Adult")
Early_escapees <- allelic_ratios %>% dplyr::mutate(name_organ= paste0(allelic_ratios$name,"_", allelic_ratios$organ)) %>% 
  dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop & median_allelic_ratio<=0.90 & age %in% Early_timepoints)
Early_escapees_name <- unique(Early_escapees$name_organ)

Age_spec_escapees <- allelic_ratios %>% dplyr::mutate(name_organ= paste0(allelic_ratios$name,"_", allelic_ratios$organ)) %>% 
  dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop & median_allelic_ratio<=0.90 & age == "Aged" & ! name_organ %in% Early_escapees_name)
Age_spec_escapees_name <- unique(Age_spec_escapees$name_organ)
#Escape genes that occur in the Aged sample but not in the other time points in the respective organ are considered age-specific escape genes
#Example: Smpx is an age-specific escape gene but only in Muscle not in Heart (where it already escape in the Adult sample)

#Create data frame for plotting by filtering for early and age-specific escapees
AR_drop_plot <- allelic_ratios %>% dplyr::mutate(name_organ= paste0(allelic_ratios$name,"_", allelic_ratios$organ)) %>% 
  dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop & name_organ %in% c(Early_escapees_name, Age_spec_escapees_name))
AR_drop_plot <- AR_drop_plot %>% dplyr::mutate(escape_status = ifelse(name_organ %in% Age_spec_escapees_name, "age-spec_escapee",
                                                        ifelse(name_organ %in% Early_escapees_name, "escapee","")))
AR_drop_plot <- AR_drop_plot %>% dplyr::select(median_allelic_ratio, organ, age, escape_status, name_organ)

#Data frame is now filtered for genes that: 1. are informative in Young, Adult and Aged and 2. show a ∆AR >= 0.1, and 3. show the lowest AR in the aged time point
Delta_filter <- pivot_wider(AR_drop_plot, names_from = age, values_from = median_allelic_ratio)

#1. Filter for genes that are informative in Young, Adult and Aged
Delta_filter <- Delta_filter[!is.na(Delta_filter$Aged), ]
Delta_filter <- Delta_filter[!is.na(Delta_filter$Adult), ]
Delta_filter <- Delta_filter[!is.na(Delta_filter$Young), ]

#2. Calculate the Delta_filter AR to filter for genes that show ∆AR >= 0.1
calculate_delta <- function(row) {
  max(row, na.rm = TRUE) - min(row, na.rm = TRUE)
}
Delta_filter$deltaAR <- apply(Delta_filter[, c("Embryonic", "Young", "Adult", "Aged")], 1, calculate_delta)
Delta_filter <- Delta_filter %>% dplyr::filter(deltaAR >= 0.1)

#3. Filter for genes that show the lowest AR in the aged time point
Delta_filter <- Delta_filter %>% dplyr::mutate(status = ifelse(Aged != pmin(Embryonic, Young, Adult, Aged, na.rm = TRUE), "drop", ""))
filter_out <- Delta_filter %>% dplyr::filter(status == "drop")
Delta_filter <- Delta_filter %>% dplyr::filter(!name_organ %in% filter_out$name_organ)
#now the Delta_filter data frame can be used to filter original AR_drop_plot data frame for column name_organ


# Filter original AR_drop_plot table for genes present in the delta table
AR_drop_plot <- AR_drop_plot %>%  dplyr::filter(name_organ %in% Delta_filter$name_organ)
# Add ∆AR column for plotting transparency
combine <- Delta_filter %>%  dplyr::select(name_organ, deltaAR)
AR_drop_plot <- merge(AR_drop_plot, combine, by="name_organ",  all=T)

# prepare for plotting
AR_drop_plot$age <- factor(AR_drop_plot$age, levels = c("Embryonic", "Young", "Adult", "Aged")) 
label <- AR_drop_plot %>% dplyr::filter(name_organ %in% c("Tlr8_Lu","Kdm6a_Mu", "Plp1_He") & age == "Aged")

# Plot dotplot with connective lines
#Figure g: age_specific escapees
Figure_2g <- ggplot(AR_drop_plot[AR_drop_plot$escape_status == "age-spec_escapee", ], aes(x = age, y = median_allelic_ratio, group = name_organ, color = escape_status, shape = organ, alpha = deltaAR)) +
  geom_point(size =3) +
  geom_line() +
  scale_color_manual(values = c("#244357")) +  
  scale_alpha_continuous(range = c(0.4, 1)) +  
  scale_shape_manual(values = c("Br" = 0, "Li" = 1, "Lu" = 2, "He" = 7, "Sp" = 6, "Ki" = 3, "Mu" = 5)) +  
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid.major.x = element_line(color = "grey", size = 0.5),
        panel.grid.major.y = element_blank(), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),) +
  labs(title = "Figure 2g",
       x = "",
       y = "Median allelic ratio",
       color = "") 

#Figure h: early escapees
Figure_2h <- ggplot(AR_drop_plot[AR_drop_plot$escape_status == "escapee", ], aes(x = age, y = median_allelic_ratio, group = name_organ, shape = organ,alpha = deltaAR)) +
  geom_point(size =3) +
  geom_line() +
  geom_text(data= label, aes(label = name_organ), hjust = -0.2, vjust = -0.0) +
  scale_alpha_continuous(range = c(0.4, 1)) + 
  scale_shape_manual(values = c("Br" = 0, "Li" = 1, "Lu" = 2, "He" = 7, "Sp" = 6, "Ki" = 3, "Mu" = 5)) +  
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid.major.x = element_line(color = "grey", size = 0.5),
        panel.grid.major.y = element_blank(), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  labs(title = "Figure 2h",
       x = "",
       y = "Median allelic ratio") 

# Function to extract legend
get_legend <- function(Figure_2g) {
  tmp <- ggplotGrob(Figure_2g)
  leg <- gtable::gtable_filter(tmp, "guide-box")
  return(leg)
}

# Extract the legend from one plot
shared_legend <- get_legend(Figure_2g)

# Remove legends from the plots
Figure_2g <- Figure_2g + theme(legend.position = "none")

# Arrange plots and shared legend
Figure_2gh <- grid.arrange(
  Figure_2g, Figure_2h, shared_legend, ncol = 2, 
  layout_matrix = rbind(c(1, 3), c(2, 3)),
  widths = c(4, 1))

######---- Figure 3i: Boxplot TPM fold change ----###### 
# Calculate fold change (pseudocount)
Foldchange <- TPM_median %>% dplyr::select(-"sample")

Foldchange$medianTPM <- Foldchange$medianTPM + 1
Foldchange <- Foldchange%>%pivot_wider(names_from = sex, values_from = medianTPM)
Foldchange$FC <-(Foldchange$female) / (Foldchange$male)
Foldchange=Foldchange[complete.cases(Foldchange), ]
Foldchange <- Foldchange %>%  dplyr::filter(male != 1, female != 1) 

# exclude Xist
TPM_Boxplot <- Foldchange %>% dplyr::filter(name != "Xist")

#Add escape status
TPM_Boxplot <- TPM_Boxplot %>% dplyr::mutate(name_organ= paste0(TPM_Boxplot$name,"_", TPM_Boxplot$organ)) %>%
  mutate(escape_status = ifelse(name_organ %in% Age_spec_escapees_name, "age-spec_escapee",
                             ifelse(name_organ %in% Early_escapees_name, "escapee", NA)))
TPM_Boxplot$escape_status[is.na(TPM_Boxplot$escape_status)] <- "non-escape"

#Prepare for plotting
TPM_Boxplot <- TPM_Boxplot %>% dplyr::mutate(age_escape_status= paste0(TPM_Boxplot$age,"_", TPM_Boxplot$escape_status))
TPM_Boxplot$age_escape_status <- factor(TPM_Boxplot$age_escape_status, levels = c("Young_non-escape", "Young_escapee", "Young_age-spec_escapee",
                                                                          "Adult_non-escape", "Adult_escapee", "Adult_age-spec_escapee",
                                                                          "Aged_non-escape", "Aged_escapee", "Aged_age-spec_escapee"))
#set colors
colors <- c("#244357", "#1C642D","#DADADA")
colors <- rep(colors, times = 3)

#plot
Figure_2i <- ggplot(TPM_Boxplot, aes(x=age_escape_status, y=FC, fill=escape_status)) +
  geom_boxplot(size = 0.5,outlier.size = 0.5,lwd=0.5) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.ticks.x=element_line(colour = "darkgrey",size=0.2),
    axis.title.x=element_text(size=12),
    axis.text.x = element_text(angle = 90, hjust=0.1),
    axis.line = element_line(colour = 'darkgrey', size = 0.2),
    axis.ticks.y=element_line(colour = "darkgrey",size=0.2),
    axis.title.y=element_text(size=12),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.border = element_rect(colour = "darkgrey", fill=NA, size=0.5),
    legend.position="right") +
  scale_fill_manual(values = colors) +
  ylim(0,2) +
  labs(title = "Figure 2i",
       x = "",
       y = "TPM fold change (female vs. male)",
       fill="")


######---- Figure 2j: TPM dot plot ----###### 
selected_genes <- c("Reps2_Ki", "Gpr34_Br", "Tlr8_Lu", "Sh3kbp1_He", "Tsc22d3_Mu")
TPM_plot <- TPM %>% dplyr::mutate(name_organ= paste0(name,"_", organ), age_sex= paste0(age,"_", sex)) %>% 
                      dplyr::filter(name_organ %in% selected_genes & age != "Embryonic")

#prepare for plotting
TPM_plot$age_sex <- factor(TPM_plot$age_sex, levels = c("Young_male", "Young_female", 
                                                        "Adult_male", "Adult_female", 
                                                        "Aged_male", "Aged_female"))
#set colors
sex_colors <- c("#901C14","#2B3186")

#plot
Figure_2j <- ggplot(TPM_plot, aes(x=age_sex, y=TPM, color =sex)) + 
  geom_point()+
  stat_summary(fun.y=mean, shape=95, size=2)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black")+
  facet_wrap(~name_organ, scales = "free", nrow = 1, ncol = 5)+
  scale_color_manual(values = sex_colors) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid.major.x = element_line(color = "grey", size = 0.5),
        panel.grid.major.y = element_blank(), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust=0.1),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  labs(title = "Figure 2j",
       x = "",
       y = "TPM") 

######---- Figure 2k: Disease associations of escape genes ----###### 
Disease_Link <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "l Disease_association", skip = 1)
Disease_Link <- Disease_Link %>% dplyr:: filter(escape == TRUE)

non_conserved_escape = c("5530601H04Rik", "LOC108167322", "F630028O10Rik", "Gm15232")

Disease_Link <- Disease_Link %>% dplyr::mutate(conserved = ifelse(name %in% non_conserved_escape, FALSE, TRUE))
Disease_Link <- Disease_Link %>% dplyr::mutate(category= ifelse(conserved == TRUE & disease_association == TRUE, "disease_association",
                                                                ifelse(conserved == FALSE ,"not_conserved", "no_disease_association")))

#Count frequencies
Disease_Link <- distinct(Disease_Link %>% dplyr::select(-disease_id, -disease_term))
count_table_disease_association <- ddply(Disease_Link,"category",nrow)
count_table_disease_association <- count_table_disease_association %>% mutate(percentage = round(100 * V1 / sum(V1), 1))

# calculate position of labels
count_table_disease_association <- count_table_disease_association[c(3,2,1),]
count_table_disease_association <- count_table_disease_association %>%
  mutate(cumulative = cumsum(V1),  # Cumulative sum of V1
         midpoint = cumulative - V1 / 2) 

count_table_disease_association$label <- paste0(count_table_disease_association$percentage, "%, n=", count_table_disease_association$V1)

#set colors
disease_colors <- c("#A62C24","#BFD7E1" ,"#DADADA")

#plot
Figure_2k <- ggplot(count_table_disease_association, aes(x = "", y = V1, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(theta = "y") +
  scale_fill_manual(values=disease_colors) +
  scale_size_identity()+
  theme_void() + 
  theme(legend.position = "bottom") +
  geom_text(aes(label = label, y = midpoint)) +
  labs(title = "Figure 2k",
       subtitle = "Escapees: 52 genes",
       fill="") 

######---- Figure 2l: Disease enrichment escape genes ----###### 
Disease_Enrichment <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "l Disease_association", skip = 1)
Disease_Enrichment <- distinct(Disease_Enrichment %>% dplyr::select(-disease_id, -disease_term))
count_table_disease_enrichment <- as.data.table(table(Disease_Enrichment$escape, Disease_Enrichment$disease_association))
colnames(count_table_disease_enrichment) <- c("escape", "disease_association", "count")
count_table_disease_enrichment <- count_table_disease_enrichment %>% dplyr::group_by(escape) %>% dplyr::mutate(percentage = round(100 * count / sum(count), 1))

# prepare for plotting
count_table_disease_enrichment$escape <- ifelse(count_table_disease_enrichment$escape == TRUE, "escape", "non-escape")
count_table_disease_enrichment$disease_association <- ifelse(count_table_disease_enrichment$disease_association == TRUE, "disease_association", "no_disease_association")
count_table_disease_enrichment$escape <- factor(count_table_disease_enrichment$escape, levels = c("non-escape", "escape"))
count_table_disease_enrichment$disease_association <- factor(count_table_disease_enrichment$disease_association, levels = c("no_disease_association","disease_association"))

# set colors
enrichment_colors <- c("#BFD7E1", "#A62C24")

#plot
Figure_2l <- ggplot(count_table_disease_enrichment, aes(x=escape, y=percentage, fill=disease_association)) +
  geom_bar(stat = "identity") +
  labs(title = "Figure 2l",
       x = "",
       y = "Disease association [%]") +
  theme_minimal()+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10))+
  scale_fill_manual(values= enrichment_colors)

# Fisher's exact test for significance
values <- count_table_disease_enrichment$count
fishers_matrix <- matrix(values, nrow = 2, byrow = TRUE)
fisher.test(fishers_matrix)

Figure_2b
Figure_2c
Figure_2d
Figure_2e <- grid.arrange(Float1_plot, Float2_plot, Float3_plot, ncol=3)
Figure_2f <- grid.arrange(Heatmap_Ki, Heatmap_He, Heatmap_Lu, Heatmap_Mu, Heatmap_Br, Heatmap_Sp, Heatmap_Li, 
                          layout_matrix = rbind(c(1, 2 , 3, 4),
                                                c(NA, 5, 6, 7)))
Figure_2gh <- grid.arrange(
  Figure_2g, Figure_2h, shared_legend, ncol = 2, 
  layout_matrix = rbind(c(1, 3), c(2, 3)),
  widths = c(4, 1))
Figure_2i
Figure_2j
Figure_2k
Figure_2l

