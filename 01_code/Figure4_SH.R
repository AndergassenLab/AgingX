##***************************************##
##  Figure 4 R script to generate plots  ##
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
library(karyoploteR)
library(gdata)
library(ggpubr)
library(ggrepel)
library(readxl)
library(gridExtra)
library(ggplotify)

######------ Set environment  ------###### 
input <- "/Users/shoelzl/Desktop/mnt_new2/ge43sik2/WP1_aging_SH/manuscript/EscaperPaper/20241216_nature_aging_revision/revised_supplemetary_tables/"
output <- "/Users/shoelzl/Desktop/mnt_new2/ge43sik2/WP1_aging_SH/10_script_submission/Figure_04/"
  
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

# combine all tables
allelic_ratios <- rbind(allelic_ratios_Embryonic, rbind(allelic_ratios_Young, rbind(allelic_ratios_Adult, allelic_ratios_Aged)))

# determine male escape genes
escapers_XY <- allelic_ratios%>%dplyr::filter(chr == "chrX" & sex == "male" & median_allelic_ratio<=0.90)
drop <- unique(escapers_XY$name)
drop <- c(drop , "Xist", "Tsix", "G530011O06Rik")

######---- Figure 4a: Karyoplot across all ages ----###### 

Karyoplot <- allelic_ratios %>% dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop & median_allelic_ratio<=0.90)

# get coordinates for all escapees
label <- Karyoplot %>%  dplyr::select("chr", "start", "end", "name")
label <- distinct(label)
label$chr <- gsub("chr","",label$chr)
genes_coord_all <- regioneR::toGRanges(as.data.frame(label)) 
seqlevelsStyle(genes_coord_all) <- "UCSC"

# separate organs and ages and save names of resulting data frames in list
ages <- c("Embryonic", "Young", "Adult", "Aged")
organs <- c("He", "Br", "Li", "Lu", "Sp", "Mu", "Ki")
plotting_vector <- list()


for(a in ages){
  for(o in organs){
    df <- Karyoplot %>% dplyr::filter(organ == o & age== a)
    #Assign name and add to vector
    df_name <- paste0(a, "_", o)
    plotting_vector[[df_name]] <- df
  }}

#remove Embryonic Muscle and Spleen, because these organs were not existing at that time point
plotting_vector[["Embryonic_Mu"]] <- NULL
plotting_vector[["Embryonic_Sp"]] <- NULL

# Create lists for each age
names_plotting_vector <- names(plotting_vector)
plotting_vector_Embryonic <- plotting_vector[grep("Embryonic", names_plotting_vector, value = TRUE)]
plotting_vector_Young <- plotting_vector[grep("Young", names_plotting_vector, value = TRUE)]
plotting_vector_Adult <- plotting_vector[grep("Adult", names_plotting_vector, value = TRUE)]
plotting_vector_Aged <- plotting_vector[grep("Aged", names_plotting_vector, value = TRUE)]

# set position for plot to separate ages
position_vector_Embryonic <- c(1:5)
position_vector_Young <- c(7:13)
position_vector_Adult <- c(15:21)
position_vector_Aged <- c(23:29)

# plot
# chromosome ideogram
Figure_4a <- as.ggplot(expression(
  kp <- plotKaryotype(genome="mm10",chromosomes="chrX",plot.type = 2),
  # add title
  kpAddMainTitle(kp, main= "Figure 4a",col="black"),
  # add scale
  kpAddBaseNumbers(kp),
  # add gene names
  kpPlotMarkers(kp, data = genes_coord_all,
                labels = genes_coord_all$name,
                text.orientation = "vertical",
                r1 = 0.5, cex = 0.9,
                ymax=5),
  # add "kpRect" for escape genes for all organs by looping through
  for(i in seq_along(plotting_vector_Embryonic)){
    plot_value <- plotting_vector_Embryonic[[i]]
    position_value <- position_vector_Embryonic[i]
    df_name <- names(plotting_vector_Embryonic)[i]
    kpRect(kp, chr= "chrX", x0=plot_value$start, x1=plot_value$end, y0=0, y1=1, 
           ymin = 0,ymax = 1, r0=autotrack(position_value,26),data.panel = 2,
           border= ifelse(grepl('He', df_name),"#8B1812",
                        ifelse(grepl('Br', df_name),"#DCB465",
                             ifelse(grepl('Li', df_name),"#C97A41",
                                 ifelse(grepl('Lu', df_name),"#A77A76",
                                      ifelse(grepl('Mu', df_name),"#555463",
                                          ifelse(grepl('Ki', df_name),"#244C51",
                                              ifelse(grepl('Sp', df_name),"#97A092",""))))))))
    kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(position_value,26),data.panel = 2)
  },
  for(i in seq_along(plotting_vector_Young)){
    plot_value <- plotting_vector_Young[[i]]
    position_value <- position_vector_Young[i]
    df_name <- names(plotting_vector_Young)[i]
    kpRect(kp, chr= "chrX", x0=plot_value$start, x1=plot_value$end, y0=0, y1=1, 
           ymin = 0,ymax = 1, r0=autotrack(position_value,26),data.panel = 2,
           border= ifelse(grepl('He', df_name),"#8B1812",
                          ifelse(grepl('Br', df_name),"#DCB465",
                                 ifelse(grepl('Li', df_name),"#C97A41",
                                        ifelse(grepl('Lu', df_name),"#A77A76",
                                               ifelse(grepl('Mu', df_name),"#555463",
                                                      ifelse(grepl('Ki', df_name),"#244C51",
                                                             ifelse(grepl('Sp', df_name),"#97A092",""))))))))
    kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(position_value,26),data.panel = 2)
  },
  for(i in seq_along(plotting_vector_Adult)){
    plot_value <- plotting_vector_Adult[[i]]
    position_value <- position_vector_Adult[i]
    df_name <- names(plotting_vector_Adult)[i]
    kpRect(kp, chr= "chrX", x0=plot_value$start, x1=plot_value$end, y0=0, y1=1, 
           ymin = 0,ymax = 1, r0=autotrack(position_value,26),data.panel = 2,
           border= ifelse(grepl('He', df_name),"#8B1812",
                          ifelse(grepl('Br', df_name),"#DCB465",
                                 ifelse(grepl('Li', df_name),"#C97A41",
                                        ifelse(grepl('Lu', df_name),"#A77A76",
                                               ifelse(grepl('Mu', df_name),"#555463",
                                                      ifelse(grepl('Ki', df_name),"#244C51",
                                                             ifelse(grepl('Sp', df_name),"#97A092",""))))))))
    kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(position_value,26),data.panel = 2)
  },
  for(i in seq_along(plotting_vector_Aged)){
    plot_value <- plotting_vector_Aged[[i]]
    position_value <- position_vector_Aged[i]
    df_name <- names(plotting_vector_Aged)[i]
    kpRect(kp, chr= "chrX", x0=plot_value$start, x1=plot_value$end, y0=0, y1=1, 
           ymin = 0,ymax = 1, r0=autotrack(position_value,26),data.panel = 2,
           border= ifelse(grepl('He', df_name),"#8B1812",
                          ifelse(grepl('Br', df_name),"#DCB465",
                                 ifelse(grepl('Li', df_name),"#C97A41",
                                        ifelse(grepl('Lu', df_name),"#A77A76",
                                               ifelse(grepl('Mu', df_name),"#555463",
                                                      ifelse(grepl('Ki', df_name),"#244C51",
                                                             ifelse(grepl('Sp', df_name),"#97A092",""))))))))
    kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(position_value,26),data.panel = 2)
  }
))


######---- Figure 4b: Enrichment plot ----###### 

#To compute whether there is an enrichment of aged escapees at distal chromosome regions, we count the number of Aged and early (Embryonic, Young and Adult combined) along the X chr.
#To count how many escape genes occur in each 20Mb sliding window across the X chromosome, an annotation was created (10Mb overlap).

# Create annotation for 20 Mb sliding windows on X as .bed file 
options(scipen = 999)
n <- 18
start = seq(0, by = 10000000, length.out = n)
end = c(seq(20000000, by = 10000000, length.out = n - 2), 171031299, 171031299)
bed_20MbSW_annotation <- data.frame(
  chr = rep("chrX", n),
  start = seq(0, by = 10000000, length.out = n),
  end = c(seq(20000000, by = 10000000, length.out = n - 2), 171031299, 171031299),
  description = paste0("chrX:", start, "-", end),
  value = rep(0, n),
  dot = rep(".", n),
  stringsAsFactors = FALSE
)

# Create .bed file for escape genes
#Early timpoints
early_timepoints <- c("Embryonic", "Young", "Adult")
esc_earlyTimepoints <- allelic_ratios %>% dplyr::filter(sex == "female" & 
                                                           chr == "chrX" & 
                                                           ! name %in% drop & 
                                                           median_allelic_ratio<=0.90 & 
                                                           age %in% early_timepoints)
esc_earlyTimepoints <- esc_earlyTimepoints %>% dplyr::select(chr, start, end, name)
esc_earlyTimepoints <- distinct(esc_earlyTimepoints)
esc_earlyTimepoints$V5 <- 1000
esc_earlyTimepoints$V6 <- "."

#Aged
esc_Aged <- allelic_ratios %>% dplyr::filter(sex == "female" & 
                                                          chr == "chrX" & 
                                                          ! name %in% drop & 
                                                          median_allelic_ratio<=0.90 & 
                                                          age == "Aged")
esc_Aged <- esc_Aged %>% dplyr::select(chr, start, end, name)
esc_Aged <- distinct(esc_Aged)
esc_Aged$V5 <- 1000
esc_Aged$V6 <- "."

#Save annotation and escape .bed files
write.table(bed_20MbSW_annotation,paste0(output,"20Mb_sliding_window_X.bed"), quote=FALSE, sep = "\t", col.names = F, row.names =F)
write.table(esc_earlyTimepoints,paste0(output,"EscapeGenesAcrossAllOrgans_earlyTimepoints.bed"), quote=FALSE, sep = "\t", col.names = F, row.names =F)
write.table(esc_Aged,paste0(output,"EscapeGenesAcrossAllOrgans_Aged.bed"), quote=FALSE, sep = "\t", col.names = F, row.names =F)

# Escape genes (saved in a .bed format) were counted within each window using intersectBed
# requirement: bedtools
# Define file paths
annotation <- paste0(output,"20Mb_sliding_window_X.bed")
Early <- paste0(output,"EscapeGenesAcrossAllOrgans_earlyTimepoints.bed")
Aged <- paste0(output,"EscapeGenesAcrossAllOrgans_Aged.bed")

# Run intersectBed command to count escape genes within sliding windows
system(paste("intersectBed -a", annotation, "-b", Early, "-wa -c >", paste0(output,"EscapeCount_EarlyTimepoints.bed")))
system(paste("intersectBed -a", annotation, "-b", Aged, "-wa -c >", paste0(output,"EscapeCount_Aged.bed")))

# Read the results back into R
Early_Count <- read.table(paste0(output,"EscapeCount_EarlyTimepoints.bed"), header = FALSE, sep = "\t")
Aged_Count <- read.table(paste0(output,"EscapeCount_Aged.bed"), header = FALSE, sep = "\t")

# Combine counts
Early_Count <- Early_Count %>% dplyr::select(-"V4", -"V5", -"V6")
Aged_Count <- Aged_Count %>% dplyr::select(-"V4", -"V5", -"V6")

Early_Count <- setNames(Early_Count, c("chr", "start", "end", "escape_count_earlyTP"))
Aged_Count <- setNames(Aged_Count, c("chr", "start", "end", "escape_count_Aged"))

SlidingWindow_counts <- merge(Early_Count, Aged_Count, by = c("chr", "start", "end"))

SlidingWindow_counts <- SlidingWindow_counts %>% arrange(end)

#Compute enrichment
SlidingWindow_counts$enrich <- ((SlidingWindow_counts$escape_count_Aged) / (SlidingWindow_counts$escape_count_Aged + SlidingWindow_counts$escape_count_earlyTP)) *100

#Binomial test to test significance
for (i in 1:nrow(SlidingWindow_counts)) {
  row <- i
  SlidingWindow_counts$binom_test[row] <- pbinom(min(SlidingWindow_counts[row,4:5]), SlidingWindow_counts[row,4]+SlidingWindow_counts[row,5], 0.5, lower.tail=TRUE, log = FALSE)
}

#Add delta of escape genes for color gradient (= escape gain in aging)
SlidingWindow_counts <- SlidingWindow_counts %>% dplyr::mutate(sum_escape_count = escape_count_earlyTP + escape_count_Aged)
SlidingWindow_counts <- SlidingWindow_counts %>% dplyr::mutate(delta = escape_count_Aged - escape_count_earlyTP)

#color for delta
# Create a color palette with 5 colors from light grey to green
color_palette <- colorRampPalette(c("gray93", "darkgreen"))(5)

# Normalize the delta column to range between 1 and 5
normalized_values <- cut(SlidingWindow_counts$delta, breaks=5, labels=FALSE)

# Assign colors based on normalized values
SlidingWindow_counts$delta_colors <- color_palette[normalized_values]

#set scales
scale_top <- 100
scale_bottom <- 40

#plot
Figure_4b <- as.ggplot(expression(
kp <- plotKaryotype(genome="mm10",chromosomes="chrX",plot.type = 4),
kpAddBaseNumbers(kp),
kpAxis(kp,ymin = scale_bottom ,ymax = scale_top, tick.pos=c(40,50,60,70,80,90,100)),
kpRect(kp, chr="chrX", x0=SlidingWindow_counts$start, x1=SlidingWindow_counts$end, y0=SlidingWindow_counts$enrich, y1=SlidingWindow_counts$enrich, 
       ymin = scale_bottom,ymax = scale_top, border= SlidingWindow_counts$delta_colors),
# Add a legend
legend_labels <- c(paste0("≤ ", round(max(SlidingWindow_counts$delta) * 1/5, 2)),
                   paste0("≤ ", round(max(SlidingWindow_counts$delta) * 2/5, 2)),
                   paste0("≤ ", round(max(SlidingWindow_counts$delta) * 3/5, 2)),
                   paste0("≤ ", round(max(SlidingWindow_counts$delta) * 4/5, 2)),
                   paste0("≤ ", max(SlidingWindow_counts$delta))),
legend("bottomright", legend=legend_labels, fill=color_palette, title="No. age-specific escape gain", cex=0.5)
))

######---- Plot Figure 4 ----###### 
grid.arrange(Figure_4a, Figure_4b,
             layout_matrix = rbind(c(1),
                                   c(1),
                                   c(1),
                                   c(2)))


