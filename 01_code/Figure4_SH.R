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
library(alluvial)
library(gridExtra)
library(ggplotify)
library(scales)

######------ Set environment  ------###### 
input <- "/Users/shoelzl/Desktop/mnt_new2/ge43sik2/WP1_aging_SH/manuscript/EscaperPaper/20241216_nature_aging_revision1/revised_supplemetary_tables/"
output <- "/Users/shoelzl/Desktop/mnt_new2/ge43sik2/WP1_aging_SH/10_script_submission/Figure_04/"
  
######------ Read in Allelic Ratios of all time points ------###### 
allelic_ratios_Embryonic <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "f AllelicRatio_Embryonic")
allelic_ratios_Young <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "h AllelicRatio_Young")
allelic_ratios_Adult <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "b AllelicRatio_Adult")
allelic_ratios_Aged <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "j AllelicRatio_Aged")

######------ Read in ATAC tables ------###### 
ATAC_count <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "p ATAC peaks (SNPsplit)")
ATAC_allelic_ratios <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "q ATAC peaks (AllelicRatio)")

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

ATAC_allelic_ratios <- ATAC_allelic_ratios %>% dplyr::filter(min_total_reads >= min)
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

######---- Plot Figure 4 part 1 ----###### 
grid.arrange(Figure_4a, Figure_4b,
             layout_matrix = rbind(c(1),
                                   c(1),
                                   c(1),
                                   c(2)))


######---- Figure 4d: Fractions of ATAC peaks ----######
Barplot <- ATAC_count[,c(1,2,3,4,7,10,13,16,19,22,25,28)]

Barplot <- pivot_longer(Barplot, cols = starts_with("mean"), names_to = "sample", values_to = "peak_count")

#add chromosome classification (Xa, Xi or autosome), age and organ
Barplot <- Barplot %>% dplyr::filter(chr != "chrY")
Barplot$chr_class <- ifelse(!grepl('chrX', Barplot$chr),"autosome",
                            ifelse(grepl('CAST', Barplot$sample) & grepl('chrX', Barplot$chr),"Xi",
                                   ifelse(grepl('BL6', Barplot$sample) & grepl('chrX', Barplot$chr),"Xa", "")))

Barplot$age <- ifelse(grepl('9w', Barplot$sample),"Adult",
                      ifelse(grepl('78w', Barplot$sample),"Aged", ""))

Barplot$organ <- ifelse(grepl('Li', Barplot$sample),"Li",
                        ifelse(grepl('Ki', Barplot$sample),"Ki", ""))

#add percentages
count_table <- Barplot %>% dplyr::select(age, chr_class, organ, peak_count) %>% 
  dplyr::group_by(age, chr_class, organ) %>% 
  dplyr::summarise(peak_count=sum(peak_count)) %>% ungroup()
summary <- count_table %>% dplyr::group_by(organ, chr_class) %>% 
  dplyr::summarise(sum=sum(peak_count)) %>% ungroup()
count_table <- left_join(count_table, summary, by = c("chr_class", "organ"))
count_table$perc <- round((count_table$peak_count/count_table$sum)*100,1)

# set order and colors
count_table$organ <- factor(count_table$organ,
                            levels = c("Li", "Ki"),ordered = TRUE)
age_colors <- c("#9AAAB4", "#234257")

#plot
Figure_4d <- ggplot(count_table, aes(x = chr_class, y = perc, fill = age)) +
  geom_bar(stat = "identity", col= "black") +
  facet_wrap(~organ) +
  scale_fill_manual(values= age_colors) +
  labs(title = "Figure 4d",
       x = "",
       y = "Fraction of ATAC peaks [%]",
       fill = "") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Fisher's exact test for significance testing 
#Liver Xa
Count_A = count_table[count_table$organ == "Li" & count_table$age == "Aged" & count_table$chr_class == "Xa", "peak_count"][[1]]
Total_A = count_table[count_table$organ == "Li" & count_table$age == "Aged" & count_table$chr_class == "autosome", "peak_count"][[1]]+Count_A
Count_B = count_table[count_table$organ == "Li" & count_table$age == "Adult" & count_table$chr_class == "Xa", "peak_count"][[1]]
Total_B = count_table[count_table$organ == "Li" & count_table$age == "Adult" & count_table$chr_class == "autosome", "peak_count"][[1]]+Count_B
contingency_table_Li_Xa <- matrix(c(Count_A, Total_A - Count_A, Count_B, Total_B - Count_B), nrow = 2)
fisher.test(contingency_table_Li_Xa)

#Liver Xi
Count_A = count_table[count_table$organ == "Li" & count_table$age == "Aged" & count_table$chr_class == "Xi", "peak_count"][[1]]
Total_A = count_table[count_table$organ == "Li" & count_table$age == "Aged" & count_table$chr_class == "autosome", "peak_count"][[1]]+Count_A
Count_B = count_table[count_table$organ == "Li" & count_table$age == "Adult" & count_table$chr_class == "Xi", "peak_count"][[1]]
Total_B = count_table[count_table$organ == "Li" & count_table$age == "Adult" & count_table$chr_class == "autosome", "peak_count"][[1]]+Count_B
contingency_table_Li_Xi <- matrix(c(Count_A, Total_A - Count_A, Count_B, Total_B - Count_B), nrow = 2)
fisher.test(contingency_table_Li_Xi)

#Kidney Xa
Count_A = count_table[count_table$organ == "Ki" & count_table$age == "Aged" & count_table$chr_class == "Xa", "peak_count"][[1]]
Total_A = count_table[count_table$organ == "Ki" & count_table$age == "Aged" & count_table$chr_class == "autosome", "peak_count"][[1]]+Count_A
Count_B = count_table[count_table$organ == "Ki" & count_table$age == "Adult" & count_table$chr_class == "Xa", "peak_count"][[1]]
Total_B = count_table[count_table$organ == "Ki" & count_table$age == "Adult" & count_table$chr_class == "autosome", "peak_count"][[1]]+Count_B
contingency_table_Ki_Xa <- matrix(c(Count_A, Total_A - Count_A, Count_B, Total_B - Count_B), nrow = 2)
fisher.test(contingency_table_Ki_Xa)

#Kidney Xi
Count_A = count_table[count_table$organ == "Ki" & count_table$age == "Aged" & count_table$chr_class == "Xi", "peak_count"][[1]]
Total_A = count_table[count_table$organ == "Ki" & count_table$age == "Aged" & count_table$chr_class == "autosome", "peak_count"][[1]]+Count_A
Count_B = count_table[count_table$organ == "Ki" & count_table$age == "Adult" & count_table$chr_class == "Xi", "peak_count"][[1]]
Total_B = count_table[count_table$organ == "Ki" & count_table$age == "Adult" & count_table$chr_class == "autosome", "peak_count"][[1]]+Count_B
contingency_table_Ki_Xi <- matrix(c(Count_A, Total_A - Count_A, Count_B, Total_B - Count_B), nrow = 2)
fisher.test(contingency_table_Ki_Xi)
#peak gain on the aged inactive X chromosome in the Kidney is highly significant: p-value = 1.502e-06

######---- Figure 4e: Rank plot of aged peaks enrichment ----######

#remove Y info and filter for Kidney mean
Rank_plot <-  ATAC_count %>% dplyr::filter(chr != "chrY") %>% dplyr::select(chr, start, end, name, contains("Ki") & contains("mean"))

# Round up only if the decimal part is 0.5 = rounded_up_x1 <- round(x1 + 0.1^8, 0)
Rank_plot$mean_count_Ki_CASTgenome_78w = round(Rank_plot$mean_count_Ki_CASTgenome_78w + 0.1^8, 0)
Rank_plot$mean_count_Ki_BL6genome_78w = round(Rank_plot$mean_count_Ki_BL6genome_78w + 0.1^8, 0)
Rank_plot$mean_count_Ki_CASTgenome_9w = round(Rank_plot$mean_count_Ki_CASTgenome_9w + 0.1^8, 0)
Rank_plot$mean_count_Ki_BL6genome_9w = round(Rank_plot$mean_count_Ki_BL6genome_9w + 0.1^8, 0)

#remove rows where in any mean colum there is a 0 
Rank_plot <- subset(Rank_plot, mean_count_Ki_CASTgenome_78w != 0)
Rank_plot <- subset(Rank_plot, mean_count_Ki_BL6genome_78w != 0)
Rank_plot <- subset(Rank_plot, mean_count_Ki_CASTgenome_9w != 0)
Rank_plot <- subset(Rank_plot, mean_count_Ki_BL6genome_9w != 0)

#set color for X windows to red
Rank_plot$col = "black"
Rank_plot[Rank_plot$chr == "chrX","col"] = "red"

#Calculate enrichment of mean columns
Rank_plot$enrich_Xi <- ((Rank_plot$mean_count_Ki_CASTgenome_78w) / (Rank_plot$mean_count_Ki_CASTgenome_78w + Rank_plot$mean_count_Ki_CASTgenome_9w)) *100
Rank_plot$enrich_Xa <- ((Rank_plot$mean_count_Ki_BL6genome_78w) / (Rank_plot$mean_count_Ki_BL6genome_78w + Rank_plot$mean_count_Ki_BL6genome_9w)) *100

#Binomial test to test significance
for (i in 1:nrow(Rank_plot)) {
  row <- i
  Rank_plot$binom_test_Xi[row] <- pbinom(min(Rank_plot[row,5:6]), Rank_plot[[row,5]]+Rank_plot[[row,6]], 0.5, lower.tail=TRUE, log = FALSE)
  Rank_plot$binom_test_Xa[row] <- pbinom(min(Rank_plot[row,7:8]), Rank_plot[[row,7]]+Rank_plot[[row,8]], 0.5, lower.tail=TRUE, log = FALSE)
}

#plot
Figure_4e <- as.ggplot(expression(
  par(mfcol = c(2,1)), 
  #order/rank matrix accoring to enrich_Xa
  plot_Xa <-  Rank_plot[order(Rank_plot$enrich_Xa),],
  plot(plot_Xa$enrich_Xa,ylab="Fraction of aged peaks [%]",xlab="Rank",main = "BL6 genome",ylim=c(45,100),pch = "_",col=plot_Xa$col),
  #order/rank matrix accorind to enrich_Xi
  plot_Xi <-  Rank_plot[order(Rank_plot$enrich_Xi),],
  plot(plot_Xi$enrich_Xi,ylab="Fraction of aged peaks [%]",xlab="Rank",main = "CAST genome",ylim=c(45,100),pch = "_",col=plot_Xi$col),
  text(x = 1:nrow(plot_Xi)-8,y = plot_Xi$enrich_Xi,labels = ifelse(plot_Xi$enrich_Xi >= 65, paste(plot_Xi$chr,plot_Xi$start/1000000,sep="_"), ""),cex = 0.5)
))

######---- Figure 4h: Bar plot of allelic ratios ----######
#extract median_allelic_ratio
median_allelic_ratios <- ATAC_allelic_ratios
median_allelic_ratios$sample <- gsub("_Rep1","",median_allelic_ratios$sample)
median_allelic_ratios$sample <- gsub("_Rep2","",median_allelic_ratios$sample)
median_allelic_ratios <- median_allelic_ratios %>% dplyr::select(-allelic_ratio, -total_reads)
median_allelic_ratios <- distinct(median_allelic_ratios)


# add colors according to allelic ratio
stacked <- median_allelic_ratios %>%
  mutate(color = ifelse(median_allelic_ratio <= 0.50, "#1c642d",
                        ifelse(median_allelic_ratio > 0.5 & median_allelic_ratio <= 0.55, "#0F7031", 
                               ifelse(median_allelic_ratio > 0.55 & median_allelic_ratio <= 0.6, "#357B30",
                                      ifelse(median_allelic_ratio > 0.6 & median_allelic_ratio <= 0.65, "#4E8330",
                                             ifelse(median_allelic_ratio > 0.65 & median_allelic_ratio <= 0.7, "#658C2D",
                                                    ifelse(median_allelic_ratio > 0.7 & median_allelic_ratio <= 0.75, "#78962A",
                                                           ifelse(median_allelic_ratio > 0.75 & median_allelic_ratio <= 0.8, "#8D9F25",
                                                                  ifelse(median_allelic_ratio > 0.8 & median_allelic_ratio <= 0.85, "#A2A71D",
                                                                         ifelse(median_allelic_ratio > 0.85 & median_allelic_ratio <= 0.9, "#A2A71D",
                                                                                ifelse(median_allelic_ratio > 0.9 & median_allelic_ratio <= 0.95, "#C97314",
                                                                                       ifelse(median_allelic_ratio > 0.95 & median_allelic_ratio <= 1, "#8b1913",""))))))))))))


#calculate percentages
stacked <- stacked %>%
  dplyr::group_by(sample, color) %>%
  dplyr::summarise(count = n())
stacked <- stacked %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(sum_count = sum(count))
stacked <- stacked %>%  dplyr:: mutate(percentage = (count/sum_count)*100)


#set order
stacked$color <- factor(stacked$color, levels = c("#1c642d", "#0F7031", "#357B30", "#4E8330", "#658C2D", "#78962A" ,"#8D9F25", "#A2A71D", "#B3B112", "#C97314" ,"#8b1913"))
stacked$sample <- factor(stacked$sample, levels = c("Ki_9w", "Ki_78w"))

#plot
Figure_4h <- ggplot(stacked, aes(x=sample, fill=color))+
  geom_bar(aes(y = percentage, fill=color), stat = "identity", position = "stack")+
  scale_fill_identity() +
  labs(title = "Figure 4h", x = "", y = "Fraction of ARs on chrX [%]") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))


######---- Figure 4i: Float chart of allelic status ----######
#Extract names of biallelic peaks for subsequent filtering
basis_peaks <- median_allelic_ratios %>% dplyr::filter(organ == "Ki", min_total_reads >= 20, median_allelic_ratio <= 0.9)
basis_peaks <- unique(basis_peaks$name)

#Extract allelic ratios for all peaks in basis_peaks
Float <- median_allelic_ratios %>% dplyr::filter(organ == "Ki", min_total_reads >= 20, name %in% basis_peaks) %>% dplyr::select(-replicate_count, -min_total_reads,-age)

#prepare for plotting
Float <- pivot_wider(Float, names_from = sample, values_from = median_allelic_ratio)
Float <- Float %>% dplyr::rename(Adult = "Ki_9w", Aged = "Ki_78w") %>% dplyr::select(Adult, Aged)
Float <- Float %>% mutate_all(~ ifelse(is.na(.), "NI", 
                                       ifelse(. <= 0.9, "escape", "non_escape")))
Float <- Float[,c("Adult", "Aged")]
Float <- Float %>% 
  dplyr::group_by(Float[,1:2]) %>%
  dplyr::mutate(freq = n()) %>%
  ungroup %>% 
  distinct(Float[,1:2], .keep_all = TRUE) 

Float$Adult <- factor(Float$Adult, levels = c("escape", "non_escape", "NI")) 
Float$Aged <- factor(Float$Aged, levels = c("escape", "non_escape", "NI")) 

#plot
Figure_4i <- as.ggplot(expression(alluvial(Float[,1:2], freq=Float$freq,cex = 0.7)))


######---- Figure 4kl: Identification and characterization of top ∆AR peaks & Regulatory element plot ----######
#Calculate delta for all peaks on chrX that appear [1] in both replicates per timepoint and [2] in both timepoints
#filter step 1: filter for peaks that appear in both replicates
DeltaAR <- median_allelic_ratios %>% dplyr::select(chr, start, end, sample, median_allelic_ratio, name, replicate_count) %>% 
  dplyr::filter(replicate_count == 2) %>% dplyr::select(-replicate_count)
DeltaAR <- pivot_wider(DeltaAR, names_from = sample, values_from = median_allelic_ratio)
#filter step 2: filter for peaks that appear in timepoints
DeltaAR <- DeltaAR[complete.cases(DeltaAR), ]
#Calculate delta
DeltaAR <- DeltaAR %>% dplyr::mutate(delta= DeltaAR$Ki_9w - DeltaAR$Ki_78w)
#Rank ordering
DeltaAR <- DeltaAR[order(DeltaAR$delta),]

#save top ∆AR peaks for subsequent plots
DeltaAR_peaks <- DeltaAR %>% dplyr::filter(delta > 0.099)

#distances to escape genes and their overlapping regulatory elements of top ∆AR peaks were manually documented 
manual_distance <- data.frame(
  chr = c("chrX", "chrX", "chrX", "chrX", "chrX", "chrX", "chrX", "chrX", "chrX", "chrX", 
          "chrX", "chrX", "chrX", "chrX", "chrX", "chrX", "chrX", "chrX", "chrX"),
  start = c(12673033, 13229951, 13236184, 13242587, 19904495, 48256259, 94012604, 102037223, 
            139774435, 139807376, 150458669, 152178844, 160765671, 161961185, 162565341, 
            164068949, 164439567, 168673415, 169879462),
  end = c(12674032, 13231096, 13236799, 13243451, 19905029, 48256956, 94013061, 102037884, 
          139775115, 139808816, 150459389, 152179965, 160766471, 161962072, 162566257, 
          164069526, 164440957, 168674812, 169880087),
  escape_gene = c("1810030O07Rik", "Ddx3x", "Ddx3x", "Ddx3x", "Kdm6a", "Utp14a", "Pdk3", "Jpx", 
                  "Cldn2", "Cldn2", "Iqsec2", "Iqsec2", "Nhs", "Nhs", "Reps2", "Cltrn", "Vegfd", 
                  "Vegfd", "Vegfd"),
  distance = c(0, -50000, -44000, -38000, 1742000, 0, 180000, -1455000, -26000, 7000, 
               -1683000, 0, -1393000, 198000, 0, -19000, 66000, 4300000, 5507000),
  element = c("Promoter", "Enhancer", "Enhancer", "Enhancer", "Enhancer", "Promoter", "Enhancer", 
              "Enhancer", "Enhancer", "Enhancer", "Enhancer", "Promoter", "Enhancer", 
              "Enhancer", "Promoter", "intergenic", "Enhancer", "Enhancer", 
              "Enhancer"),
  stringsAsFactors = FALSE
)


Distance <- merge(DeltaAR_peaks, manual_distance, by = c("chr", "start", "end"))

#assign colors to regulatory elements
Distance <- Distance %>% dplyr::mutate(color = ifelse(element == "Enhancer", '#F8B813', 
                                                      ifelse(element == "Promoter", '#E84B1F',
                                                             ifelse(element == "intergenic", '#B1B0AF',""))))

# create pie plot for Figure 4i (fraction of near and far ∆AR peaks)
near_peaks <- Distance %>% dplyr::filter(distance >= -200000 & distance <= 200000)
near <- nrow(near_peaks)
far_peaks <- Distance %>% dplyr::filter(distance < -200000 | distance > 200000)
far <- nrow(far_peaks)

pie_4k <- data.frame(category= c("Escapee-distance < 200 kb", "Escapee-distance > 200 kb"),
                     number= c(near,far),
                     percentage= c(round(((near/(sum(near,far)))*100), digits = 2), round(((far/(sum(near,far)))*100),digits = 2)))
pie_4k$label <- paste0(pie_4k$category, " ", pie_4k$percentage, "%")

Figure_4k_pie <- ggplot(pie_4k, aes(x = "", y = number, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y") +
  scale_fill_manual(values = c("#532A85", "#BFBDBE")) +
  theme_void() +
  theme(legend.position = "none") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            color= "white") +
  labs(title = "")

#for Rank plot, add colors for distance
DeltaAR <- DeltaAR %>% mutate(color = ifelse(name %in% near_peaks$name & delta > 0.099, '#532A85', 
                                             ifelse(name %in% far_peaks$name & delta > 0.099, '#BFBDBE','black')))
DeltaAR$Rank <- rownames(DeltaAR)
DeltaAR$Rank <- as.numeric(DeltaAR$Rank)

#set colors
distance_colors <- c("#532A85", "#BFBDBE", "black")

#plot
Figure_4k_rank <- ggplot(DeltaAR, aes(x= Rank, y= delta, col = color)) +
  geom_point(size=3)+
  scale_color_manual(values= distance_colors)+
  geom_hline(yintercept=0,  color = "black", lwd=0.5)+
  geom_hline(yintercept=0.1, linetype="dashed", color = "grey", lwd=0.5)+
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.ticks.x=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x = element_text(angle = 0, hjust=0.1),
    axis.line = element_line(colour = 'black', size = 0.2),
    axis.ticks.y=element_line(colour = "black",size=0.2),
    axis.title.y=element_text(size=12),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    legend.position="none")+
  labs(title = "Figure 4k",
       x = "Rank",
       y = "∆AR (Adult - Aged)")

# embed Pie plot in Rank plot
Figure_4k_pie_grob <- ggplotGrob(Figure_4k_pie)
Figure_4k <- Figure_4k_rank +
  annotation_custom(
    Figure_4k_pie_grob,
    xmin = 0,
    xmax = 100,
    ymin = 0.1,
    ymax = 0.5)

#prepare Figure 4l
Distance <- Distance %>% dplyr::filter(-200000 <= distance, distance <= 200000)
Distance$plottingDistance <- Distance$distance /1000

pie_4l <- ddply(Distance,c("element"),nrow)
colnames(pie_4l)[2] <- "number"
pie_4l <- pie_4l %>%
  dplyr::mutate(sum_count = sum(number)) %>% dplyr:: mutate(percentage = round(((number/sum_count)*100),2))

#set colors
element_colors <- c("#F8B915", "#B1B1B0", "#E94C20")

#plot pie 4l
Figure_4l_pie <- ggplot(pie_4l, aes(x = "", y = number, fill = element)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y") +
  scale_fill_manual(values = element_colors) +
  theme_void() +
  theme(legend.position = "none") +
  geom_text(aes(label = percentage),
            position = position_stack(vjust = 0.5),
            color= "black") +
  labs(title = "")

#plot scatter 4l
Figure_4l_scatter <- ggplot(Distance, aes(x= plottingDistance, y= delta, color= element)) +
  geom_point(size=3) +
  scale_color_manual(values= element_colors) +
  scale_x_continuous(limits= c(-200,200), labels = label_comma()) +
  theme_minimal() +
  labs(title = "Figure 4l",
       x = "Distance to escape gene [kb]",
       y = "∆AR (Adult - Aged)",
       color = "")

# embed Pie plot in Rank plot
Figure_4l_pie_grob <- ggplotGrob(Figure_4l_pie)
Figure_4l <- Figure_4l_scatter +
  annotation_custom(
    Figure_4l_pie_grob,
    xmin = 0,
    xmax = -200,
    ymin = 0.3,
    ymax = 0.5)

######---- Figure 4f,g,j: Karyoplot plot ----######

#Figure 4f window analysis
Karyoplot_window <- Rank_plot %>% dplyr::filter(chr== "chrX")
#Figure 4g Escape peaks
Adult_escapePeak <- median_allelic_ratios %>% dplyr::filter(age== "Adult" & median_allelic_ratio <= 0.9)
Aged_escapePeak <- median_allelic_ratios %>% dplyr::filter(age== "Aged" & median_allelic_ratio <= 0.9)
#Escape gene labels
drop <- c("Xist", "Gm14719", "Mid1")
allelic_ratios_Adult <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "b AllelicRatio_Adult")
allelic_ratios_Aged <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "j AllelicRatio_Aged")

allelic_ratios_Adult <- allelic_ratios_Adult %>%  dplyr::filter(sex== "female" & organ == "Ki" & chr== "chrX" & median_allelic_ratio <=0.9 & ! name %in% drop)
allelic_ratios_Aged <- allelic_ratios_Aged %>%  dplyr::filter(sex== "female" & organ == "Ki" & chr== "chrX" & median_allelic_ratio <=0.9 & ! name %in% drop)

label <- rbind(allelic_ratios_Adult, allelic_ratios_Aged)
label <- label %>%  dplyr::select("chr", "start", "end", "name")
label <- distinct(label)
label$chr <- gsub("chr","",label$chr)
genes_coord_all <- regioneR::toGRanges(as.data.frame(label)) 
seqlevelsStyle(genes_coord_all) <- "UCSC"

#Escape genes
allelic_ratios_Adult <- distinct(allelic_ratios_Adult %>% dplyr::select("chr", "start", "end", "name"))
allelic_ratios_Aged <- distinct(allelic_ratios_Aged %>% dplyr::select("chr", "start", "end", "name"))


#set scales
scale_top <- 100
scale_bottom <- 0

#plot
Figure_4fgj <- as.ggplot(expression(
  kp <- plotKaryotype(genome="mm10",chromosomes="chrX",plot.type = 2),
  kpAddBaseNumbers(kp),
  kpAxis(kp,ymin = scale_bottom ,ymax = scale_top),
  kpRect(kp, chr="chrX", x0=Karyoplot_window$start, x1=Karyoplot_window$end, y0=Karyoplot_window$enrich_Xi, y1=Karyoplot_window$enrich_Xi, 
         ymin = scale_bottom,ymax = scale_top, border="darkgrey",data.panel = 1),
  kpRect(kp, chr="chrX", x0=Karyoplot_window$start, x1=Karyoplot_window$end, y0=Karyoplot_window$enrich_Xa, y1=Karyoplot_window$enrich_Xa, 
         ymin = scale_bottom,ymax = scale_top, border="#CD9C25",data.panel = 1),
  kpAddLabels(kp, labels = "Fraction of aged peaks", srt = 90,data.panel = 1, pos=1),
  
  #Escape genes
  kpPlotMarkers(kp, data = genes_coord_all,
                labels = genes_coord_all$name,
                text.orientation = "vertical",
                r1 = 0.5, cex = 0.6,
                ymax=5),
  
  #Escape genes
  kpRect(kp, chr= "chrX", x0=allelic_ratios_Adult$start, x1=allelic_ratios_Adult$end, y0=0, y1=1, 
         ymin = 0,ymax = 1, border="#1E652D",r0=autotrack(1,7),data.panel = 2),
  kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(1,7),data.panel = 2),
  kpAddLabels(kp, labels = "Adult escape genes", r0=autotrack(1,7), label.margin = 0.01,data.panel = 2),
  
  kpRect(kp, chr= "chrX", x0=allelic_ratios_Aged$start, x1=allelic_ratios_Aged$end, y0=0, y1=1, 
         ymin = 0,ymax = 1, border="#1E652D",r0=autotrack(2,7),data.panel = 2),
  kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(2,7),data.panel = 2),
  kpAddLabels(kp, labels = "Aged escape genes", r0=autotrack(2,7), label.margin = 0.01,data.panel = 2),
  
  #Escape peaks
  kpRect(kp, chr= "chrX", x0=Adult_escapePeak$start, x1=Adult_escapePeak$end, y0=0, y1=1, 
         ymin = 0,ymax = 1, border="#1E652D",r0=autotrack(4,7),data.panel = 2),
  kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(4,7),data.panel = 2),
  kpAddLabels(kp, labels = "Adult escape peaks", r0=autotrack(4,7), label.margin = 0.01,data.panel = 2),
  
  kpRect(kp, chr= "chrX", x0=Aged_escapePeak$start, x1=Aged_escapePeak$end, y0=0, y1=1, 
         ymin = 0,ymax = 1, border="#1E652D",r0=autotrack(5,7),data.panel = 2),
  kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(5,7),data.panel = 2),
  kpAddLabels(kp, labels = "Aged escape peaks", r0=autotrack(5,7), label.margin = 0.01,data.panel = 2),
  
  #Top ∆AR peaks
  kpRect(kp, chr= "chrX", x0=DeltaAR_peaks$start, x1=DeltaAR_peaks$end, y0=0, y1=1, 
         ymin = 0,ymax = 1, border="#532A85",r0=autotrack(7,7),data.panel = 2),
  kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(7,7),data.panel = 2),
  kpAddLabels(kp, labels = "∆AR >= 0.1", r0=autotrack(7,7), label.margin = 0.01,data.panel = 2)
))


######---- Plot Figure 4 part 2 ----###### 
grid.arrange(Figure_4d, Figure_4e,Figure_4fgj, Figure_4h,Figure_4i,Figure_4k,Figure_4l,
             layout_matrix = rbind(c(NA,1,2,2),
                                   c(3,3,4,5),
                                   c(3,3,6,6),
                                   c(3, 3, 7,7)))

