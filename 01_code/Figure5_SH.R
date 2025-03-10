##***************************************##
##   R script to generate main Figure 5  ##
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
input <- "/Users/shoelzl/Desktop/mnt_new2/ge43sik2/WP1_aging_SH/manuscript/EscaperPaper/20241216_nature_aging_revision/revised_supplemetary_tables/"

######------ Read in ATAC tables ------###### 
ATAC_count <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "p ATAC peaks (SNPsplit)")
allelic_ratios <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "q ATAC peaks (AllelicRatio)")

# set and apply min_total_reads cutoff
min=20
allelic_ratios <- allelic_ratios %>% dplyr::filter(min_total_reads >= min)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#


######---- Figure 5b: Fractions of ATAC peaks ----######
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
Figure_5b <- ggplot(count_table, aes(x = chr_class, y = perc, fill = age)) +
  geom_bar(stat = "identity", col= "black") +
  facet_wrap(~organ) +
  scale_fill_manual(values= age_colors) +
  labs(title = "Figure 5a",
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

######---- Figure 5c: Rank plot of aged peaks enrichment ----######

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
Figure_5c <- as.ggplot(expression(
par(mfcol = c(2,1)), 
#order/rank matrix accoring to enrich_Xa
plot_Xa <-  Rank_plot[order(Rank_plot$enrich_Xa),],
plot(plot_Xa$enrich_Xa,ylab="Fraction of aged peaks [%]",xlab="Rank",main = "BL6 genome",ylim=c(45,100),pch = "_",col=plot_Xa$col),
#order/rank matrix accorind to enrich_Xi
plot_Xi <-  Rank_plot[order(Rank_plot$enrich_Xi),],
plot(plot_Xi$enrich_Xi,ylab="Fraction of aged peaks [%]",xlab="Rank",main = "CAST genome",ylim=c(45,100),pch = "_",col=plot_Xi$col),
text(x = 1:nrow(plot_Xi)-8,y = plot_Xi$enrich_Xi,labels = ifelse(plot_Xi$enrich_Xi >= 65, paste(plot_Xi$chr,plot_Xi$start/1000000,sep="_"), ""),cex = 0.5)
))

######---- Figure 5f: Bar plot of allelic ratios ----######
#extract median_allelic_ratio
median_allelic_ratios <- allelic_ratios
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
Figure_5f <- ggplot(stacked, aes(x=sample, fill=color))+
  geom_bar(aes(y = percentage, fill=color), stat = "identity", position = "stack")+
  scale_fill_identity() +
  labs(title = "Figure 5f", x = "", y = "Fraction of ARs on chrX [%]") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))


######---- Figure 5g: Float chart of allelic status ----######
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
Figure_5g <- as.ggplot(expression(alluvial(Float[,1:2], freq=Float$freq,cex = 0.7)))


######---- Figure 5ij: Identification and characterization of top ∆AR peaks & Regulatory element plot ----######
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

# create pie plot for Figure 5i (fraction of near and far ∆AR peaks)
near_peaks <- Distance %>% dplyr::filter(distance >= -200000 & distance <= 200000)
near <- nrow(near_peaks)
far_peaks <- Distance %>% dplyr::filter(distance < -200000 | distance > 200000)
far <- nrow(far_peaks)

pie_5i <- data.frame(category= c("Escapee-distance < 200 kb", "Escapee-distance > 200 kb"),
                      number= c(near,far),
                      percentage= c(round(((near/(sum(near,far)))*100), digits = 2), round(((far/(sum(near,far)))*100),digits = 2)))
pie_5i$label <- paste0(pie_5i$category, " ", pie_5i$percentage, "%")

Figure_5i_pie <- ggplot(pie_5i, aes(x = "", y = number, fill = category)) +
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
Figure_5i_rank <- ggplot(DeltaAR, aes(x= Rank, y= delta, col = color)) +
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
  labs(title = "Figure 5i",
       x = "Rank",
       y = "∆AR (Adult - Aged)")

# embed Pie plot in Rank plot
Figure_5i_pie_grob <- ggplotGrob(Figure_5i_pie)
Figure_5i <- Figure_5i_rank +
  annotation_custom(
    Figure_5i_pie_grob,
    xmin = 0,
    xmax = 100,
    ymin = 0.1,
    ymax = 0.5)

#prepare Figure 5j
Distance <- Distance %>% dplyr::filter(-200000 <= distance, distance <= 200000)
Distance$plottingDistance <- Distance$distance /1000

pie_5j <- ddply(Distance,c("element"),nrow)
colnames(pie_5j)[2] <- "number"
pie_5j <- pie_5j %>%
  dplyr::mutate(sum_count = sum(number)) %>% dplyr:: mutate(percentage = round(((number/sum_count)*100),2))

#set colors
element_colors <- c("#F8B915", "#B1B1B0", "#E94C20")

#plot pie 5j
Figure_5j_pie <- ggplot(pie_5j, aes(x = "", y = number, fill = element)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y") +
  scale_fill_manual(values = element_colors) +
  theme_void() +
  theme(legend.position = "none") +
  geom_text(aes(label = percentage),
            position = position_stack(vjust = 0.5),
            color= "black") +
  labs(title = "")
  
#plot scatter 5j
Figure_5j_scatter <- ggplot(Distance, aes(x= plottingDistance, y= delta, color= element)) +
  geom_point(size=3) +
  scale_color_manual(values= element_colors) +
  scale_x_continuous(limits= c(-200,200), labels = label_comma()) +
  theme_minimal() +
  labs(title = "Figure 5j",
       x = "Distance to escape gene [kb]",
       y = "∆AR (Adult - Aged)",
       color = "")
  
# embed Pie plot in Rank plot
Figure_5j_pie_grob <- ggplotGrob(Figure_5j_pie)
Figure_5j <- Figure_5j_scatter +
  annotation_custom(
    Figure_5j_pie_grob,
    xmin = 0,
    xmax = -200,
    ymin = 0.3,
    ymax = 0.5)

######---- Figure 5d,e,h: Karyoplot plot ----######

#Figure 5d window analysis
Karyoplot_window <- Rank_plot %>% dplyr::filter(chr== "chrX")
#Figure 5e Escape peaks
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
Figure_5deh <- as.ggplot(expression(
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


######---- Plot Figure 5 ----###### 
grid.arrange(Figure_5b, Figure_5c,Figure_5deh, Figure_5f,Figure_5g,Figure_5i,Figure_5j,
             layout_matrix = rbind(c(NA,1,2,2),
                                   c(3,3,4,5),
                                   c(3,3,6,6),
                                   c(3, 3, 7,7)))


