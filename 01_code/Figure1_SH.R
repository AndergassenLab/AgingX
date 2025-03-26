##***************************************##
##  R script to generate main Figure 1   ##
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
library(readxl)
library(gridExtra)
library(ggplotify)
library(karyoploteR)
library(eulerr)

######------ Set environment  ------###### 
input <- "/Users/shoelzl/Desktop/mnt_new2/ge43sik2/WP1_aging_SH/manuscript/EscaperPaper/20250321_nature_aging_revision3/"

######------ Read in Allelic Ratios  ------###### 
allelic_ratios <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "b AllelicRatio_Adult")

# set and apply min_total_reads cutoff
min=20
allelic_ratios <- allelic_ratios %>% dplyr::filter(min_total_reads >= min)

# extract median data
allelic_ratios$sample <- gsub("_1","",allelic_ratios$sample)
allelic_ratios$sample <- gsub("_2","",allelic_ratios$sample)
allelic_ratios$sample <- gsub("_3","",allelic_ratios$sample)
allelic_ratios <- allelic_ratios %>% dplyr::select(-"allelic_ratio", -"total_reads")
allelic_ratios <- distinct(allelic_ratios)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#

######------ Read in TMPs  ------###### 
TPM <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "c TPM_Adult")

# extract median data
TPM$sample <- gsub("_1","",TPM$sample)
TPM$sample <- gsub("_2","",TPM$sample)
TPM$sample <- gsub("_3","",TPM$sample)
TPM <- TPM %>% dplyr::select(-"TPM")
TPM <- distinct(TPM)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#

######---- Figure 1b: Violin Plot for Allelic ratios of X linked genes across organs----###### 

# determine male escape genes
escapers_XY <- allelic_ratios%>%dplyr::filter(chr == "chrX" & sex == "male" & median_allelic_ratio<=0.90)
drop <- unique(escapers_XY$name)

Violin <- allelic_ratios %>% dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop)

# set organ colors
organ_colors <- c("#DCB465", "#8B1812","#244C51", "#C97A41", "#A77A76", "#555463", "#97A092")

# plot
Figure_1b <- ggplot(Violin, aes(x=organ, y=median_allelic_ratio, fill=organ)) + 
  geom_violin(width = 1, size = 0.1, scale = "width", col = "black") +
  labs(x="Organ",y="Allelic ratio") +
  geom_boxplot(width = 0.3, size=0.2, alpha = 0.2, color = "white", outlier.shape=NA) +
  scale_fill_manual(values = organ_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0,0.5,0.9,1)) +
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
    legend.position="none") +
  ggtitle("Figure 1b") +
  geom_text(aes(label = ifelse(median_allelic_ratio<=0.90, name,"")), size=2, na.rm = TRUE, hjust = -0.3) +
  geom_point(data= Violin %>% dplyr::filter(median_allelic_ratio<=0.90), size= 1,col="black") +
  geom_hline(yintercept=0.9, linetype="dashed", color = "#8b0000", lwd=0.5)

######---- Figure 1c: Stacked Barplot for informative X-linked genes and their escape status ----###### 

#From now on Xist will not be considered as an escape gene
drop <- c(drop , "Xist")

StackedBar <- allelic_ratios %>% dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop)

StackedBar <- StackedBar %>% dplyr::mutate(escape_status = ifelse(median_allelic_ratio<=0.90, "escape","non-escape"))

#Count X-linked genes
StackedBar <- ddply(StackedBar,c("organ","sex","escape_status"),nrow)
colnames(StackedBar)[4] <- "Xlinked_genes"

#Add percentages for labels
count_table_sum <- StackedBar %>% dplyr::group_by(sex, organ) %>% dplyr::summarise(sum=sum(Xlinked_genes))
StackedBar <- merge(StackedBar, count_table_sum, by = c("sex", "organ"))
StackedBar <- StackedBar %>%  dplyr::mutate(percentage= (Xlinked_genes/sum)*100)
StackedBar$percentage <- round(StackedBar$percentage, 1)

#set barcolors
bar_colors <- c("#1C642D", "#DADADA")

#plot
Figure_1c <- ggplot(StackedBar, aes(x=reorder(organ, -Xlinked_genes), y=Xlinked_genes)) +
  geom_bar(stat="identity",width=0.5, aes(fill=escape_status), col= "black") +
  scale_fill_manual(values = bar_colors) +
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
  geom_text(aes(label=ifelse(escape_status == "escape", percentage,"")), vjust=0, size=4)+
  ggtitle("Figure 1c")

######---- Figure 1d: Barplot for number of escape genes per organ ----###### 

Barplot <- allelic_ratios %>% dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop & median_allelic_ratio<=0.90)

Barplot <- ddply(Barplot,c("organ"),nrow)
colnames(Barplot)[2] <- "Number_of_escape_genes"

Figure_1d <- ggplot(Barplot, aes(x=reorder(organ, -Number_of_escape_genes), y=Number_of_escape_genes, fill = organ)) +
  geom_bar(stat="identity",width=0.5) +
  scale_fill_manual(values = organ_colors) +
  geom_text(aes(label=Number_of_escape_genes), vjust=-0.3, size=4) +
  ylim(0,15) +
  labs(x="",y="Number of escape genes") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.ticks.x=element_line(colour = "black",size=0.2),
    axis.title.x=element_text(size=12),
    axis.text.x = element_text(angle = 0, hjust=0.1),
    axis.line = element_line(colour = 'black', size = 0.2),
    axis.ticks.y=element_line(colour = "black",size=0.2),
    axis.title.y=element_text(size=12),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    legend.position="right")+
  ggtitle("Figure 1d")


######---- Figure 1e: Escape heatmap ----###### 
Escapees <- allelic_ratios %>% dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop & median_allelic_ratio<=0.90)

#Extract escape genes to plot in heatmap
Escapees <- unique(Escapees$name)

#Extract allelic ratios for escape genes
Heatmap <- allelic_ratios %>% dplyr::filter(sex == "female" & name %in% Escapees)
Heatmap <- Heatmap %>% dplyr::select(name, median_allelic_ratio, organ)

#Bring into heatmap format
Heatmap <- Heatmap %>% pivot_wider(names_from = organ, values_from = median_allelic_ratio)
Heatmap <- Heatmap %>% remove_rownames %>% column_to_rownames(var="name")

# define Row and Column order
Heatmap <-  Heatmap[c("Kdm6a","Kdm5c","Ddx3x","Eif2s3x","Ftx","Jpx", "Pbdc1","5530601H04Rik",
                                       "4930578C19Rik","Plp1","Slc16a2","Utp14a", "Smpx", "Avpr2","Cfp","Cybb",
                                       "Tlr8","Xpnpep2","Dmrtc1a","LOC108167322",
                                       "Nrk"),]
Heatmap <- Heatmap[, c("Br", "He", "Ki", "Li", "Lu", "Mu", "Sp")]
Heatmap <- as.matrix(Heatmap)


#set colors and breaks
my_colors <- c("#1c642d", "#0F7031", "#357B30", "#4E8330", "#658C2D", "#78962A" ,"#8D9F25", "#A2A71D", "#B3B112", "#C97314" ,"#8b1913")
my_breaks <- seq(0.5,1,by = 0.05)

#plot
Figure_1e <- as.ggplot(function() {heatmap.2(Heatmap, Rowv = NA, Colv = NA, tracecol=NULL, sepcolor="black", keysize=1, dendrogram = "none",
                        cexCol=1,cexRow = 1,na.color = "white", breaks = my_breaks,
                        sepwidth=c(0.005,0.01), colsep=1:ncol(Heatmap), rowsep=1:nrow(Heatmap),
                        col=colorRampPalette(my_colors),main = "Figure 1e")})


######---- Figure 1f: Correlation Escape - TPM fold change ----###### 

# Calculate fold change (pseudocount)
Foldchange <- TPM %>% dplyr::select(-"sample")

Foldchange$medianTPM <- Foldchange$medianTPM + 1
Foldchange <- Foldchange%>%pivot_wider(names_from = sex, values_from = medianTPM)
Foldchange$FC <-(Foldchange$female) / (Foldchange$male)
Foldchange=Foldchange[complete.cases(Foldchange), ]
Foldchange <- Foldchange %>%  dplyr::filter(male != 1, female != 1)

# Extract escape genes
Escapees <- allelic_ratios %>% dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop & median_allelic_ratio<=0.90) %>% 
  dplyr::select("name", "organ", "median_allelic_ratio")

Correlation <- join(Foldchange,Escapees, by=c("name","organ"))
Correlation=Correlation[complete.cases(Correlation), ]


# Perform the Shapiro-Wilk test for normality on both variables
shapiro.test(Correlation$median_allelic_ratio)
shapiro.test(Correlation$FC)
## pvalue < 0.05 --> no normal distribution --> Spearman

# caluclate Correlation for all organs combined and individually
overall <- cor(Correlation$FC, Correlation$median_allelic_ratio, method = "spearman")

organs <- c("He", "Sp", "Br", "Li", "Lu", "Ki", "Mu")
# Loop through each organ
for (organ in organs) {
  # filter for organ
  table_correlation <- Correlation %>% dplyr::filter(organ == !!organ)
  # calculate correlation
  correlation <- cor(table_correlation$FC, table_correlation$median_allelic_ratio, method = "spearman")
  # assign correlation to organ
  assign(paste0("correlation_", organ), correlation)
}

# plot
Figure_1f <- ggplot(Correlation, aes(x=median_allelic_ratio, y=FC, aes(col=organ))) +
  geom_point(aes(col=organ), size = 2) +
  theme(aspect.ratio=1,
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_line(colour = "black",size=0.2),
        axis.title.x=element_text(size=12),
        axis.text.x = element_text(angle = 0, hjust=0.1),
        axis.line = element_line(colour = 'black', size = 0.2),
        axis.ticks.y=element_line(colour = "black",size=0.2),
        axis.title.y=element_text(size=12),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position="right") +
  ylab("TPM fold change (female vs. male)") +
  xlab("Allelic ratio") +
  xlim(0.5,0.92) +
  ylim(0.5,2.0) +
  scale_color_manual(values = organ_colors) +
  geom_smooth(method="lm", aes(), se=FALSE,size = 0.5, color = "black") +
  annotate("text", x=0.55, y=1.2, label= paste("overall R =",round(overall,digits = 4))) +  
  annotate("text", x=0.55, y=1.1, label= paste("He R =",round(correlation_He,digits = 4))) +
  annotate("text", x=0.55, y=1.0, label= paste("Sp R =",round(correlation_Sp,digits = 4))) +
  annotate("text", x=0.55, y=0.9, label= paste("Lu R =",round(correlation_Lu,digits = 4))) +
  annotate("text", x=0.55, y=0.8, label= paste("Br R =",round(correlation_Br,digits = 4))) +
  annotate("text", x=0.55, y=0.7, label= paste("Ki R =",round(correlation_Ki,digits = 4))) +
  annotate("text", x=0.55, y=0.6, label= paste("Mu R =",round(correlation_Mu,digits = 4))) +
  annotate("text", x=0.55, y=0.5, label= paste("Li R =",round(correlation_Li,digits = 4))) +
  ggtitle("Figure 1f") 


  

######---- Figure 1g: Boxplot TPM fold change ----###### 

# Add escape status to fold changes
Escape_info <- allelic_ratios %>% dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop & median_allelic_ratio<=0.90)
Escape_info <- Escape_info%>%dplyr::select("name","organ")
Escape_info$escape_status <- "escape"
Boxplot <- join(Foldchange,Escape_info, by=c("name","organ"))
Boxplot[is.na(Boxplot )] <- "non-escape"

# exclude Xist
Boxplot <- Boxplot %>% dplyr::filter(name != "Xist")

#set comparisons and order
my_comparisons <- list(c("non-escape", "escape"))
Boxplot$escape_status <- factor(Boxplot$escape_status, levels = c("non-escape", "escape"))

# plot
Figure_1g <- ggplot(Boxplot, aes(x=escape_status, y=FC, fill=escape_status)) +
  geom_boxplot(size = 0.5,outlier.size = 0.5,lwd=0.5) +
  facet_wrap(~organ, ncol = 7) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.ticks.x=element_line(colour = "black",size=0.2),
    axis.title.x=element_text(size=12),
    axis.text.x = element_text(angle = 90, hjust=0.1),
    axis.line = element_line(colour = 'black', size = 0.2),
    axis.ticks.y=element_line(colour = "black",size=0.2),
    axis.title.y=element_text(size=12),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    legend.position="right") +
  scale_fill_manual(values = c("#DADADA", "#1C642D")) +
  labs(x="",y="TPM fold change (female vs. male)") +
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", aes(label=..p.format..),label.y =2) +
  ggtitle("Figure 1g") #+
  #ylim(0,2) 


######---- Figure 1h: Functional roles of escape genes ----###### 
Function <- read_excel(paste0(input,"SupplementaryTable_2.xlsx"))
Function <- Function[-1, ]

Function <- Function[grep("Adult", Function$found_in), ]
Function <- Function %>%  dplyr::select("Gene_symbol", "Function")
Function <- Function %>%
  mutate(Function = sub("lncRNA.*", "lncRNA", Function))
Function <- Function %>%
  mutate(Function = sub("neuronal_maintenance.*", "neuronal_maintenance", Function))

# count functional roles of escape genes
Function <- ddply(Function,c("Function"),nrow)
colnames(Function)[2] <- "count"

# calculate percentage
Function$Perc <- (Function$count / sum(Function$count))*100

# reorder
Function$Function <- factor(Function$Function, levels = c("unknown_function", "neuronal_maintenance", "muscle_structure", "translational_regulator",
                                                      "ribosome_biogenesis","proteolysis", "hormone_transporter", "hormone_receptor", "kinase", "transcriptional_regulator", "epigenetic_regulator", "immune_system", "lncRNA"))
# set colors
function_colors <- c("#91CEC4", "#FBF4B5", "#BEBADA", "#EF7E73", "#80B0D3", "#F8B366", "#8FA068", "#B2D06D", "#F9CDE2", "#D9D9D9", "#BB80B6", "#CCE3C3", "#FFEE70")

# plot
Figure_1h <- ggplot(Function, aes(x="", y=Perc, fill=Function)) +
  geom_bar(stat = "identity") +
  labs(title = "Figure 1h",
       x = "",
       y = "Function [%]") +
  theme_minimal()+
  scale_fill_manual(values= function_colors) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10))


######---- Figure 1i: Karyoplot for all escape genes across organs ----###### 
Karyoplot <- allelic_ratios %>% dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop & median_allelic_ratio<=0.90)

# get coordinates for all escapees
label <- Karyoplot %>%  dplyr::select("chr", "start", "end", "name")
label <- distinct(label)
label$chr <- gsub("chr","",label$chr)
genes_coord_all <- regioneR::toGRanges(as.data.frame(label)) 
seqlevelsStyle(genes_coord_all) <- "UCSC"

# separate organs
Karyoplot_He <- Karyoplot %>%  dplyr::filter(organ == "He")
Karyoplot_Br <- Karyoplot %>%  dplyr::filter(organ == "Br")
Karyoplot_Li <- Karyoplot %>%  dplyr::filter(organ == "Li")
Karyoplot_Lu <- Karyoplot %>%  dplyr::filter(organ == "Lu")
Karyoplot_Sp <- Karyoplot %>%  dplyr::filter(organ == "Sp")
Karyoplot_Mu <- Karyoplot %>%  dplyr::filter(organ == "Mu")
Karyoplot_Ki <- Karyoplot %>%  dplyr::filter(organ == "Ki")

# plot
# chromosome ideogram
Figure_1i <- as.ggplot(expression(
  kp <- plotKaryotype(genome="mm10",chromosomes="chrX",plot.type = 2),
  # add scale
  kpAddBaseNumbers(kp),
  # add gene names
  kpPlotMarkers(kp, data = genes_coord_all,
                labels = genes_coord_all$name,
                text.orientation = "vertical",
                r1 = 0.5, cex = 0.9,
                ymax=5),
  # add "kpRect" for escape genes for all organs
  kpRect(kp, chr= "chrX", x0=Karyoplot_He$start, x1=Karyoplot_He$end, y0=0, y1=1, 
         ymin = 0,ymax = 1, border="#8B1812",r0=autotrack(7,7),data.panel = 2),
  kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(7,7),data.panel = 2),
  kpAddLabels(kp, labels = "Heart", r0=autotrack(7,7), label.margin = 0.01,data.panel = 2),
  kpRect(kp, chr= "chrX", x0=Karyoplot_Sp$start, x1=Karyoplot_Sp$end, y0=0, y1=1, 
         ymin = 0,ymax = 1, border="#97A092",r0=autotrack(6,7),data.panel = 2),
  kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(6,7),data.panel = 2),
  kpAddLabels(kp, labels = "Spleen", r0=autotrack(6,7), label.margin = 0.01,data.panel = 2),
  kpRect(kp, chr= "chrX", x0=Karyoplot_Br$start, x1=Karyoplot_Br$end, y0=0, y1=1, 
         ymin = 0,ymax = 1, border="#DCB465",r0=autotrack(5,7),data.panel = 2),
  kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(5,7),data.panel = 2),
  kpAddLabels(kp, labels = "Brain", r0=autotrack(5,7), label.margin = 0.01,data.panel = 2),
  kpRect(kp, chr= "chrX", x0=Karyoplot_Lu$start, x1=Karyoplot_Lu$end, y0=0, y1=1, 
         ymin = 0,ymax = 1, border="#A77A76",r0=autotrack(4,7),data.panel = 2),
  kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(4,7),data.panel = 2),
  kpAddLabels(kp, labels = "Lung", r0=autotrack(4,7), label.margin = 0.01,data.panel = 2),
  kpRect(kp, chr= "chrX", x0=Karyoplot_Li$start, x1=Karyoplot_Li$end, y0=0, y1=1, 
         ymin = 0,ymax = 1, border="#C97A41",r0=autotrack(3,7),data.panel = 2),
  kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(3,7),data.panel = 2),
  kpAddLabels(kp, labels = "Liver", r0=autotrack(3,7), label.margin = 0.01,data.panel = 2),
  kpRect(kp, chr= "chrX", x0=Karyoplot_Ki$start, x1=Karyoplot_Ki$end, y0=0, y1=1, 
         ymin = 0,ymax = 1, border="#244C51",r0=autotrack(2,7),data.panel = 2),
  kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(2,7),data.panel = 2),
  kpAddLabels(kp, labels = "Kidney", r0=autotrack(2,7), label.margin = 0.01,data.panel = 2),
  kpRect(kp, chr= "chrX", x0=Karyoplot_Mu$start, x1=Karyoplot_Mu$end, y0=0, y1=1, 
         ymin = 0,ymax = 1, border="#555463",r0=autotrack(1,7),data.panel = 2),
  kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(1,7),data.panel = 2),
  kpAddLabels(kp, labels = "Muscle", r0=autotrack(1,7), label.margin = 0.01,data.panel = 2),
  # add title
  kpAddMainTitle(kp, main= "Figure 1i",col="black")
))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------#

######------ Read in Allelic Ratios Cardiac Cell Types (CCT) ------###### 
allelic_ratios_CCT <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "d AllelicRatio_cardiac_celltype")

# set and apply min_total_reads cutoff
min_CCT=10
allelic_ratios_CCT <- allelic_ratios_CCT %>% dplyr::filter(min_total_reads >= min_CCT)

# extract median data
allelic_ratios_CCT$sample <- gsub("_1","",allelic_ratios_CCT$sample)
allelic_ratios_CCT$sample <- gsub("_2","",allelic_ratios_CCT$sample)
allelic_ratios_CCT$sample <- gsub("_3","",allelic_ratios_CCT$sample)
allelic_ratios_CCT <- allelic_ratios_CCT %>% dplyr::select(-"allelic_ratio", -"total_reads")
allelic_ratios_CCT <- distinct(allelic_ratios_CCT)

######------ Read in TMPs  ------###### 
TPM_CCT <- read_excel(paste0(input,"SupplementaryTable_1.xlsx"), sheet = "e TPM_cardiac_celltype")

# extract median data
TPM_CCT$sample <- gsub("_1","",TPM_CCT$sample)
TPM_CCT$sample <- gsub("_2","",TPM_CCT$sample)
TPM_CCT$sample <- gsub("_3","",TPM_CCT$sample)
TPM_CCT <- TPM_CCT %>% dplyr::select(-"TPM")
TPM_CCT <- distinct(TPM_CCT)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#


######---- Figure 1j: CCT Escape heatmap----###### 

# determine CCT male escape genes
escapers_XY <- allelic_ratios_CCT %>%dplyr::filter(chr == "chrX" & sex == "male" & median_allelic_ratio<=0.90)
drop_CCT <- c(unique(escapers_XY$name), "Xist", "Tsix")
  
escapees_He <- allelic_ratios %>% dplyr::filter(organ == "He" & sex == "female" & chr == "chrX" & ! name %in% drop & median_allelic_ratio<=0.90)
escapees_CCT <- allelic_ratios_CCT %>% dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop_CCT & median_allelic_ratio<=0.90)

# Extract escape genes to plot in heatmap
escapees_He <- unique(escapees_He$name)
escapees_CCT <- unique(escapees_CCT$name)
Escapees <- unique(c(escapees_He, escapees_CCT))
  
# Extract allelic ratios for escape genes
# whole heart
Heatmap_He <- allelic_ratios %>% dplyr::filter(sex == "female" & name %in% Escapees & organ == "He")
Heatmap_He <- Heatmap_He %>% dplyr::select(name, median_allelic_ratio, organ)
# CCT
Heatmap_CCT <- allelic_ratios_CCT %>% dplyr::filter(sex == "female" & name %in% Escapees)
Heatmap_CCT <- Heatmap_CCT %>% dplyr::select(name, median_allelic_ratio, celltype)
colnames(Heatmap_CCT)[colnames(Heatmap_CCT) == "celltype"] <- "organ"
  
# combine whole heart and CCT
esc_Heatmap <- rbind(Heatmap_He, Heatmap_CCT)
  
#prepare data frame for plotting
esc_Heatmap <- esc_Heatmap %>% pivot_wider(names_from = organ, values_from = median_allelic_ratio)
esc_Heatmap <- esc_Heatmap %>% column_to_rownames(var="name")
  
# arrange row and column order
esc_Heatmap <- esc_Heatmap[, c("He", "CM", "MP", "FB", "EC")]
esc_Heatmap <- esc_Heatmap[ c("Kdm5c", "Kdm6a", "Ddx3x", "Pbdc1", "Utp14a", "5530601H04Rik","9530027J09Rik", "Plp1", "Eif2s3x", "Jpx", "Smpx",  "Ftx", "4930578C19Rik",
                                "Cask","Iqsec2", "Usp11","Gpr34", "Gm15232", "1810030O07Rik", "Med14", "Slc16a2"),]
esc_Heatmap <- as.matrix(esc_Heatmap)
# set color gradient
my_colors <- c("#1c642d", "#1c642d", "#0F7031", "#357B30", "#4E8330", "#658C2D", "#78962A" ,"#8D9F25", "#A2A71D", "#B3B112", "#C97314" ,"#8b1913")
my_breaks <- seq(0.5,1,by = 0.05)
  
# plot
Figure_1j <- as.ggplot(function() {heatmap.2(esc_Heatmap, Rowv = NA, Colv = NA, tracecol=NULL, sepcolor="black", keysize=1, dendrogram = "none",
                                               cexCol=1,cexRow = 1,na.color = "white", breaks = my_breaks,
                                               sepwidth=c(0.005,0.01), colsep=1:ncol(esc_Heatmap), rowsep=1:nrow(esc_Heatmap),
                                               col=colorRampPalette(my_colors),main = "Figure 1j")})


######---- Figure 1k: Karyoplot for all escape genes across CCT ----###### 
Karyoplot_CCT <- allelic_ratios_CCT %>% dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop_CCT & median_allelic_ratio<=0.90)
Karyoplot_He <- allelic_ratios %>% dplyr::filter(organ == "He" & sex == "female" & chr == "chrX" & ! name %in% drop & median_allelic_ratio<=0.90)

# get coordinates for all escapees
label_CCT <- Karyoplot_CCT %>%  dplyr::select("chr", "start", "end", "name")
label_He <- Karyoplot_He %>% dplyr::select("chr", "start", "end", "name")
label <- rbind(label_CCT, label_He)
label <- distinct(label)
label$chr <- gsub("chr","",label$chr)
genes_coord_all <- regioneR::toGRanges(as.data.frame(label)) 
seqlevelsStyle(genes_coord_all) <- "UCSC"

# separate celltypes
Karyoplot_CM <- Karyoplot_CCT %>%  dplyr::filter(celltype == "CM")
Karyoplot_FB <- Karyoplot_CCT %>%  dplyr::filter(celltype == "FB")
Karyoplot_EC <- Karyoplot_CCT %>%  dplyr::filter(celltype == "EC")
Karyoplot_MP <- Karyoplot_CCT %>%  dplyr::filter(celltype == "MP")

# plot
# chromosome ideogram
Figure_1k <- as.ggplot(expression(
  kp <- plotKaryotype(genome="mm10",chromosomes="chrX",plot.type = 2),
  # add scale
  kpAddBaseNumbers(kp),
  # add gene names
  kpPlotMarkers(kp, data = genes_coord_all,
                labels = genes_coord_all$name,
                text.orientation = "vertical",
                r1 = 0.5, cex = 0.9,
                ymax=5),
  # add "kpRect" for escape genes for all organs
  kpRect(kp, chr= "chrX", x0=Karyoplot_He$start, x1=Karyoplot_He$end, y0=0, y1=1, 
         ymin = 0,ymax = 1, border="#8B1812",r0=autotrack(1,5),data.panel = 2),
  kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(1,5),data.panel = 2),
  kpAddLabels(kp, labels = "Heart", r0=autotrack(1,5), label.margin = 0.01,data.panel = 2),
  
  kpRect(kp, chr= "chrX", x0=Karyoplot_CM$start, x1=Karyoplot_CM$end, y0=0, y1=1, 
         ymin = 0,ymax = 1, border="#E69A8F",r0=autotrack(2,5),data.panel = 2),
  kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(2,5),data.panel = 2),
  kpAddLabels(kp, labels = "CM", r0=autotrack(2,5), label.margin = 0.01,data.panel = 2),
  kpRect(kp, chr= "chrX", x0=Karyoplot_MP$start, x1=Karyoplot_MP$end, y0=0, y1=1, 
         ymin = 0,ymax = 1, border="#9B95C9",r0=autotrack(3,5),data.panel = 2),
  kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(3,5),data.panel = 2),
  kpAddLabels(kp, labels = "MP", r0=autotrack(3,5), label.margin = 0.01,data.panel = 2),
  kpRect(kp, chr= "chrX", x0=Karyoplot_FB$start, x1=Karyoplot_FB$end, y0=0, y1=1, 
         ymin = 0,ymax = 1, border="#9D9D9C",r0=autotrack(4,5),data.panel = 2),
  kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(4,5),data.panel = 2),
  kpAddLabels(kp, labels = "FB", r0=autotrack(4,5), label.margin = 0.01,data.panel = 2),
  kpRect(kp, chr= "chrX", x0=Karyoplot_EC$start, x1=Karyoplot_EC$end, y0=0, y1=1, 
         ymin = 0,ymax = 1, border="#7D9DCC",r0=autotrack(5,5),data.panel = 2),
  kpRect(kp, chr="chrX", x0=0, x1=171031299, y0=0, y1=1, border="black",r0=autotrack(5,5),data.panel = 2),
  kpAddLabels(kp, labels = "EC", r0=autotrack(5,5), label.margin = 0.01,data.panel = 2),
  # add title
  kpAddMainTitle(kp, main= "Figure 1k",col="black")
))

######---- Figure 1l: Venn diagram for whole heart and CCT escape genes ----###### 

#create euler
Venn_He_CCT <- list(escapees_He = escapees_He, escapees_CCT = escapees_CCT)

#plot
Figure_1l <- plot(euler(Venn_He_CCT), shape = "circle", quantities = TRUE, input = "disjoint", fills = "transparent", main = "Figure 1l")


######---- Figure 1m: Distance plot for CCT and whole He escapees ----###### 
Distance_CCT <- allelic_ratios_CCT %>% dplyr::filter(sex == "female" & chr == "chrX" & ! name %in% drop_CCT & median_allelic_ratio<=0.90)
Distance_He <- allelic_ratios %>% dplyr::filter(organ == "He" & sex == "female" & chr == "chrX" & ! name %in% drop & median_allelic_ratio<=0.90)

Distance_He <- Distance_He %>%  dplyr::select(name, start, end)
Distance_CCT <- Distance_CCT %>%  dplyr::select(name, start, end)
Distance_CCT <- distinct(Distance_CCT)
Distance_CCT <- Distance_CCT %>%dplyr:: filter(! name %in% Distance_He$name)


# Function to calculate shortest distance
calculate_shortest_distance <- function(CCT_row, Distance_He) {
  distances_start <- abs(CCT_row$start - Distance_He$start)
  distances_end <- abs(CCT_row$end - Distance_He$end)
  distances_3 <- abs(CCT_row$start - Distance_He$end)
  distances_4 <- abs(CCT_row$end - Distance_He$start)
  min_distance <- min(distances_start, distances_end, distances_3, distances_4)
  return(min_distance)
}

# Create a new column to store the shortest distance in Distance_CCT
Distance_CCT$distance <- NA

# LOOP to calculate the shortest distance
for (i in 1:nrow(Distance_CCT)) {
  #loop through rows of CCT-specific escapees
  CCT_row <- Distance_CCT[i, ]
  #calculate the shortest distance with formula above and loop through each heart row
  shortest_distances <- sapply(1:nrow(Distance_He), function(j) {
    calculate_shortest_distance(CCT_row, Distance_He[j, ])
  })
  #define which values of the 12 rows is the smallest and assign to new column
  shortest_distance <- min(shortest_distances)
  Distance_CCT$distance[i] <- shortest_distance
}

# filter for distance <= 2500000
Distance_CCT <- Distance_CCT %>% dplyr::arrange(Distance_CCT$distance) %>% dplyr::filter(distance <= 2500000) #8/9 are <= 2500000: use later for pie plot
# manually add direction and closest whole heart escape gene
Distance_CCT$position <- c("positive", "negative", "negative", "positive", "positive", "negative", "negative", "positive")
Distance_CCT$He_escapee <- c("Utp14a", "Kdm5c", "Pbdc1","Ddx3x", "Ddx3x", "Ddx3x", "Ddx3x", "4930578C19Rik" )


Distance_CCT <- Distance_CCT %>%
  mutate(plot_distance = ifelse(grepl('negative', position), (distance*(-1)), 
                                ifelse(grepl('positive', position), distance, ''))) %>% dplyr::arrange(distance)

# prepare for plotting
Distance_CCT$scatter <- rownames(Distance_CCT)
Distance_CCT$plot_distance <- as.numeric(Distance_CCT$plot_distance)
Distance_CCT$scatter <- as.numeric(Distance_CCT$scatter)
Distance_CCT$label <- paste0(Distance_CCT$He_escapee,"_",Distance_CCT$name,"_",Distance_CCT$plot_distance)

# plot
Figure_1m_scatter <- ggplot(Distance_CCT, aes(x= plot_distance, y= scatter)) +
  geom_point(position = "identity") +
  scale_x_continuous(limits= c(-2500000,2500000))+
  scale_y_continuous(limits= c(0,9))+
  labs(title = "Figure 1m", x = "Distance between CCT-specific and Heart escapees [Mb]", y = "") +
  geom_text(aes(label = label))+
  theme_minimal()

# create pie chart
pie <- data.frame(
  category = c(">2.5MB_apart", "<2.5Mb_apart"),
  value = c(1, 8))
# add percentage
pie <- pie %>%
  mutate(percentage = value / sum(value) * 100)

# set colors
pie_colors <- c("#DADADA", "white")

# plot pie
Figure_1m_pie <- ggplot(pie, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y") +
  scale_fill_manual(values = pie_colors) +
  theme_void() +
  theme(legend.position = "none") +
  geom_text(aes(label = paste0(round(percentage), "%")),
            position = position_stack(vjust = 0.5)) +
  labs(title = "")

# embed Pie plot in Distance plot
Figure_1m_pie_grob <- ggplotGrob(Figure_1m_pie)
Figure_1m <- Figure_1m_scatter +
  annotation_custom(
    Figure_1m_pie_grob,
    xmin = -2600000,
    xmax = -800000,
    ymin = 1.8,
    ymax = 6)

Figure_1b
Figure_1c
Figure_1d
Figure_1e
Figure_1f
Figure_1g
Figure_1h
Figure_1i
Figure_1j
Figure_1k
Figure_1l
Figure_1m

