library(tidyverse)
library(MSstatsTMT)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(ggrepel)
library(ggforce)
library(ggfortify)

#### Load Input Arguments ####
source <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[1])) # get the script file location
args <- commandArgs(trailingOnly = TRUE) # get all other arguments

# the below row is for manually loading arguments and going through the script interactively
source <- "/home/emil/projects/pisa"
args <- c("20230216_ascorbate_arabidopsis", 3, 0.05, 8)

#### Specify input information ####
exp <- args[1] # the first argument is the path to the folder containing the relevant data files. The rest of the script presume that the folder name format is "date_metabolite_organism", e.g. 20240312_ppGpp_pcc7942
nbr_conc <- as.numeric(args[2]) # the second argument is the number of concentrations in the PISA experiment, excluding the blank, typically 1 or 3
cutoff <- as.numeric(args[3]) # the third argument is the desired adjusted p-value cutoff for statistical significance
outliers <- as.numeric(args[4:length(args)]) # if any clear outliers have been identified they can be added as arguments one by one after the third to exclude them from the analysis
date <- str_split(exp, "_")[[1]][1] # extract the date from the path argument
metabolite <- str_split(exp, "_")[[1]][2] # extract the metabolite used from the path argument
organism <- str_split(exp, "_")[[1]][3] # extract the organism name from the path argument

path <- paste(source, "data", organism, exp, sep = "/") # paste the path to the folder with input data
uniprot_path <- paste(source, "resources", organism, "", sep = "/") # paste the path to the uniprot annotation file. This script presumes an organization of one folder per organism, but it can be changed if desired

#### Load Data ####
save_dir <- paste(source, "results", sep = "/") # specify the directory in which files are to be saved.
runs <- unlist(list.files(path, pattern = ".raw")) # list all individual raw files in the path folder

uniprot_files <- unlist(list.files(uniprot_path, pattern = paste("uniprot_", organism, "_([0-9]*)\\.tab", sep = ""))) # list all files that follows the format "uniprot_organism_date.tab"

uniprot <- read_tsv(paste(uniprot_path, uniprot_files[order(uniprot_files, decreasing = TRUE)][1], sep = "")) # paste the path with the latest file and load

evidence <- read_tsv(paste(path, "/combined/txt/evidence.txt", sep = "")) # load the evidence.txt output file from maxquant
proteinGroups <-  read_tsv(paste(path, "/combined/txt/proteinGroups.txt", sep = "")) # load the proteinGroups.txt output file from maxquant
annotation.mq <- read_csv(paste(path, "/annotation.mq", nbr_conc, ".csv", sep = "")) # load the annotation.mq input file from maxquant
annotation.mq$Run <- rep(substr(runs, 1, nchar(runs) - 4), each = 16)

#### Convert MaxQuant output to MSstats input ####
input.mq <- MaxQtoMSstatsTMTFormat(evidence, proteinGroups, annotation.mq) # apply the MSstats function to convert maxquant output to MSstats input
input.mq <- subset(input.mq, Condition != "Empty") # if not all TMT channels are used (number of concentrations in experiment are less than 3), adjust unused channels to empty
input.mq <- separate_rows(input.mq, ProteinName, sep = ";") # separate rows for ambiguous peptides. In the output file, the protein that ambigiuous peptides map to are annotated as "protein_a;protein_b" 
input.mq$ProteinName <- gsub(";", "", input.mq$ProteinName) # remove the residual ";" from the protein strings. The result is multiple rows for each ambiguous peptide, one for each protein it maps to

#### outlier removal ####
input.mq$Channel <- as.numeric(gsub(".*?([0-9]+).*", "\\1", input.mq$Channel)) # turn channel column entirely numerical
input.mq <- input.mq %>% filter(!(Channel %in% outliers)) # remove any channels in input argument
input.mq$Channel <- paste("channel", input.mq$Channel, sep = "") # paste back the channel text

#### Perform protein summarization and normalization ####
quant <- proteinSummarization(input.mq, reference_norm = FALSE) # apply protein summarization function from MSstats

#### Perform group comparison ####
comp.mtx <- cbind(matrix(rep(-1, nbr_conc)),diag(nbr_conc)) # create a comparison matrix that compares all metabolite concentrations against the blank
colnames(comp.mtx) <- levels(quant$ProteinLevelData$Condition) # name the columns of the comparison matrix with the appropriate condition numbers
rownames(comp.mtx) <- levels(quant$ProteinLevelData$Condition)[-1] # name the rows of the comparison matrix with the appropriate condition numbers

comp_result <- groupComparisonTMT(quant, contrast.matrix = comp.mtx)[[1]] # apply MSstats group comparison function to calculate significance and fold changes

#### Combine result with uniprot file ####
comp_result$Protein <- unlist(lapply(1:nrow(comp_result), function(x) { # adjust uniprot code format so that results can be merged with uniprot annotation file
                                entry <- str_split(comp_result$Protein[x], "\\|")[[1]][2]
                              }))

colnames(comp_result)[1] <- "Entry" # change the name of the uniprot code column to the same as in the uniprot annotation file

comp_result <- merge(comp_result, uniprot, by = "Entry") # merge the result file with the uniprot annotation file

#### Annotate significant proteins ####
comp_result <- comp_result %>% mutate(Sign = case_when(adj.pvalue < cutoff ~ "Sign", TRUE ~ "Unsign")) # create a new column with values "Sign" or "Unsign" depending on adjusted p value and cutoff argument

#### P-value Distribution ####
pval_dist <- map(.x = levels(quant$ProteinLevelData$Condition)[-1], .f = \(bb){ # plot the p value distribution for each comparison
                df <- subset(comp_result, Label == bb)
                ggplot(df, aes(x = pvalue)) + geom_histogram(bins = 100) + scale_x_continuous(breaks = seq(0, 1, by = 0.2)) + theme_bw() + ylab(paste("conc = ", bb, sep = ""))
             })

if (nbr_conc == 1) {
  pval_dist <- grid.arrange(pval_dist[[1]], nrow = 1, top = "P-value Distribution") # create a plot with just one distribution if number of concentrations is one
}

if (nbr_conc == 3) {
  pval_dist <- grid.arrange(pval_dist[[1]], pval_dist[[2]], pval_dist[[3]], nrow = 3, top =  "P-value Distribution") # create a stacked plot of three distributions if number of concentrations is three
}

#### QQ ####
qq_plot <- ggplot(input.mq, aes(sample = log2(Intensity), colour = Channel)) + stat_qq(size = 0.5) + theme_bw() + geom_abline(intercept = 13, slope = 1) + 
  ggtitle("QQ-Plot") + xlab("Theoretical Quantiles") + ylab("Sample Quantiles") # make a qq plot from untransformed input data for qc purposes

#### Intensity distribution curve ####
dist_plot <- map(.x = list(input.mq), .f = \(cc){
  df <- cc %>% group_by(Channel) %>% group_split() # split the untransformed raw data by TMT channel
  df <- map(.x = df, .f = \(dd){
    tmp <- dd[order(dd$Intensity, na.last = FALSE),] # arrange all peptides by intensity
    tmp <- tibble::rowid_to_column(tmp, "Fraction") # add row numbers
    tmp$Fraction <- tmp$Fraction / max(tmp$Fraction) # transform row numbers into fractions
    return(tmp)
    }) %>% bind_rows
  
  g <- ggplot(df, aes(x = log2(Intensity), y = Fraction, group = Channel)) + geom_point(aes(color = Channel), size = 0.1) +  # plot the ordered intensities against the fractions
    theme_bw() + guides(colour = guide_legend(override.aes = list(size = 5))) + 
    ggtitle("Intensity distribution plot") + xlab("log2 Intenisty") + ylab("Fraction")
  return(g)
})
 
#### Peptide Count ####
length <- map(.x = input.mq$Channel %>% as.factor() %>% levels(), .f = \(ee){
  df <- input.mq %>% subset(Channel == ee) %>% na.omit() %>% nrow() %>% as.data.frame() # count peptides (rows) in the raw data for each channel
  df$Channel <- ee # specify which channel was used
  colnames(df) <- c("Length", "Channel") # change column names
  return(df)
}) %>% bind_rows()

peptide_count <- ggplot(length, aes(x = Channel, y = Length)) + geom_bar(stat = "identity", width = 0.5, color = "black") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + theme(aspect.ratio = 4/3) + 
  ggtitle("Peptide Count") + xlab("Channel") + ylab("Nr of peptides") # plot counted peptides by channel

#### PCA plots ####
proc <- quant$FeatureLevelData # select processed, pre-comparison data
proc <- proc %>% distinct(Channel, PeptideSequence, .keep_all = TRUE) # retain only one charge state for each peptide
proc <- proc %>% mutate(log2Intensity = replace(log2Intensity, is.na(log2Intensity), 0)) # replace NAs in intensity column with 0
proc$Peptide_ID <- paste("p_", as.numeric(as.factor(proc$PSM)), sep = "") # create IDs for each peptides
proc_wide <- proc %>% select(Peptide_ID, Channel, log2Intensity) %>% # transform data to wide format
  arrange(Channel) %>% pivot_wider(names_from = Channel, values_from = log2Intensity)
pca_matrix <- proc_wide %>% column_to_rownames("Peptide_ID") %>% t # transpose the matrix and add channels as row names
pca <- prcomp(pca_matrix, center = T, scale = F) # calculate principal components

palette <- brewer.pal(length(levels(proc$Condition)),"Set1")[-6] # choose colours for a PCA biplot
groups <- proc %>% select(Channel, Condition) %>% distinct # extract the groups

pca_plot <- ggplot(pca, aes(x = PC1, y = PC2, colour = groups$Condition)) + geom_point() + theme_bw() + scale_color_manual(values = palette) + 
  geom_mark_ellipse(data = pca, aes(colour = groups$Condition)) + geom_text_repel(label = groups$Channel)

#### Volcano plots ####
vol_full <- ggplot(comp_result, aes(x = log2FC, y = -log10(adj.pvalue), label = paste(`Gene Names (primary)`, `Gene Names (ordered locus)`, sep = " - "))) + 
geom_point(mapping = aes(color = Sign), size = 1, alpha = 0.5) + 
geom_text_repel(data = subset(comp_result, Sign == "Sign"), size = 2.5, color =  "darkcyan", box.padding = 0.5, point.padding = 0.2, 
                 force = 10, segment.size = 0.1, segment.color = "black", segment.alpha = 0.3, max.iter = 1000, max.overlaps = 50) + 
scale_color_manual(values = c("red", "gray"), name = element_blank()) +
geom_hline(yintercept = -log10(cutoff), linetype="dashed", size = 0.2 , alpha = 0.5) + theme_bw() + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(), legend.position = "none") +
xlab("Log2(FC)") + ylab("-Log10(q-value)") + facet_wrap(vars(Label)) + ggtitle("Volcano Plot")

#### Fold Change Correlation Plot ####
sign_prot <- subset(comp_result, Sign == "Sign")$Entry # extract all significant proteins from comparison result
sign_comp <- comp_result %>% filter(Entry %in% sign_prot) # extract all instances of all proteins that is ever significant

fc_plots <- map(.x = 1:round(nrow(sign_comp)/(10*nbr_conc), digits = 0), .f = \(ff){ # split the significant proteins into manageable chunks
  rows <- ff*10*nbr_conc
  g <- ggplot(sign_comp[(rows-nbr_conc*10+1):rows,], aes(x = as.numeric(Label), y = log2FC, color = Entry, group = Entry, shape = Sign, # for each chunk of significant proteins, plot fold change against concentration
                             label = paste(`Gene Names (primary)`, `Gene Names (ordered locus)`, sep = " - "))) + 
    geom_text_repel(size = 3, color =  "darkcyan", box.padding = 0.5, point.padding = 0.2, 
                    force = 10, segment.size = 0.1, segment.color = "black", segment.alpha = 0.3, max.iter = 1000, max.overlaps = 50) +
    geom_line() + geom_point(size = 5) + theme_bw()  + guides(color = "none") + 
    xlab("Concentration (1 = low, 2 = medium, 3 = high)") + ylab("log2FC") + ggtitle(paste("Fold change across concentrations - ", ff))
})

#### Save Data ####
write_tsv(comp_result, paste(save_dir, "/tables/pisa_msstats_", date, "_", metabolite, "_", organism,".tsv", sep = "")) # save comparison data in table format
write_tsv(quant$ProteinLevelData, paste(save_dir, "/tables/pisa_msstats_", date, "_", metabolite, "_", organism, "_protein_quant", ".tsv", sep = "")) # save processed data pre-comparison

if (nrow(subset(comp_result, Sign == "Sign")) != 0) { # save all plots as a pdf, disregarding the fold change correlation plot if there are no significant proteins
  ggsave(filename = paste(save_dir, "/plots/pisa_msstats_", date, "_", metabolite, "_", organism, ".pdf", sep = ""), 
         plot = marrangeGrob(c(list(pval_dist, qq_plot, peptide_count, pca_plot, vol_full), dist_plot, fc_plots), nrow=1, ncol=1), width = 15, height = 9)
} else{
  ggsave(filename = paste(save_dir, "/plots/pisa_msstats_", date, "_", metabolite, "_", organism, ".pdf", sep = ""), 
         plot = marrangeGrob(c(list(pval_dist, qq_plot, peptide_count, pca_plot, vol_full), dist_plot), nrow=1, ncol=1), width = 15, height = 9)
}
