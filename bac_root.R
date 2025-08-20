#Analysis of root bacterial microbiome - rarefied----
library(tidyverse)
library(RColorBrewer)
library(vegan)
library(ggsci)

grDevices::windows()
#prep----
metadata <- as.data.frame(read_tsv('metadata_bagen.txt', lazy = FALSE))
tot_tax <- read_tsv('tot_tax_table.tsv', lazy = FALSE)

BaGen <- pull(metadata, seq_name)

roots <- filter(metadata, seq_name%in%BaGen & Type=="Root")%>%
  pull(seq_name)
soils <- filter(metadata, seq_name%in%BaGen & Type=="Soil"& Cultivar != "ntc")%>%
  pull(seq_name)

bac <- filter(tot_tax, Kingdom == "Bacteria")%>%
  pull(ASV)
pos <- read_tsv("banan_position.txt", lazy=FALSE)%>%
  right_join(metadata, by="Cultivar")

theme_dbox <- function(){
  theme_bw() %+replace%
    theme(
      axis.text.x = element_text(size = 12, color="black"),
      axis.text.y = element_text(size = 12, color="black"),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(angle = 90,size = 14),
      title = element_text(size = 14),
      legend.text = element_text(size = 12), 
      legend.position = 'right',
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
      strip.background = element_blank(),
      strip.text = element_text(color = "black", face = "bold", size = 12)
    )
}

pal_ramp <- function(values) {
  force(values)
  function(n) {
    if (n <= length(values)) {
      values[seq_len(n)]
    } else {
      colorRampPalette(values, alpha = TRUE)(n)
    }
  }
}
pal_adaptive <- function(name, palette, alpha = 1) {
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")
  
  raw_cols <- ggsci:::ggsci_db[[name]][[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  
  pal_ramp(unname(alpha_cols))
}
scale_color_adaptive <- function(name, palette, alpha = 1, ...) {
  ggplot2::discrete_scale("colour", name, pal_adaptive(name, palette, alpha), ...)
}

scale_fill_adaptive <- function(name, palette, alpha = 1, ...) {
  ggplot2::discrete_scale("fill", name, pal_adaptive(name, palette, alpha), ...)
}

#prepare data table----
#_Loading and removing lowest abundant ASVs----
ASV_table <- read_tsv('BaGen_ASV_count.tsv',lazy=FALSE)%>%
  filter(seq_name%in%roots)%>%
  select(seq_name,any_of(bac))%>%
  tibble::column_to_rownames("seq_name")

ASV_table <- ASV_table[,colSums(ASV_table)>0]
#removing the 10% least abundant ASVs

# Calculate relative abundance for each taxa (as a ratio per sample)
relative_abundance <- ASV_table / rowSums(ASV_table)

# Calculate the mean relative abundance for each taxa
mean_relative_abundance <- colMeans(relative_abundance, na.rm = TRUE)

# Determine the threshold for the 10% least abundant taxa
threshold <- quantile(mean_relative_abundance, 0.05)

# Identify columns to keep (those with abundance above the threshold)
columns_to_keep <- names(mean_relative_abundance[mean_relative_abundance > threshold])

# Subset the data frame to keep only the selected columns
f_ASV_table <- ASV_table[, columns_to_keep]

#keep taxa that are present in over 10% of samples

# Calculate the count of non-zero values for each taxa
bac_count_nonzero <- colSums(f_ASV_table > 0)

# Determine the threshold for taxa count (10% of the number of samples)
count_threshold <- 4

# Identify columns to keep (those with count above the threshold)
columns_to_keep_count <- names(bac_count_nonzero[bac_count_nonzero > count_threshold])

# Subset the data frame to keep only the selected columns
root_ASV_cts <- f_ASV_table[, columns_to_keep_count]
rootf_ASV_cts <- root_ASV_cts[rowSums(root_ASV_cts)>3000,]

#_subsampling----
num_sub <- min(rowSums(rootf_ASV_cts))
num_sub
set.seed(11092021)
# Function to perform subsampling
subsample_otu <- function(otu_table, reads_per_sample) {
  # Initialize a matrix to store subsampled counts, same dimensions as otu_table
  subsampled_matrix <- matrix(0, nrow = nrow(otu_table), ncol = ncol(otu_table))
  rownames(subsampled_matrix) <- rownames(otu_table)
  colnames(subsampled_matrix) <- colnames(otu_table)
  
  # Apply subsampling for each sample (column)
  for (sample_idx in 1:ncol(otu_table)) {
    sample_counts <- otu_table[, sample_idx]
    
    # Filter out zero values to avoid subsampling from OTUs with zero reads
    non_zero_counts <- sample_counts[sample_counts > 0]
    
    # If there are non-zero OTUs and enough reads, perform proportional subsampling
    if (length(non_zero_counts) > 0) {
      # Proportional sampling of reads_per_sample based on relative abundance
      otu_rep <- rep(rownames(otu_table)[sample_counts > 0], non_zero_counts)
      sampled_reads <- sample(otu_rep, size = reads_per_sample, replace = TRUE)
      
      # Convert the sampled reads to a table and update the matrix
      subsampled_counts <- table(sampled_reads)
      subsampled_matrix[names(subsampled_counts), sample_idx] <- as.numeric(subsampled_counts)
    }
  }
  
  return(subsampled_matrix)
}

# Parameters
reads_per_sample <- num_sub  # Number of reads to subsample from each sample
iterations <- 5000        # Number of subsampling iterations

input <- t(rootf_ASV_cts)
# Perform repeated subsampling
subsampled_results <- replicate(iterations, subsample_otu(input, reads_per_sample), simplify = FALSE)

# Calculate the mean number of reads for each OTU in each sample across iterations
mean_subsampled_table <- Reduce(`+`, subsampled_results) / iterations

# Convert result to a data frame for easier handling
mean_subsampled_df <- as.data.frame(t(mean_subsampled_table))
rowSums(mean_subsampled_df)
colSums(mean_subsampled_df)

mean_subsampled_df%>%
  tibble::rownames_to_column("seq_name")%>%
  write_tsv("rare_bacRoot_cts.tsv")

#_clr transformation----
mean_subsampled_df <- read_tsv("rare_bacRoot_cts.tsv", lazy=FALSE)%>%
  tibble::column_to_rownames("seq_name")
checkNumZerosCol <- apply(mean_subsampled_df,2,function(x) sum(x==0))
cases <- which(checkNumZerosCol == (nrow(mean_subsampled_df) - 1))
length(cases)
#cts <-BaGen_ASV_count[,-cases]
cts.no0<- zCompositions::cmultRepl(mean_subsampled_df, output = "p-counts",z.warning = 0.9)
Root_clr <- as.data.frame(compositions::clr(cts.no0))
Root_clr %>%
  tibble::rownames_to_column("seq_name")%>%
  write_tsv("rare_bacRoot_clr.tsv")

#PERMANOVA----
Root_clr <- read_tsv("rare_bacRoot_clr.tsv", lazy = FALSE)%>%
  tibble::column_to_rownames("seq_name")

sel_metadata <- filter(metadata, seq_name %in% unlist(rownames(Root_clr)))
sel_metadata <- sel_metadata[match(rownames(Root_clr), sel_metadata$seq_name), ]

adonis2(formula = Root_clr ~ Ext.Date, data = sel_metadata, method = "euclidean")

#Remove Ext.Date bias----
#__find ext.Date sig ASVs----
library(Maaslin2)
Root_clr <- as.data.frame(read_tsv("rare_bacRoot_clr.tsv", lazy = FALSE))
rownames(Root_clr)<-Root_clr$seq_name
names(Root_clr)[1] <- "Samples"

names <- pull(Root_clr, Samples)
input_metadata <- filter(metadata, seq_name%in%names)%>%
  select(seq_name, Genome,Cultivar, Ext.Date)
row.names(input_metadata) <- input_metadata$seq_name
names(input_metadata)[1] <- "Samples"

fit_data = Maaslin2(input_data     = Root_clr, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Ext.Date_13-9_Cult", 
                    fixed_effects  = c("Ext.Date"),
                    random_effects = c("Cultivar"),
                    reference      = c("Ext.Date,13/09/2021"))

fit_data = Maaslin2(input_data     = Root_clr, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Root_Ext.Date_04-10_Cult", 
                    fixed_effects  = c("Ext.Date"),
                    random_effects = c("Cultivar"),
                    reference      = c("Ext.Date,04/10/2021"))

fit_data = Maaslin2(input_data     = Root_clr, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Root_Ext.Date_13-01_Cult", 
                    fixed_effects  = c("Ext.Date"),
                    random_effects = c("Cultivar"),
                    reference      = c("Ext.Date,13/01/2023"))

#_remove ASVs----
ASV_table <- as.data.frame(read_tsv('BaGen_ASV_count.tsv',lazy=FALSE))%>%
    filter(seq_name %in% roots)%>%
  select(seq_name,any_of(bac))
  
Ext.ASVs1 <- read_tsv('Root_Ext.Date_04-10_Cult/significant_results.tsv')%>%
  filter(qval<0.2)%>%
  pull(feature)%>%
  unique()

Ext.ASVs2 <- read_tsv('Ext.Date_13-9_Cult/significant_results.tsv')%>%
  filter(qval<0.2)%>%
  pull(feature)%>%
  unique()

Ext.ASVs3 <- read_tsv('Root_Ext.Date_13-01_Cult/significant_results.tsv')%>%
  filter(qval<0.2)%>%
  pull(feature)%>%
  unique()

Ext.ASVs <- c(Ext.ASVs1, Ext.ASVs2,Ext.ASVs3)%>%
  unique()

f_ASV_table <- select(ASV_table,-any_of(Ext.ASVs))%>%
  tibble::column_to_rownames("seq_name")

#_removing lowest abundant ASVs----

f_ASV_table <- f_ASV_table[,colSums(f_ASV_table)>0]
#removing the 10% least abundant ASVs

# Calculate relative abundance for each taxa (as a ratio per sample)
relative_abundance <- f_ASV_table / rowSums(f_ASV_table)

# Calculate the mean relative abundance for each taxa
mean_relative_abundance <- colMeans(relative_abundance, na.rm = TRUE)

# Determine the threshold for the 10% least abundant taxa
threshold <- quantile(mean_relative_abundance, 0.05)

# Identify columns to keep (those with abundance above the threshold)
columns_to_keep <- names(mean_relative_abundance[mean_relative_abundance > threshold])

# Subset the data frame to keep only the selected columns
ff_ASV_table <- f_ASV_table[, columns_to_keep]

#keep taxa that are present in over 10% of samples

# Calculate the count of non-zero values for each taxa
bac_count_nonzero <- colSums(ff_ASV_table > 0)

# Determine the threshold for taxa count (10% of the number of samples)
count_threshold <- 4

# Identify columns to keep (those with count above the threshold)
columns_to_keep_count <- names(bac_count_nonzero[bac_count_nonzero > count_threshold])

# Subset the data frame to keep only the selected columns
root_ASV_cts <- ff_ASV_table[, columns_to_keep_count]
rootf_ASV_cts <- root_ASV_cts[rowSums(root_ASV_cts)>3000,]
rootf_ASV_cts%>%
  tibble::rownames_to_column("seq_name")%>%
  write_tsv("Root_bacASV_Ects.tsv")

#_subsampling----
num_sub <- min(rowSums(rootf_ASV_cts))
num_sub
set.seed(11092021)
# Function to perform subsampling
subsample_otu <- function(otu_table, reads_per_sample) {
  # Initialize a matrix to store subsampled counts, same dimensions as otu_table
  subsampled_matrix <- matrix(0, nrow = nrow(otu_table), ncol = ncol(otu_table))
  rownames(subsampled_matrix) <- rownames(otu_table)
  colnames(subsampled_matrix) <- colnames(otu_table)
  
  # Apply subsampling for each sample (column)
  for (sample_idx in 1:ncol(otu_table)) {
    sample_counts <- otu_table[, sample_idx]
    
    # Filter out zero values to avoid subsampling from OTUs with zero reads
    non_zero_counts <- sample_counts[sample_counts > 0]
    
    # If there are non-zero OTUs and enough reads, perform proportional subsampling
    if (length(non_zero_counts) > 0) {
      # Proportional sampling of reads_per_sample based on relative abundance
      otu_rep <- rep(rownames(otu_table)[sample_counts > 0], non_zero_counts)
      sampled_reads <- sample(otu_rep, size = reads_per_sample, replace = TRUE)
      
      # Convert the sampled reads to a table and update the matrix
      subsampled_counts <- table(sampled_reads)
      subsampled_matrix[names(subsampled_counts), sample_idx] <- as.numeric(subsampled_counts)
    }
  }
  
  return(subsampled_matrix)
}

# Parameters
reads_per_sample <- num_sub  # Number of reads to subsample from each sample
iterations <- 1000        # Number of subsampling iterations

input <- t(rootf_ASV_cts)
# Perform repeated subsampling
subsampled_results <- replicate(iterations, subsample_otu(input, reads_per_sample), simplify = FALSE)

# Calculate the mean number of reads for each OTU in each sample across iterations
mean_subsampled_table <- Reduce(`+`, subsampled_results) / iterations

# Convert result to a data frame for easier handling
mean_subsampled_df <- as.data.frame(t(mean_subsampled_table))
rowSums(mean_subsampled_df)
colSums(mean_subsampled_df)

mean_subsampled_df%>%
  tibble::rownames_to_column("seq_name")%>%
  write_tsv("rare_bacRoot_Ects.tsv")

#_clr transformation----
mean_subsampled_df <- read_tsv("rare_bacRoot_Ects.tsv", lazy=FALSE)%>%
  tibble::column_to_rownames("seq_name")
checkNumZerosCol <- apply(mean_subsampled_df,2,function(x) sum(x==0))
cases <- which(checkNumZerosCol == (nrow(mean_subsampled_df) - 1))
length(cases)
#cts <-BaGen_ASV_count[,-cases]
cts.no0<- zCompositions::cmultRepl(mean_subsampled_df, output = "p-counts",z.warning = 0.9)
Root_clr <- as.data.frame(compositions::clr(cts.no0))
Root_clr %>%
  tibble::rownames_to_column("seq_name")%>%
  write_tsv("rare_bacRoot_Eclr.tsv")

#PERMANOVA filtered table----
Root_clr <- read_tsv("rare_bacRoot_Eclr.tsv", lazy = FALSE)%>%
  tibble::column_to_rownames("seq_name")

sel_metadata <- filter(metadata, seq_name %in% unlist(rownames(Root_clr)))
sel_metadata <- sel_metadata[match(rownames(Root_clr), sel_metadata$seq_name), ]

adonis2(formula = Root_clr ~ Ext.Date, data = sel_metadata, method = "euclidean")

#no ext date
adonis2(formula = Root_clr ~ Cultivar, data = sel_metadata, method = "euclidean",by="term")
adonis2(formula = Root_clr ~ Genome+Cultivar, data = sel_metadata, method = "euclidean",by="term")
adonis2(formula = Root_clr ~ Genome, data = sel_metadata, method = "euclidean",by="term")
adonis2(formula = Root_clr ~ `A/B`+Genome+Cultivar, data = sel_metadata, method = "euclidean",by="term")
adonis2(formula = Root_clr ~ Ploidity+`A/B`+Genome+Cultivar, data = sel_metadata, method = "euclidean",by="term")
adonis2(formula = Root_clr ~ Ploidity+`A/B`+Genome, data = sel_metadata, method = "euclidean",by="term")
adonis2(formula = Root_clr ~ Ploidity*`A/B`+Cultivar, data = sel_metadata, method = "euclidean",by="term")

#chck plot position
sel_pos <- filter(pos, seq_name %in% unlist(rownames(Root_clr)))
sel_pos <- sel_pos[match(rownames(Root_clr), sel_pos$seq_name), ]
adonis2(formula = Root_clr ~ plot_Column, data = sel_pos, method = "euclidean")
adonis2(formula = Root_clr ~ plot_Row, data = sel_pos, method = "euclidean")
adonis2(formula = Root_clr ~ plot_Row+plot_Column, data = sel_pos, method = "euclidean", by="term")

#PCA filtered----
Root_clr <- read_tsv("rare_bacRoot_Eclr.tsv", lazy = FALSE)%>%
  tibble::column_to_rownames("seq_name")
pca.clr <- prcomp(Root_clr, scale. = FALSE)
summary(pca.clr)
loadings <- as.data.frame(pca.clr$rotation)
PCA <- as.data.frame(pca.clr[["x"]])
PCA$seq_name <- unlist(row.names(PCA))
PCA <- inner_join(PCA, metadata)
write_tsv(PCA, 'ERoot_clr_PCA.tsv')
loadings%>%
  tibble::rownames_to_column("ASV")%>%
  write_tsv('ERoot_clr_PCA_load.tsv')

PCA <- read_tsv('ERoot_clr_PCA.tsv', lazy = FALSE)
PCA$Ploidity <- as.factor(PCA$Ploidity)
PCA$Column <- as.factor(PCA$Column)
PCA$Plate <- as.factor(PCA$Plate)


getPalette = colorRampPalette(brewer.pal(9, "Set1"))

colourCount = length(unique(PCA$Cultivar))
PCA%>%
  ggplot(aes(x=PC1, y=PC2, color=Cultivar)) +
  geom_point(aes(color = Cultivar, shape = `A/B`), size = 4, na.rm = TRUE) +
  scale_shape_discrete(name='Presence of B-type Genome')+
  #scale_color_brewer(palette = 'Set2') +
  scale_color_manual(values = getPalette(colourCount))+
  coord_fixed()+#(xlim = c(-50,50),ylim = c(-50,50)) +
  #facet_grid(time ~ soil) +
  #facet_grid(time ~ .) +
  labs(title="Bacterial communities of Banana Roots",
       x="PC 1 (8% variance)",
       y="PC 2 (6% variance)") +
  theme_dbox()+
  geom_hline(yintercept = 0, color = "black")+
  geom_vline(xintercept = 0, color = "black")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA, linewidth=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()#,
        #        text = element_text(size=20)
  )
ggsave("ERoot_PCA12_treat.png", height = 180, width = 180, units = "mm", dpi=300)

PCA%>%
  ggplot(aes(x=PC2, y=PC3, color=Cultivar)) +
  geom_point(aes(color = Cultivar, shape = `A/B`), size = 4, na.rm = TRUE) +
  scale_shape_discrete(name='Presence of B-type Genome')+
  #scale_color_brewer(palette = 'Set2') +
  scale_color_manual(values = getPalette(colourCount))+
  coord_fixed()+#(xlim = c(-50,50),ylim = c(-50,50)) +
  #facet_grid(time ~ soil) +
  #facet_grid(time ~ .) +
  labs(title="Bacterial communities of Banana Roots",
       x="PC 2 (6% variance)",
       y="PC 3 (5% variance)") +
  theme_dbox()+
  geom_hline(yintercept = 0, color = "black")+
  geom_vline(xintercept = 0, color = "black")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA, linewidth=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()#,
        #        text = element_text(size=20)
  )
ggsave("ERoot_PCA23_treat.png", height = 180, width = 180, units = "mm", dpi=300)

PCA%>%
  ggplot(aes(x=PC1, y=PC2, color=Genome)) +
  geom_point(aes(color = Genome, shape = `A/B`), size = 4, na.rm = TRUE) +
  scale_shape_discrete(name='Subgenome')+
  scale_color_jco()+
  #scale_color_manual(values = getPalette(colourCount))+
  coord_fixed()+#(xlim = c(-50,50),ylim = c(-50,50)) +
  #facet_grid(time ~ soil) +
  #facet_grid(time ~ .) +
  labs(title=NULL,
       x="PC 1 (8% variance)",
       y="PC 2 (6% variance)") +
  theme_dbox()+
  geom_hline(yintercept = 0, color = "black")+
  geom_vline(xintercept = 0, color = "black")+
  theme(panel.background = element_blank(),
        legend.position = "bottom",
  #      panel.border = element_rect(colour = "black", 
  #                                  fill=NA, linewidth=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()#,
        #        text = element_text(size=20)
  )
ggsave("ERoot_PCA12_AvB_genome.png", height = 180, width = 180, units = "mm", dpi=300)

#cultivar DAs----
library(Maaslin2)
ERoot_clr_bac <- as.data.frame(read_tsv("rare_bacRoot_Eclr.tsv", lazy = FALSE))
rownames(ERoot_clr_bac) <- ERoot_clr_bac$seq_name

names(ERoot_clr_bac)[1] <- "Samples"
names <- pull(ERoot_clr_bac, Samples)

input_metadata <- filter(metadata, seq_name%in%names)%>%
  select(seq_name, Genome,Cultivar, Ext.Date)%>%
  as.data.frame()

row.names(input_metadata) <- input_metadata$seq_name
names(input_metadata)[1] <- "Samples"

fit_data = Maaslin2(input_data     = ERoot_clr_bac, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    standardize = FALSE,
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Root_Treat_blue.E", 
                    fixed_effects  = c("Cultivar"),
                    random_effects = c("Ext.Date"),
                    reference      = c("Cultivar,Blue-banana"))

fit_data = Maaslin2(input_data     = ERoot_clr_bac, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    standardize = FALSE,
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Root_Treat_blugg.E", 
                    fixed_effects  = c("Cultivar"),
                    random_effects = c("Ext.Date"),
                    reference      = c("Cultivar,Bluggoe"))

fit_data = Maaslin2(input_data     = ERoot_clr_bac, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    standardize = FALSE,
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Root_Treat_lak.E", 
                    fixed_effects  = c("Cultivar"),
                    random_effects = c("Ext.Date"),
                    reference      = c("Cultivar,Lakatan"))

fit_data = Maaslin2(input_data     = ERoot_clr_bac, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    standardize = FALSE,
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Root_Calc4.E", 
                    fixed_effects  = c("Cultivar"),
                    random_effects = c("Ext.Date"),
                    reference      = c("Cultivar,Calcutta4"))

fit_data = Maaslin2(input_data     = ERoot_clr_bac, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    standardize = FALSE,
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Root_arab.E", 
                    fixed_effects  = c("Cultivar"),
                    random_effects = c("Ext.Date"),
                    reference      = c("Cultivar,Pissang-Awak"))

fit_data = Maaslin2(input_data     = ERoot_clr_bac, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    standardize = FALSE,
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Root_grand.E", 
                    fixed_effects  = c("Cultivar"),
                    random_effects = c("Ext.Date"),
                    reference      = c("Cultivar,Grand-Naine-5/37"))

fit_data = Maaslin2(input_data     = ERoot_clr_bac, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    standardize = FALSE,
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Root_rose.E", 
                    fixed_effects  = c("Cultivar"),
                    random_effects = c("Ext.Date"),
                    reference      = c("Cultivar,Aacv-Rose"))

fit_data = Maaslin2(input_data     = ERoot_clr_bac, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    standardize = FALSE,
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "malac.E", 
                    fixed_effects  = c("Cultivar"),
                    random_effects = c("Ext.Date"),
                    reference      = c("Cultivar,M.Malaccensis"))

fit_data = Maaslin2(input_data     = ERoot_clr_bac, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    standardize = FALSE,
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Root_balb.E", 
                    fixed_effects  = c("Cultivar"),
                    random_effects = c("Ext.Date"),
                    reference      = c("Cultivar,Balbisiana"))

fit_data = Maaslin2(input_data     = ERoot_clr_bac, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    standardize = FALSE,
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Root_natan.E", 
                    fixed_effects  = c("Cultivar"),
                    random_effects = c("Ext.Date"),
                    reference      = c("Cultivar,Natan"))

fit_data = Maaslin2(input_data     = ERoot_clr_bac, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    standardize = FALSE,
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Root_apple.E", 
                    fixed_effects  = c("Cultivar"),
                    random_effects = c("Ext.Date"),
                    reference      = c("Cultivar,Apple-banana"))

fit_data = Maaslin2(input_data     = ERoot_clr_bac, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    standardize = FALSE,
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Root_Hom.E", 
                    fixed_effects  = c("Cultivar"),
                    random_effects = c("Ext.Date"),
                    reference      = c("Cultivar,Hom-Thong-Mokho"))

fit_data = Maaslin2(input_data     = ERoot_clr_bac, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    standardize = FALSE,
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Root_Plantain.E", 
                    fixed_effects  = c("Cultivar"),
                    random_effects = c("Ext.Date"),
                    reference      = c("Cultivar,Plantain"))

fit_data = Maaslin2(input_data     = ERoot_clr_bac, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    standardize = FALSE,
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Root_Prata.E", 
                    fixed_effects  = c("Cultivar"),
                    random_effects = c("Ext.Date"),
                    reference      = c("Cultivar,Prata-Ana"))

#_list Cultivar DA bacteria----
DA_arab <- read_tsv('Root_arab.E/significant_results.tsv', lazy = FALSE)%>%
  filter(qval <0.05)%>%
  mutate(ref="Pissang-Awak")
DA_calc4 <- read_tsv('Root_Calc4.E/significant_results.tsv', lazy = FALSE)%>%
  filter(qval <0.05)%>%
  mutate(ref="Calcutta4")
DA_grand <- read_tsv('Root_grand.E/significant_results.tsv', lazy = FALSE)%>%
  filter(qval <0.05)%>%
  mutate(ref="Grand-Naine-5/37")
DA_rose <- read_tsv('Root_rose.E/significant_results.tsv', lazy = FALSE)%>%
  filter(qval <0.05)%>%
  mutate(ref="Aacv-Rose")
DA_blue <- read_tsv('Root_Treat_blue.E/significant_results.tsv', lazy = FALSE)%>%
  filter(qval <0.05)%>%
  mutate(ref="Blue-banana")
DA_malac <- read_tsv('malac.E/significant_results.tsv', lazy = FALSE)%>%
  filter(qval <0.05)%>%
  mutate(ref="M.Malaccensis")
DA_balb <- read_tsv('Root_balb.E/significant_results.tsv', lazy = FALSE)%>%
  filter(qval <0.05)%>%
  mutate(ref="Balbisiana")
DA_natan <- read_tsv('Root_natan.E/significant_results.tsv', lazy = FALSE)%>%
  filter(qval <0.05)%>%
  mutate(ref="Natan")
DA_apple <- read_tsv('Root_apple.E/significant_results.tsv', lazy = FALSE)%>%
  filter(qval <0.05)%>%
  mutate(ref="Apple-banana")
DA_bluggoe <- read_tsv('Root_Treat_blugg.E/significant_results.tsv', lazy = FALSE)%>%
  filter(qval <0.05)%>%
  mutate(ref="Bluggoe")
DA_hom <- read_tsv('Root_Hom.E/significant_results.tsv', lazy = FALSE)%>%
  filter(qval <0.05)%>%
  mutate(ref="Hom-Thong-Mokho")
DA_plant <- read_tsv('Root_Plantain.E/significant_results.tsv', lazy = FALSE)%>%
  filter(qval <0.05)%>%
  mutate(ref="Plantain")
DA_prata <- read_tsv('Root_Prata.E/significant_results.tsv', lazy = FALSE)%>%
  filter(qval <0.05)%>%
  mutate(ref="Prata-Ana")
DA_lak <- read_tsv('Root_Treat_lak.E/significant_results.tsv', lazy = FALSE)%>%
  filter(qval <0.05)%>%
  mutate(ref="Lakatan")


Cultivar_DA <- rbind(DA_arab,DA_calc4)%>%
  rbind(DA_grand)%>%
  rbind(DA_rose)%>%
  rbind(DA_blue)%>%
  rbind(DA_lak)%>%
  rbind(DA_malac)%>%
  rbind(DA_balb)%>%
  rbind(DA_natan)%>%
  rbind(DA_apple)%>%
  rbind(DA_prata)%>%
  rbind(DA_plant)%>%
  rbind(DA_hom)%>%
  rbind(DA_bluggoe)%>%
  inner_join(tot_tax, by=c("feature"="ASV"))%>%
  write_tsv('cultivar_DA_RootE.tsv')

length(unique(Cultivar_DA$feature))

#_plot DAs----
Cultivar_DA <- read_tsv('cultivar_DA_RootE.tsv', lazy=FALSE)

refcult_subgenome <- select(metadata, Cultivar, `A/B`)%>%
  unique()
colnames(refcult_subgenome)<- c("ref","ref_subgenome")

cult_DA_subgen <- select(Cultivar_DA,-c(metadata,stderr,pval,qval, N,stderr,`N.not.0` ))%>%
  filter(coef >0)%>%
  select(-c(coef,ref))%>%
  unique()%>%
  left_join(select(metadata, Cultivar, `A/B`), by=c("value"="Cultivar"))%>%
  unique()%>%
  write_tsv('cult_DA_subgenome.tsv')

cult_DA_subgen%>%
  ggplot(aes( x=value, fill=Phylum))+
  geom_bar()+
  scale_fill_adaptive(name = "jama", palette = "default")+
  labs(fill="Phylum (bacteria)")+
  facet_grid(.~`A/B`,scales = "free")+
  theme_dbox()+
  labs(#color = 'Sample',size='alr',
    x="Banana cultivar",
    y="Positively associated ASV")+
  theme(legend.position = "right",
        legend.margin = margin(0, 0, 0, 0),
        legend.key.spacing.y = unit(0.005, "cm"),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1),
        axis.text.x =  element_text(size = 12,angle =45),
        strip.text.y = element_text(size=12, angle = 0, lineheight=1),
        strip.text.x = element_text(size=12))
ggsave("DA_BAR_P.png", width = 180, height = 180, units = "mm", dpi=300)

cult_DA_subgen%>%
  ggplot(aes( x=value, fill=Class))+
  geom_bar()+
  scale_fill_adaptive(name = "jama", palette = "default")+
  labs(fill="Class (bacteria)")+
  facet_grid(.~`A/B`,scales = "free")+
  theme_dbox()+
  labs(#color = 'Sample',size='alr',
    x="Banana cultivar",
    y="Positively associated ASV")+
  theme(legend.position = "right",
        legend.margin = margin(0, 0, 0, 0),
        legend.key.spacing.y = unit(0.005, "cm"),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1),
        #axis.text.x = element_blank(),
        axis.text.x =  element_text(size = 12,angle =45),
        strip.text.y = element_text(size=12, angle = 0, lineheight=1),
        strip.text.x = element_text(size=12))
ggsave("DA_BAR_C.png", width = 180, height = 180, units = "mm", dpi=300)
