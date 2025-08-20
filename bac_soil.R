#bac soil banana cultivar analysis
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
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(angle = 90,size = 14),
      title = element_text(size = 14),
      legend.text = element_text(size = 12), 
      legend.position = 'right',
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1)  
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
  filter(seq_name%in%soils)%>%
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
soil_ASV_cts <- f_ASV_table[, columns_to_keep_count]
soilf_ASV_cts <- soil_ASV_cts[rowSums(soil_ASV_cts)>2000,]

#_subsampling----
num_sub <- min(rowSums(soilf_ASV_cts))
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

input <- t(soilf_ASV_cts)
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
  write_tsv("rare_bacSoil_cts.tsv")

#_clr transformation----
mean_subsampled_df <- read_tsv("rare_bacSoil_cts.tsv", lazy=FALSE)%>%
  tibble::column_to_rownames("seq_name")
checkNumZerosCol <- apply(mean_subsampled_df,2,function(x) sum(x==0))
cases <- which(checkNumZerosCol == (nrow(mean_subsampled_df) - 1))
length(cases)
#cts <-BaGen_ASV_count[,-cases]
cts.no0<- zCompositions::cmultRepl(mean_subsampled_df, output = "p-counts",z.warning = 0.9)
Soil_clr <- as.data.frame(compositions::clr(cts.no0))
Soil_clr %>%
  tibble::rownames_to_column("seq_name")%>%
  write_tsv("rare_bacSoil_clr.tsv")

#PERMANOVA----
Soil_clr <- read_tsv("rare_bacSoil_clr.tsv", lazy = FALSE)%>%
  tibble::column_to_rownames("seq_name")

sel_metadata <- filter(metadata, seq_name %in% unlist(rownames(Soil_clr)))
sel_metadata <- sel_metadata[match(rownames(Soil_clr), sel_metadata$seq_name), ]

adonis2(formula = Soil_clr ~ Ext.Date, data = sel_metadata, method = "euclidean")

#Remove Ext.Date bias----
#__find ext.Date sig ASVs----
library(Maaslin2)
Soil_clr <- as.data.frame(read_tsv("rare_bacSoil_clr.tsv", lazy = FALSE))
rownames(Soil_clr)<-Soil_clr$seq_name
names(Soil_clr)[1] <- "Samples"

names <- pull(Soil_clr, Samples)
input_metadata <- filter(metadata, seq_name%in%names)%>%
  select(seq_name, Genome,Cultivar, Ext.Date)
row.names(input_metadata) <- input_metadata$seq_name
names(input_metadata)[1] <- "Samples"

fit_data = Maaslin2(input_data     = Soil_clr, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Ext.Date_28-3_Soil", 
                    fixed_effects  = c("Ext.Date"),
                    random_effects = c("Cultivar"),
                    reference      = c("Ext.Date,28/03/2023"))

fit_data = Maaslin2(input_data     = Soil_clr, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Root_Ext.Date_16-3_Cult", 
                    fixed_effects  = c("Ext.Date"),
                    random_effects = c("Cultivar"),
                    reference      = c("Ext.Date,16/03/2023"))

fit_data = Maaslin2(input_data     = Soil_clr, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "Root_Ext.Date_16-03_Cult", 
                    fixed_effects  = c("Ext.Date"),
                    random_effects = c("Cultivar"),
                    reference      = c("Ext.Date,06/03/2023"))

#_remove ASVs----
ASV_table <- as.data.frame(read_tsv('BaGen_ASV_count.tsv',lazy=FALSE))%>%
  filter(seq_name %in% soils)%>%
  select(seq_name,any_of(bac))

Ext.ASVs1 <- read_tsv('Root_Ext.Date_16-03_Cult/significant_results.tsv')%>%
  filter(qval<0.2)%>%
  pull(feature)%>%
  unique()

Ext.ASVs2 <- read_tsv('Root_Ext.Date_16-3_Cult/significant_results.tsv')%>%
  filter(qval<0.2)%>%
  pull(feature)%>%
  unique()

Ext.ASVs3 <- read_tsv('Ext.Date_28-3_Soil/significant_results.tsv')%>%
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
soil_ASV_cts <- ff_ASV_table[, columns_to_keep_count]
soilf_ASV_cts <- soil_ASV_cts[rowSums(soil_ASV_cts)>2000,]
soilf_ASV_cts%>%
  tibble::rownames_to_column("seq_name")%>%
  write_tsv("Soil_bacASV_Ects.tsv")

#_subsampling----
num_sub <- min(rowSums(soilf_ASV_cts))
num_sub
set.seed(11092021)

# Parameters
reads_per_sample <- num_sub  # Number of reads to subsample from each sample
iterations <- 5000        # Number of subsampling iterations

input <- t(soilf_ASV_cts)
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
  write_tsv("rare_bacSoil_Ects.tsv")

#_clr transformation----
mean_subsampled_df <- read_tsv("rare_bacSoil_Ects.tsv", lazy=FALSE)%>%
  tibble::column_to_rownames("seq_name")
checkNumZerosCol <- apply(mean_subsampled_df,2,function(x) sum(x==0))
cases <- which(checkNumZerosCol == (nrow(mean_subsampled_df) - 1))
length(cases)
#cts <-BaGen_ASV_count[,-cases]
cts.no0<- zCompositions::cmultRepl(mean_subsampled_df, output = "p-counts",z.warning = 0.9)
Soil_clr <- as.data.frame(compositions::clr(cts.no0))
Soil_clr %>%
  tibble::rownames_to_column("seq_name")%>%
  write_tsv("rare_bacSoil_Eclr.tsv")

#PERMANOVA filtered table----
Soil_clr <- read_tsv("rare_bacSoil_Eclr.tsv", lazy = FALSE)%>%
  tibble::column_to_rownames("seq_name")

sel_metadata <- filter(metadata, seq_name %in% unlist(rownames(Soil_clr)))
sel_metadata <- sel_metadata[match(rownames(Soil_clr), sel_metadata$seq_name), ]


adonis2(formula = Soil_clr ~ Ext.Date, data = sel_metadata, method = "euclidean")

#chck plot position
sel_pos <- filter(pos, seq_name %in% unlist(rownames(Soil_clr)))
sel_pos <- sel_pos[match(rownames(Soil_clr), sel_pos$seq_name), ]
adonis2(formula = Soil_clr ~ plot_Column, data = sel_pos, method = "euclidean")
adonis2(formula = Soil_clr ~ plot_Row, data = sel_pos, method = "euclidean")
adonis2(formula = Soil_clr ~ plot_Row+plot_Column, data = sel_pos, method = "euclidean",by="term")


#no ext date
adonis2(formula = Soil_clr ~ Cultivar, data = sel_metadata, method = "euclidean",by="term")
adonis2(formula = Soil_clr ~ Genome+Cultivar, data = sel_metadata, method = "euclidean",by="term")
adonis2(formula = Soil_clr ~ Genome, data = sel_metadata, method = "euclidean",by="term")
adonis2(formula = Soil_clr ~ `A/B`+Genome+Cultivar, data = sel_metadata, method = "euclidean",by="term")
adonis2(formula = Soil_clr ~ Ploidity+`A/B`+Genome+Cultivar, data = sel_metadata, method = "euclidean",by="term")
adonis2(formula = Soil_clr ~ Ploidity+`A/B`+Genome, data = sel_metadata, method = "euclidean",by="term")
adonis2(formula = Soil_clr ~ Ploidity*`A/B`+Cultivar, data = sel_metadata, method = "euclidean",by="term")

#PCA filtered----
Soil_clr <- read_tsv("rare_bacSoil_Eclr.tsv", lazy = FALSE)%>%
  tibble::column_to_rownames("seq_name")
pca.clr <- prcomp(Soil_clr, scale. = FALSE)
summary(pca.clr)
loadings <- as.data.frame(pca.clr$rotation)
PCA <- as.data.frame(pca.clr[["x"]])
PCA$seq_name <- unlist(row.names(PCA))
PCA <- inner_join(PCA, metadata)
write_tsv(PCA, 'ESoil_clr_PCA.tsv')
loadings%>%
  tibble::rownames_to_column("ASV")%>%
  write_tsv('ESoil_clr_PCA_load.tsv')

PCA <- read_tsv('ESoil_clr_PCA.tsv', lazy = FALSE)
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
  labs(title="Bacterial communities of Banana Soil",
       x="PC 1 (6% variance)",
       y="PC 2 (5% variance)") +
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
ggsave("ESoil_PCA12_treat.png", height = 180, width = 180, units = "mm", dpi=300)

PCA%>%
  ggplot(aes(x=PC2, y=PC3, color=Cultivar)) +
  geom_point(aes(color = Cultivar, shape = `A/B`), size = 4, na.rm = TRUE) +
  scale_shape_discrete(name='Presence of B-type Genome')+
  #scale_color_brewer(palette = 'Set2') +
  scale_color_manual(values = getPalette(colourCount))+
  coord_fixed()+#(xlim = c(-50,50),ylim = c(-50,50)) +
  #facet_grid(time ~ soil) +
  #facet_grid(time ~ .) +
  labs(title="Bacterial communities of Banana Soil",
       x="PC 2 (5% variance)",
       y="PC 3 (4% variance)") +
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
ggsave("ESoil_PCA23_treat.png", height = 180, width = 180, units = "mm", dpi=300)

PCA%>%
  ggplot(aes(x=PC1, y=PC2, color=Genome)) +
  geom_point(aes(color = Genome, shape = `A/B`), size = 4, na.rm = TRUE) +
  scale_shape_discrete(name='Type-B Genome presence')+
  scale_color_jco()+
  #scale_color_manual(values = getPalette(colourCount))+
  coord_fixed()+#(xlim = c(-50,50),ylim = c(-50,50)) +
  #facet_grid(time ~ soil) +
  #facet_grid(time ~ .) +
  labs(title="Bacterial communities of Banana Soil",
       x="PC 1 (5% variance)",
       y="PC 2 (5% variance)") +
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
ggsave("ESoil_PCA12_AvB_genome.png", height = 180, width = 180, units = "mm", dpi=300)

