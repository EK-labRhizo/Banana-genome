#Root network analysis
library(tidyverse)
library(RColorBrewer)
library(vegan)
library(ggsci)
library(NetCoMi)

grDevices::windows()

#prep----
metadata <- as.data.frame(read_tsv('metadata_bagen.txt', lazy = FALSE))
tot_tax <- read_tsv('tot_tax_table.tsv', lazy = FALSE)%>%
  mutate(Taxa=paste(Phylum,Class,Order, Family,Genus, sep="_"))
remove_fung <- filter(tot_tax, Phylum=="Fungi_phy_Incertae_sedis")%>%
  pull(ASV)

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
scale_fill_adaptive <- function(name, palette, alpha = 1, ...) {
  ggplot2::discrete_scale("fill", name, pal_adaptive(name, palette, alpha), ...)
}


#A/B network----
remove_fung <- filter(tot_tax, Phylum=="Fungi_phy_Incertae_sedis")%>%
  pull(ASV)
Root_Ects<- as.data.frame(read_tsv("Root_ASV_Ects.tsv",lazy = FALSE))%>%
  select(-any_of(remove_fung))%>%
  tibble::column_to_rownames("seq_name")
sample_order <- row.names(Root_Ects)
sample_types <- metadata$`A/B`[match(sample_order, metadata$seq_name)]
sample_types_f <- factor(sample_types)
levels(sample_types_f)

Root_Ects <- as.matrix(Root_Ects)
AvB_net_var <- netConstruct(Root_Ects,
                              group = sample_types,
                              dataType="counts",
                              filtTax = "highestVar",
                              filtTaxPar = list(highestVar = 1000),
                              filtSamp = "none",
                              measure = "spieceasi",
                              measurePar = list(nlambda=20, 
                                                pulsar.params=list(rep.num=50)),
                              normMethod = "none",
                              zeroMethod = "none",
                              sparsMethod = "none", 
                              dissFunc = "signed",
                              verbose = 3,
                              seed = 11092003)

AvB_props_var <- netAnalyze(AvB_net_var, 
                              centrLCC = FALSE,
                              clustMethod = "cluster_louvain",
                              hubPar = c("eigenvector"))
summary(AvB_props_var)

#_compare network----
comp_AvB_perm <- netCompare(AvB_props_var, 
                            permTest = TRUE, 
                            nPerm = 50,
                            verbose = FALSE,
                            adjust = "adaptBH",
                            seed = 11092003)

summary(comp_AvB_perm, 
        groupNames = c("A", "B"),
        showCentr = c("degree", "between", "closeness","eigenvector"), 
        numbNodes = 5)

summary_output <- capture.output(summary(comp_AvB_perm, 
                                         groupNames = c("A", "B"),
                                         showCentr = c("degree", "between", "closeness","eigenvector")))

output_file <- "AvB_network_summary.txt"
writeLines(summary_output, output_file)
cat("Summary saved to", output_file)

#_Extract parameters----
tot_net_tax <- as.data.frame(AvB_props_var[["lccNames1"]])%>%
  rename("ASV"=`AvB_props_var[["lccNames1"]]`)%>%
  left_join(tot_tax)%>%
  write_tsv("AvB_var_net_tax.tsv")
#A
louvain_root <- as.data.frame(AvB_props_var[["clustering"]][[("clust1")]])%>%
  tibble::rownames_to_column("ASV")%>%
  inner_join(tot_tax)%>%
  unique()%>%
  write_tsv("A_louvain_rootVar_ASV.tsv")

central_root <- as.data.frame(AvB_props_var[["centralities"]][["degree1"]])%>%
  tibble::rownames_to_column("ASV")%>%
  inner_join(tot_tax)%>%
  unique()%>%
  write_tsv("A_central_rootVar_ASV.tsv")

between_root <- as.data.frame(AvB_props_var[["centralities"]][["between1"]])%>%
  tibble::rownames_to_column("ASV")%>%
  inner_join(tot_tax)%>%
  unique()%>%
  write_tsv("A_between_rootVar_ASV.tsv")

close_root <- as.data.frame(AvB_props_var[["centralities"]][["close1"]])%>%
  tibble::rownames_to_column("ASV")%>%
  inner_join(tot_tax)%>%
  unique()%>%
  write_tsv("A_close_rootVar_ASV.tsv")


eigen_root <- as.data.frame(AvB_props_var[["centralities"]][["eigenv1"]])%>%
  tibble::rownames_to_column("ASV")%>%
  inner_join(tot_tax)%>%
  unique()%>%
  write_tsv("A_eigen_rootVar_ASV.tsv")

hubs_root <- as.data.frame(AvB_props_var[["hubs"]][["hubs1"]])%>%
  inner_join(tot_tax, by=c(`AvB_props_var[["hubs"]][["hubs1"]]`="ASV"))%>%
  write_tsv("A_hubs_rootVar_ASV.tsv")

#B
louvain_soil <- as.data.frame(AvB_props_var[["clustering"]][["clust2"]])%>%
  tibble::rownames_to_column("ASV")%>%
  inner_join(tot_tax)%>%
  unique()%>%
  write_tsv("B_louvain_rootVar_ASV.tsv")

central_soil <- as.data.frame(AvB_props_var[["centralities"]][["degree2"]])%>%
  tibble::rownames_to_column("ASV")%>%
  inner_join(tot_tax)%>%
  unique()%>%
  write_tsv("B_central_rootVar_ASV.tsv")

between_soil <- as.data.frame(AvB_props_var[["centralities"]][["between2"]])%>%
  tibble::rownames_to_column("ASV")%>%
  inner_join(tot_tax)%>%
  unique()%>%
  write_tsv("B_between_rootVar_ASV.tsv")

close_soil <- as.data.frame(AvB_props_var[["centralities"]][["close2"]])%>%
  tibble::rownames_to_column("ASV")%>%
  inner_join(tot_tax)%>%
  unique()%>%
  write_tsv("B_close_rootVar_ASV.tsv")


eigen_soil <- as.data.frame(AvB_props_var[["centralities"]][["eigenv2"]])%>%
  tibble::rownames_to_column("ASV")%>%
  inner_join(tot_tax)%>%
  unique()%>%
  write_tsv("B_eigen_rootVar_ASV.tsv")

hubs_root <- as.data.frame(AvB_props_var[["hubs"]][["hubs2"]])%>%
  inner_join(tot_tax, by=c(`AvB_props_var[["hubs"]][["hubs2"]]`="ASV"))%>%
  write_tsv("B_hubs_rootVar_ASV.tsv")

#_view----
plot(AvB_props_var, 
     repulsion = 1.7,
     nodeColor = "cluster",
     highlightHubs = FALSE,
     nodeSize = "clr",
     edgeTranspHigh = 20,
     hubBorderCol  = "gray40")

#_plot----
#__color by taxa----
#___phylum----
set.seed(11092003)
net <- plot(AvB_props_var, 
            sameLayout = TRUE, 
            layoutGroup = "union", 
            nodeFilter = "highestClose",
            nodeFilterPar = 25,
            nodeSize = "eigenvector",
            rmSingles = "inboth",
            groupNames = c("A", "B"),
            hubBorderCol  = "gray40")

phylum_mapping <- setNames(taxonomy_df1$Phylum, taxonomy_df1$ASV)
phylum_labels <- phylum_mapping[current_labels1]

class_mapping <- setNames(taxonomy_df1$Class, taxonomy_df1$ASV)
class_labels <- class_mapping[current_labels1]

Order_mapping <- setNames(taxonomy_df1$Order, taxonomy_df1$ASV)
Order_labels <- Order_mapping[current_labels1]

shape <- select(taxonomy_df1, ASV, Kingdom)%>%
  mutate(shape=ifelse(Kingdom=="Bacteria", 16, 17))

node_shape <- factor(setNames(shape$shape,shape$ASV))

#Subset ASVs in the network
node_names_in_network <- net[["labels"]][["labels1"]] # Extract node names
asv_to_phylum_filtered <- phylum_labels[node_names_in_network]   # Subset ASVs

#Get unique phyla in the filtered network
unique_phyla <- sort(unique(asv_to_phylum_filtered))

#Reorder the color vector to match phyla
n_phyla <- length(unique_phyla)
extended_paletten1 <- pal_adaptive("lancet", "lanonc")(17)
noCol <- c( "#0099B4FF", "#497BA9FF","#D5575DFF","#AD002AFF")
extended_paletten <- setdiff(extended_paletten1,noCol)
legend_colors <- setNames(extended_paletten,unique_phyla)
legend_colors_reordered <- legend_colors[unique_phyla]


png("AvB_ASV_Var_PhylCol25_circle.png", width = 183, height = 89,units = "mm",res = 300)
set.seed(11092003)
net <- plot(AvB_props_var, 
            mar = c(7,3,5,3),
            layout = 'layout_in_circle',
            sameLayout = TRUE, 
            layoutGroup = "union", 
            nodeFilter = "highestClose",
            nodeFilterPar = 25,
            nodeColor = "feature",
            featVecCol = asv_to_phylum_filtered,
            nodeShape = c("circle","triangle"),
            featVecShape = node_shape,
            colorVec =legend_colors_reordered,
            sameFeatCol =TRUE,
            nodeSize = "eigenvector",
            labels = list(new_labels1,new_labels2),
            labelScale = FALSE,
            rmSingles = "inboth",
            cexNodes = 0.9, 
            cexLabels = 0.5,
            nodeTransp = 0,
            posCol = "darkturquoise", 
            negCol = "orange",
            cexTitle = 1.5,
            showTitle = TRUE,
            groupNames = c("A", "B"),
            hubBorderCol  = "gray40")

legend(bty="n",
  "bottom", 
  legend = unique_phyla, 
  fill = legend_colors_reordered, 
  cex=0.75,
  text.width = 0.35, 
  xpd = TRUE, 
  ncol = 5,
  y.intersp = 0.5
  
)
mtext("Phylum", side = 1, line = 0, cex = 0.75)

dev.off()

#____reorder----
adj_mat <- AvB_props_var[["input"]][["adjaMat1"]]
node_names <- colnames(adj_mat)
taxonomy_info <- tot_tax[match(node_names, tot_tax$ASV), "Phylum"]
taxonomy_info <- pull(taxonomy_info, Phylum)
new_order <- node_names[order(taxonomy_info)]

#____reconstruct ordered network----
Root_Ects<- as.data.frame(read_tsv("Root_ASV_Ects.tsv",lazy = FALSE))%>%
  select(-any_of(remove_fung))%>%
  tibble::column_to_rownames("seq_name")
sample_order <- row.names(Root_Ects)
sample_types <- metadata$`A/B`[match(sample_order, metadata$seq_name)]
sample_types_f <- factor(sample_types)
levels(sample_types_f)

f_Root_Ects <- select(as.data.frame(Root_Ects), any_of(node_names))
of_Root_Ects <- f_Root_Ects[,new_order]
of_Root_Ects <-as.matrix(of_Root_Ects)

AvB_net_var_or <- netConstruct(of_Root_Ects,
                            group = sample_types,
                            dataType="counts",
                            filtTax = "none",
                            filtSamp = "none",
                            measure = "spieceasi",
                            measurePar = list(nlambda=20, 
                                              pulsar.params=list(rep.num=50)),
                            normMethod = "none",#"clr", #treated internally in SPIEC-EASI
                            zeroMethod = "none",#"multRepl", #treated internally in SPIEC-EASI
                            sparsMethod = "none", 
                            dissFunc = "signed",
                            verbose = 3,#notifying on each step
                            seed = 11092003)

AvB_props_var_or <- netAnalyze(AvB_net_var_or, 
                            centrLCC = FALSE,
                            clustMethod = "cluster_louvain",
                            hubPar = c("eigenvector"))
summary(AvB_props_var_or)

#____plot reordered color by taxa----
#change labels
current_labels1 <- AvB_props_var_or[["lccNames1"]]
current_labels2 <- AvB_props_var_or[["lccNames2"]]

taxonomy_df1 <- tot_tax %>%
  filter(ASV %in% current_labels1) %>%
  select(ASV,Kingdom,Phylum,Class, Order,Family)%>%  #%>% # Select only the necessary columns
  mutate(
    Class = ifelse(is.na(Class), paste0("Unknown_", Phylum), Class),
    Order = ifelse(is.na(Order), paste0("Unknown_", Class), Order),
    Family = ifelse(is.na(Family), paste0("Unknown_", Order), Family))%>%
  unique()

taxonomy_df2 <- tot_tax %>%
  filter(ASV %in% current_labels2) %>%
  select(ASV, Phylum)%>%
  unique()

# Create a named vector for easy replacement
label_mapping1 <- setNames(taxonomy_df1$Phylum, taxonomy_df1$ASV)
label_mapping2 <- setNames(taxonomy_df2$Phylum, taxonomy_df2$ASV)

# Replace the labels in the network
new_labels1 <- label_mapping1[current_labels1]
new_labels2 <- label_mapping2[current_labels2]

phylum_mapping <- setNames(taxonomy_df1$Phylum, taxonomy_df1$ASV)
phylum_labels <- phylum_mapping[current_labels1]

class_mapping <- setNames(taxonomy_df1$Class, taxonomy_df1$ASV)
class_labels <- class_mapping[current_labels1]

Order_mapping <- setNames(taxonomy_df1$Order, taxonomy_df1$ASV)
Order_labels <- Order_mapping[current_labels1]

shape <- select(taxonomy_df1, ASV, Kingdom)%>%
  mutate(shape=ifelse(Kingdom=="Bacteria", 16, 17))

node_shape <- factor(setNames(shape$shape,shape$ASV))

#Subset ASVs in the network
node_names_in_network <- net[["labels"]][["labels1"]] # Extract node names
asv_to_phylum_filtered <- phylum_labels[node_names_in_network]   # Subset ASVs

#Get unique phyla in the filtered network
unique_phyla <- sort(unique(asv_to_phylum_filtered))

#Reorder the color vector to match phyla
n_phyla <- length(unique_phyla)
extended_paletten1 <- pal_adaptive("lancet", "lanonc")(19)
noCol <- c( "#0099B4FF", "#497BA9FF","#D5575DFF","#AD002AFF")
extended_paletten <- setdiff(extended_paletten1,noCol)
legend_colors <- setNames(extended_paletten1,unique_phyla)
legend_colors_reordered <- legend_colors[unique_phyla]


png("AvB_ASV_Var_PhylCol25_circle_or.png", width = 183, height = 89,units = "mm",res = 300)
set.seed(11092003)
net <- plot(AvB_props_var_or, 
            mar = c(8,3,5,3),
            layout = 'layout_in_circle',
            sameLayout = TRUE, 
            layoutGroup = "union", 
            nodeFilter = "highestClose",
            nodeFilterPar = 25,
            nodeShape = c("circle","triangle"),
            featVecShape = node_shape,
            nodeSize = "eigenvector",
            labels = FALSE,
            labelScale = FALSE,
            rmSingles = "inboth",
            cexNodes = 0.9, 
            cexLabels = 0.5,
            nodeTransp = 0,
            #cexHubLabels = 1,
            posCol = "darkturquoise", 
            negCol = "orange",
            cexTitle = 0.75,
            showTitle = TRUE,
            #edge.curved = -10,
            groupNames = c("A", "B"),
            hubBorderCol  = "gray40")


# Define bacterial and fungal phyla
bacterial_phyla <- filter(taxonomy_df1, Kingdom == "Bacteria" & Phylum %in% unique_phyla)%>%
  pull(Phylum)%>%
  unique()
fungal_phyla <- filter(taxonomy_df1, Kingdom == "Fungi" & Phylum %in% unique_phyla)%>%
  pull(Phylum)%>%
  unique()

# Separate colors for bac and fung
bacteria_colors <- legend_colors_reordered[names(legend_colors_reordered) %in% bacterial_phyla]
fungi_colors <- legend_colors_reordered[names(legend_colors_reordered) %in% fungal_phyla]

legend(bty="n","bottomleft", 
       legend = names(bacteria_colors), 
       fill = bacteria_colors, 
       text.width = 0.35,
       #x=-0.5,y=-1.1,
       ncol = 4,y.intersp = 0.51,
       xpd=TRUE,
       title = "Bacterial Phyla", cex = 0.75)

# Add Fung Legend
legend("bottomright", 
       legend = names(fungi_colors), 
       fill = fungi_colors, 
       text.width = 0.35,
       #x=0.1,y=-1.1,
       y.intersp = 0.51,
       xpd=TRUE,
       title = "Fungal Phyla", cex = 0.75, bty = "n")
dev.off()

legend(bty="n",
              "bottom", # Position of the legend
               legend = unique_phyla, # Phyla names
               fill = legend_colors_reordered, # Corresponding colors
               cex=0.75,
               #  pt.cex = 0.01, 
               text.width = 0.35, 
               xpd = TRUE, 
               ncol = 5,
               y.intersp = 0.5
               
        )
