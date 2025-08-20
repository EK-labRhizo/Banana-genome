#Root vs soil bacterial microbiome
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

Ext.ASVs4 <- read_tsv('Root_Ext.Date_16-03_Cult/significant_results.tsv')%>%
  filter(qval<0.2)%>%
  pull(feature)%>%
  unique()

Ext.ASVs5 <- read_tsv('Root_Ext.Date_16-3_Cult/significant_results.tsv')%>%
  filter(qval<0.2)%>%
  pull(feature)%>%
  unique()

Ext.ASVs6 <- read_tsv('Ext.Date_28-3_Soil/significant_results.tsv')%>%
  filter(qval<0.2)%>%
  pull(feature)%>%
  unique()

Ext.ASVs <- c(Ext.ASVs1, Ext.ASVs2,Ext.ASVs3,Ext.ASVs4,Ext.ASVs5,Ext.ASVs6)%>%
  unique() 

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

#find DAs----
#_Upload clr tables----
library(Maaslin2)
Root_clr <- read_tsv("rare_bacRoot_clr.tsv", lazy=FALSE)
Soil_clr <- read_tsv("rare_bacSoil_clr.tsv", lazy=FALSE)

full_clr <- as.data.frame(full_join(Root_clr,Soil_clr))
rownames(full_clr) <- full_clr$seq_name

input_metadata <- filter(metadata, seq_name %in% (full_clr$seq_name))%>%
  mutate(Sample = gsub("^[A-Za-z]", "", Sample))%>%
  select(seq_name, Type, Cultivar, Sample)
rownames(input_metadata) <- input_metadata$seq_name  

#_per sampled tree----

fit_data = Maaslin2(input_data     = full_clr_no.na, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    standardize = FALSE,
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "soilVroot_DA", 
                    fixed_effects  = c("Type"),
                    random_effects = c("Sample"))

DA_rootVsoil <- read_tsv("soilVroot_DA/significant_results.tsv", lazy=FALSE)%>%
  filter(qval < 0.05)%>%
  inner_join(tot_tax, by=c("feature"="ASV"))%>%
  write_tsv("DA_rootVsoil.tsv")

#Plot RvS DAs----
#_rootVsoil PCA loadings----
Root_clr <- read_tsv("rare_bacRoot_clr.tsv", lazy=FALSE)%>%
  tibble::column_to_rownames("seq_name")%>%
  t(.)%>%
  as.data.frame()%>%
  tibble::rownames_to_column("ASV")
Soil_clr <- read_tsv("rare_bacSoil_clr.tsv", lazy=FALSE)%>%
  tibble::column_to_rownames("seq_name")%>%
  t(.)%>%
  as.data.frame()%>%
  tibble::rownames_to_column("ASV")

full_clr_no.na <- inner_join(Root_clr,Soil_clr)%>%
  tibble::column_to_rownames("ASV")%>%
  t(.)%>%
  as.data.frame()

#_PERMANOVA----
sel_metadata <- filter(metadata, seq_name %in% unlist(rownames(full_clr_no.na)))
sel_metadata <- sel_metadata[match(rownames(full_clr_no.na), sel_metadata$seq_name), ]

adonis2(formula = full_clr_no.na ~ Type, data = sel_metadata, method = "euclidean", strata = sel_metadata$Cultivar)

#chck plot position
sel_pos <- filter(pos, seq_name %in% unlist(rownames(full_clr_no.na)))
sel_pos <- sel_pos[match(rownames(full_clr_no.na), sel_pos$seq_name), ]
adonis2(formula = full_clr_no.na ~ plot_Column, data = sel_pos, method = "euclidean")
adonis2(formula = full_clr_no.na ~ plot_Row, data = sel_pos, method = "euclidean")


pca.clr <- prcomp(full_clr_no.na, scale. = FALSE)
summary(pca.clr)
loadings <- as.data.frame(pca.clr$rotation)

PCA <- as.data.frame(pca.clr[["x"]])
PCA$seq_name <- unlist(row.names(PCA))
PCA <- inner_join(PCA, metadata)

PCA$Ploidity <- as.factor(PCA$Ploidity)
PCA$Column <- as.factor(PCA$Column)
PCA$Plate <- as.factor(PCA$Plate)
write_tsv(PCA, "bac_RvS_PCA.tsv")
root_DA <- read_tsv('DA_rootVsoil.tsv', lazy=FALSE)%>%
  filter(abs(coef)>1)%>%
  pull(feature)%>%
  unique()

RvS_loading <- tibble::rownames_to_column(loadings, "ASV")%>%
  filter(ASV %in%root_DA)%>%
  select(ASV:PC3)%>%
  left_join(tot_tax)%>%
  write_tsv('bac_RvS_PCA_loading.tsv')

RvS_loading <- read_tsv('bac_RvS_PCA_loading.tsv')

ggplot(data=PCA,aes(x=PC1, y=PC2)) +
  geom_point(aes(shape = Type), size = 4, na.rm = TRUE) +
  geom_point(data = filter(RvS_loading, abs(PC1)>0.03|abs(PC2)>0.03), aes(x = 1000*PC1, y = 1000*PC2,fill=Class),
             shape = 22, size=3, color="white") +
  scale_shape_discrete(name='Sample Type',labels = c("Rhizosphere", "Soil"))+
  scale_fill_adaptive(name = "lancet", palette = "lanonc")+
  #scale_color_jco()+
  #scale_fill_brewer(palette = "BrBG")+
  #scale_fill_manual(values = getPalette(11))+
  coord_fixed()+#(xlim = c(-50,50),ylim = c(-50,50)) +
  #facet_grid(time ~ soil) +
  #facet_grid(time ~ .) +
  labs(title=NULL,
       x="PC 1 (15% variance)",
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
ggsave("RvS_PCA12_load03.png", height = 180, width = 360, units = "mm", dpi=300)


#A/B only root V soil----

A_list <- filter(metadata, `A/B`=="A")%>%
  pull(seq_name)

B_list <- filter(metadata, `A/B`=="B")%>%
  pull(seq_name)

A_clr <- as.data.frame(filter(full_clr, seq_name%in%A_list))
B_clr <- as.data.frame(filter(full_clr, seq_name%in%B_list))

#_A list----
library(Maaslin2)

row.names(A_clr) <- A_clr$seq_name
names <- pull(A_clr, seq_name)
input_metadata <- as.data.frame(filter(metadata, seq_name%in%names))
row.names(input_metadata) <- input_metadata$seq_name

fit_data = Maaslin2(input_data     = A_clr, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    standardize = FALSE,
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "RvS_A_batch", 
                    fixed_effects  = c("Type"),
                    random_effects = c("Ext.Date"))

#_B list----
library(Maaslin2)

row.names(B_clr) <- B_clr$seq_name
names <- pull(B_clr, seq_name)
input_metadata <- as.data.frame(filter(metadata, seq_name%in%names))
row.names(input_metadata) <- input_metadata$seq_name

fit_data = Maaslin2(input_data     = B_clr, 
                    input_metadata = input_metadata, 
                    min_prevalence = -Inf,
                    min_abundance = -Inf,
                    normalization  = "NONE",
                    transform = "NONE",
                    standardize = FALSE,
                    plot_heatmap = FALSE,
                    plot_scatter = FALSE,
                    output         = "RvS_B_batch", 
                    fixed_effects  = c("Type"),
                    random_effects = c("Ext.Date"))

#_compare A to B----
A_sig_root <- read_tsv("RvS_A_batch/significant_results.tsv")%>%
  filter(qval < 0.05)%>%
  select(feature, coef)%>%
  left_join(tot_tax, by=c("feature"="ASV"))%>%
  mutate(`A/B`="A")%>%
  write_tsv("A_sig_root_batch.tsv")

B_sig_root <- read_tsv("RvS_B_batch/significant_results.tsv")%>%
  filter(qval < 0.05)%>%
  select(feature, coef)%>%
  left_join(tot_tax, by=c("feature"="ASV"))%>%
  mutate(`A/B`="B")%>%
  write_tsv("B_sig_root_batch.tsv")

#host selection strength----
library(vegan)
A_list <- filter(metadata, `A/B`=="A")%>%
  pull(seq_name)
B_list <- filter(metadata, `A/B`=="B")%>%
  pull(seq_name)
BaGen_cts <- read_tsv("rare_BaGen_cts.tsv", lazy=FALSE)
A_cts <- filter(BaGen_cts, seq_name%in%A_list)%>%
  tibble::column_to_rownames( "seq_name")
B_cts <- filter(BaGen_cts, seq_name%in%B_list)%>%
  tibble::column_to_rownames( "seq_name")

A_changed_list <- read_tsv("A_sig_root_batch.tsv", lazy = FALSE)%>%
  pull(feature)
B_changed_list <- read_tsv("B_sig_root_batch.tsv", lazy = FALSE)%>%
  pull(feature)

#_A host selection----
changed_ASVs <- select(A_cts, any_of(A_changed_list))

D_changed <-  diversity(changed_ASVs)
D_changed <-as.data.frame(D_changed)%>%
  tibble::rownames_to_column("seq_name")
A_R_changed <- as.data.frame(A_cts/rowSums(A_cts))%>%
  select(any_of(A_changed_list))%>%
  mutate(tot_relabund = rowSums(.))%>%
  select(tot_relabund)

A_R_cons <- as.data.frame(A_cts/rowSums(A_cts))%>%
  select(!any_of(A_changed_list))%>%
  mutate(tot_relabund = rowSums(.))%>%
  select(tot_relabund)
A_R_cons <-tibble::rownames_to_column(A_R_cons,"seq_name")
colnames(A_R_cons)[2] <- "tot_relbund_cons"

host_sel_A <- tibble::rownames_to_column(A_R_changed,"seq_name")%>%
  left_join(D_changed)%>%
  left_join(A_R_cons)%>%
  mutate(HSS = length(A_changed_list)*tot_relabund*D_changed/((ncol(A_cts)-ncol(changed_ASVs))*tot_relbund_cons))

host_sel_A <- 
  left_join(host_sel_A,select(metadata, seq_name,Type,`A/B`))

#_B host selection----
changed_ASVs <- select(B_cts, any_of(B_changed_list))

D_changed <-  diversity(changed_ASVs)
D_changed <-as.data.frame(D_changed)%>%
  tibble::rownames_to_column("seq_name")

B_R_changed <- as.data.frame(B_cts/rowSums(B_cts))%>%
  select(any_of(B_changed_list))%>%
  mutate(tot_relabund = rowSums(.))%>%
  select(tot_relabund)

B_R_cons <- as.data.frame(B_cts/rowSums(B_cts))%>%
  select(!any_of(B_changed_list))%>%
  mutate(tot_relabund = rowSums(.))%>%
  select(tot_relabund)

B_R_cons <-tibble::rownames_to_column(B_R_cons,"seq_name")
colnames(B_R_cons)[2] <- "tot_relbund_cons"

host_sel_B <- tibble::rownames_to_column(B_R_changed,"seq_name")%>%
  left_join(D_changed)%>%
  left_join(B_R_cons)%>%
  mutate(HSS = length(B_changed_list)*tot_relabund*D_changed/((ncol(B_cts)-ncol(changed_ASVs))*tot_relbund_cons))

host_sel_B <- 
  left_join(host_sel_B,select(metadata, seq_name,Type,`A/B`))
host_sel <- rbind(host_sel_A,host_sel_B)#%>%
#filter(Type == "Root")

#_statistical analysis----
kruskal.test(HSS~`A/B`, data=host_sel)

library(ggpubr)
host_sel%>%
  filter(Type=="Root")%>%
  ggplot(aes(x=`A/B`, y=HSS)) +
  geom_boxplot() +
  geom_point()+
  labs(title=NULL,
       x=NULL,
       y="Host Selection Strength") +
  theme_dbox()+
  stat_compare_means()+
  theme(axis.text.x = element_text(size = 12))
ggsave("hostSelection_AvB.png", height = 180, width = 180, units = "mm", dpi=300 )

host_sel%>%
  write_tsv("host_sel_stre.tsv")

