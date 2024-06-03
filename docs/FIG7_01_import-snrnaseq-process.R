# ## ######################################## ## #
#            IMPORT, PROCESS SNRNA-SEQ           #
# ## ######################################## ## #

# Date: Fri May 31 13:17:16 2024 ------------------
# Updated by: Daniela M. Roth

# Initial analysis utilized PopsicleR


# Load packages, functions and set directories ---------------------------------------
# Clean up code, add save scripts

source(here::here("docs/packages.R"))
source(here::here("docs/functions.R"))
source(here::here("docs/themes.R"))

seq_dir <- here::here("raw-data/snRNA-seq")
dir.create(here(fig7 <-
                  here("Figure-7")), recursive = TRUE) 


ins.data <- Read10X(data.dir = seq_dir)
ins <- CreateSeuratObject(counts = ins.data, project = "ins",
                          min.cells = 3, min.features = 200)
ins

# Preprocessing plots
prepro.plots(data = ins, output_dir = fig7)

# Based on violin plots, selected nFeature_RNA 200-4000, nCount_RNA > 200, percent.mito <0.15
ins <- add_percent_mito(ins)

filtered_ins <- subset(x = ins,
                       subset = (nFeature_RNA > 200) &
                         (nFeature_RNA < 4000) &
                         (nCount_RNA > 200) &
                         (percent.mito < 0.15))

filtered_ins <- genelvlfilt(filtered_ins)



# Pre-process object
filtered_ins <- NormalizeData(filtered_ins)
filtered_ins <- FindVariableFeatures(filtered_ins)
filtered_ins <- ScaleData(filtered_ins)
filtered_ins <- RunPCA(filtered_ins, nfeatures.print = 10)

# Find significant PCs
stdv <- filtered_ins[["pca"]]@stdev
sum.stdv <- sum(filtered_ins[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv)*100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) -1 ] -
                     percent.stdv[2:length(percent.stdv)])>0.1),
                  decreasing = T)[1] +1
min.pc <- min(co1, co2)
min.pc
# min pc calculated to be 10

# finish pre-processing
filtered_ins <- RunUMAP(filtered_ins, dims = 1:min.pc)
filtered_ins <- FindNeighbors(object = filtered_ins, dims = 1:min.pc)
filtered_ins <- FindClusters(object = filtered_ins, resolution = 0.1)

# pK identification (no ground-truth)
sweep.list <- paramSweep(filtered_ins, PCs = 1:min.pc)
sweep.stats <- summarizeSweep(sweep.list)
bcmvn <- find.pK(sweep.stats)

# Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

# Homotypic doublet proportion estimate
annotations <- filtered_ins@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp.poi <- round(optimal.pk * nrow(filtered_ins@meta.data))
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

# run DoubletFinder
doubletfinder_ins <- doubletFinder(seu = filtered_ins,
                              PCs = 1:min.pc,
                              pK = optimal.pk,
                              nExp = nExp.poi.adj)


colnames(doubletfinder_ins@meta.data)[8] <- "doublet_finder"

# subset and save
ins.singlets <- subset(doubletfinder_ins, doublet_finder == "Singlet")

# compare original and filtered
filtered_ins
ins.singlets

ins <- ins.singlets

# # Regress cell cycle
# s.genes = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
# g2m.genes = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
# 
# ins <- CellCycleScoring(ins, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# 
# ins <- RunPCA(ins, features = c(s.genes, g2m.genes))
# DimPlot(ins)
# 
# ins <- ScaleData(ins, vars.to.regress = c("S.Score","G2M.Score"),features = rownames(ins))
# 
# ins <- RunPCA(ins, features = c(s.genes, g2m.genes))
# DimPlot(ins)

ins <- FindNeighbors(ins, dims = 1:10)
ins <- FindClusters(ins, resolution = 0.8)
ins <- RunUMAP(ins, dims = 1:10)
DimPlot(ins, reduction = "umap",alpha = 0.5)

# Percent Difference in Expression
# Basic FindAllMarkers DE test
all_markers_pct <- FindAllMarkers(ins,verbose = T) %>% 
  Add_Pct_Diff()

all_markers_pct <- all_markers_pct %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  arrange(cluster)

write.csv(all_markers_pct, file = here(fig7, "all_markers_pct.csv"))

# Extract the top N marker genes per cluster for plotting
top_5 <- Extract_Top_Markers(marker_dataframe = all_markers_pct, num_genes = 5, rank_by = "avg_log2FC")


top50_markers_pct <- all_markers_pct %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>% 
  top_n(n=50, wt = avg_log2FC)%>%
  arrange(cluster)

write.csv(top50_markers_pct, file = here(fig7, "top50_markers_pct.csv"))


(plot <- DotPlot(
  object = ins,
  features =   top_5,
  scale.by = "radius",
  dot.scale = 8,
  split.by = NULL,
  cluster.idents = FALSE,
) + scale_colour_gradient2(low = "dodgerblue",
                           mid = "floralwhite",
                           high = "red2") +  custom_dotplot_theme() +RotatedAxis())


# List of genes for feature plots
progenitors <- c("Axin2","Gli1","Prrx1","Six2")
osteogenic <- c("Crabp1","Runx2","Sp7","Dmp1")
chondrogenic <- c("Col2a1","Acan","Mgp","Sox9")
osteoclasts <- c("Ctsk","Mmp9","Pheta1","Cd44")
vascular <- c("Mcam","Vwf","Pecam1","Pdgfrb")
myeloid_lymphocyte <- c("Pou2f2","Il1rl1","Gata2")
neurons_gli1 <- c("Neurod1","Cplx3","Otx2","Gfra3","Sox10","Foxd3")
erythrocytes <- c("Hba-a1","Hba-a2","Hbb-bs","Gypa","Gybp","Alas2","Klf1","Slc25a37","Slc2a1")
smoothmuscle <- c("Acta2","Tagln","Myh11","Des")


(plot <-
    FeaturePlot_scCustom(
      seurat_object = ins,
      reduction = "umap",
      na_cutoff = 0,
      features = progenitors,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = ins,
      reduction = "umap",
      na_cutoff = 0,
      features = osteogenic,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = ins,
      reduction = "umap",
      na_cutoff = 0,
      features = chondrogenic,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = ins,
      reduction = "umap",
      na_cutoff = 0,
      features = osteoclasts,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = ins,
      reduction = "umap",
      na_cutoff = 0,
      features = vascular,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = ins,
      reduction = "umap",
      na_cutoff = 0,
      features = myeloid_lymphocyte,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = ins,
      reduction = "umap",
      na_cutoff = 0,
      features = neurons_gli1,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = ins,
      reduction = "umap",
      na_cutoff = 0,
      features = erythrocytes,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = ins,
      reduction = "umap",
      na_cutoff = 0,
      features = smoothmuscle,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

dittoDimPlot(ins, "ident",opacity=1,do.ellipse=FALSE,do.label = TRUE,labels.size = 2.8,labels.repel = TRUE,labels.highlight = FALSE,do.raster=TRUE,raster.dpi = 600,size = 0.7,order = "randomize",legend.size =3,main=NULL,
             color.panel = c("darkolivegreen2", "orange", "purple", "lightcoral", "skyblue","maroon2","slateblue","gold","dodgerblue3","plum1","darkseagreen3","pink"))


# Create simple annotation files
# Create_Cluster_Annotation_File(file_path = fig7, file_name = "cluster_annotation")

annotation_info <- Pull_Cluster_Annotation(annotation = here(fig7,"cluster_annotation.csv"))

# Rename clusters
obj_annot <- Rename_Clusters(seurat_object = ins, new_idents = annotation_info$new_cluster_idents, meta_col_name = "mes_annotation")

Idents(obj_annot) <- factor(x = Idents(obj_annot), levels = sort(levels(obj_annot)))
obj_annot$mesenchyme <- Idents(obj_annot)


levels(x=obj_annot)=c("mes.prog.1","mes.prog.2","osteog.1","osteog.2","osteog.3","osteog.4","chondro","hemato","vasc","epith")

dittoDimPlot(obj_annot, "ident",opacity=1,do.ellipse=FALSE,do.label = TRUE,labels.size = 2.8,labels.repel = TRUE,labels.highlight = FALSE,do.raster=TRUE,raster.dpi = 600,size = 0.7,order = "randomize",legend.size =3,main=NULL,
             color.panel = c("darkolivegreen2", "orange", "lightcoral", "skyblue","maroon2","slateblue","dodgerblue3","plum1","darkseagreen3","pink"))

# The following are based on previous calculation of markers: repeat and replace with actual markers from this new analysis

ins <- obj_annot


# Basic FindAllMarkers DE test
all_markers_pct <- FindAllMarkers(ins,verbose = T) %>% 
  Add_Pct_Diff()

all_markers_pct <- all_markers_pct %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  arrange(cluster)

write.csv(all_markers_pct, file = here(fig7, "all_markers_pct_annot.csv"))

# Extract the top N marker genes per cluster for plotting
top_5 <- Extract_Top_Markers(marker_dataframe = all_markers_pct, num_genes = 5, rank_by = "avg_log2FC")


top50_markers_pct <- all_markers_pct %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>% 
  top_n(n=50, wt = avg_log2FC)%>%
  arrange(cluster)

write.csv(top50_markers_pct, file = here(fig7, "top50_markers_pct_annot.csv"))


(plot <- DotPlot(
  object = ins,
  features =   top_5,
  scale.by = "size",
  dot.scale = 5,
  cluster.idents = FALSE,
) + scale_colour_gradient2(low = "dodgerblue",
                           mid = "floralwhite",
                           high = "red2")+ theme(legend.position="bottom",
                                                legend.text = element_text(size = 6),
                                                legend.title = element_text(size=6),
                                                legend.justification = "left",
                                                axis.text = element_text(size = 8),
                                                axis.text.x = element_text(face = "italic",angle=45,hjust = 1),
                                                axis.title = element_text(size = 8),
                                                axis.title.x = element_text(hjust = 0.5)) )
  



# Chondrogenic markers

chondro.markers <- FindMarkers(ins,ident.1 = c("chondro"))
head(chondro.markers,n=5)

my_genes <- c("Col2a1","Col9a1","Meg3","Lama2")

plots <- FeaturePlot(ins,features=c(my_genes),pt.size = 0.5,ncol = 1,combine=FALSE,cols = c("azure","firebrick2"))

plots <- lapply(X = plots, FUN = function(x) x + theme(plot.title = element_text(size = 5,face = "italic"),text = element_text(size=5),axis.text = element_text(size=4)))

CombinePlots(plots = plots)

# Osteogenic markers

osteog.markers <- FindMarkers(ins,ident.1 = c("osteog.1","osteog.2","osteog.3","osteog.4"))
head(osteog.markers,n=5)

my_genes <- c("Runx2","Satb2","Cdh2","Col1a1")

plots <- FeaturePlot(ins,features=c(my_genes),pt.size = 0.5,ncol = 1,combine=FALSE,cols = c("azure","firebrick2"))

plots <- lapply(X = plots, FUN = function(x) x + theme(plot.title = element_text(size = 5,face = "italic"),text = element_text(size=5),axis.text = element_text(size=4)))

CombinePlots(plots = plots)

# Osteogenic markers first 2 clusters

markers <- FindMarkers(ins,ident.1 = c("osteog.1","osteog.2"))
head(markers,n=5)

my_genes <- c("Mmp13","Cp","Vcan","Tnc")

plots <- FeaturePlot(ins,features=c(my_genes),pt.size = 0.5,ncol = 1,combine=FALSE,cols = c("azure","firebrick2"))

plots <- lapply(X = plots, FUN = function(x) x + theme(plot.title = element_text(size = 5,face = "italic"),text = element_text(size=5),axis.text = element_text(size=4)))

CombinePlots(plots = plots)

# Last 2 osteogenic clusters

markers <- FindMarkers(ins,ident.1 = c("osteog.3","osteog.4"))
head(markers,n=5)


my_genes <- c("Col1a2","Col1a1","Col24a1","Smpd3")

plots <- FeaturePlot(ins,features=c(my_genes),pt.size = 0.5,ncol = 1,combine=FALSE,cols = c("azure","firebrick2"))

plots <- lapply(X = plots, FUN = function(x) x + theme(plot.title = element_text(size = 5,face = "italic"),text = element_text(size=5),axis.text = element_text(size=4)))

CombinePlots(plots = plots)

# Mesenchymal progenitor cluster markers

markers <- FindMarkers(ins,ident.1 = c("mes.prog.1","mes.prog.2"))
head(markers,n=5)


my_genes <- c("Ptn","Dlg2","Col3a1","Pdzrn4")

plots <- FeaturePlot(ins,features=c(my_genes),pt.size = 0.5,ncol = 1,combine=FALSE,cols = c("azure","firebrick2"))

plots <- lapply(X = plots, FUN = function(x) x + theme(plot.title = element_text(size = 5,face = "italic"),text = element_text(size=5),axis.text = element_text(size=4)))

CombinePlots(plots = plots)

# Mes prog 1 markers

markers <- FindMarkers(ins,ident.1 = c("mes.prog.1"))
head(markers,n=5)

my_genes <- c("Col14a1","Meg3","Frem1","Rarb")

plots <- FeaturePlot(ins,features=c(my_genes),pt.size = 0.5,ncol = 1,combine=FALSE,cols = c("azure","firebrick2"))

plots <- lapply(X = plots, FUN = function(x) x + theme(plot.title = element_text(size = 5,face = "italic"),text = element_text(size=5),axis.text = element_text(size=4)))

CombinePlots(plots = plots)

# Mes prog 2 markers

markers <- FindMarkers(ins,ident.1 = c("mes.prog.2"))
head(markers,n=5)

my_genes <- c("Ptn","Pdzrn4","Postn","Dlg2")

plots <- FeaturePlot(ins,features=c(my_genes),pt.size = 0.5,ncol = 1,combine=FALSE,cols = c("azure","firebrick2"))

plots <- lapply(X = plots, FUN = function(x) x + theme(plot.title = element_text(size = 5,face = "italic"),text = element_text(size=5),axis.text = element_text(size=4)))

CombinePlots(plots = plots)

# Hematopoietic markers

markers <- FindMarkers(ins,ident.1 = c("hemato"))
head(markers,n=5)

my_genes <- c("St18","Slc9b2","Csf1r","Atp6v0d2")

plots <- FeaturePlot(ins,features=c(my_genes),pt.size = 0.5,ncol = 1,combine=FALSE,cols = c("azure","firebrick2"))

plots <- lapply(X = plots, FUN = function(x) x + theme(plot.title = element_text(size = 5,face = "italic"),text = element_text(size=5),axis.text = element_text(size=4)))

CombinePlots(plots = plots)

# Vascular markers
markers <- FindMarkers(ins,ident.1 = c("vasc"))
head(markers,n=5)

my_genes <- c("Ptprb","Pecam1","Egfl7","Kdr")

plots <- FeaturePlot(ins,features=c(my_genes),pt.size = 0.5,ncol = 1,combine=FALSE,cols = c("azure","firebrick2"))

plots <- lapply(X = plots, FUN = function(x) x + theme(plot.title = element_text(size = 5,face = "italic"),text = element_text(size=5),axis.text = element_text(size=4)))

CombinePlots(plots = plots)

# Olfactory epithelium markers
markers <- FindMarkers(ins,ident.1 = c("epith"))
head(markers,n=5)

my_genes <- c("Dnah12","Dnah5","Dnah3","Ccdc180")

plots <- FeaturePlot(ins,features=c(my_genes),pt.size = 0.5,ncol = 1,combine=FALSE,cols = c("azure","firebrick2"))

plots <- lapply(X = plots, FUN = function(x) x + theme(plot.title = element_text(size = 5,face = "italic"),text = element_text(size=5),axis.text = element_text(size=4)))

CombinePlots(plots = plots)


# dot plot of osteogenic/meaningful markers:
markerlist <- c("Col14a1","Frem1","Axin2","Muc5b","Gli1","Dlg2","Ptn","Postn","Col12a1","Prrx1","Runx2","Alpl","Dlx5","Sparc","Sp7","Bglap","Dmp1","Phex","Sox9","Acan","Col2a1","Slc9b2","Acp5","Ocstamp","Kdr","Pecam1","Egfl7","Dnah3","Dnah12","Rgs22")

DotPlot(object = ins, features =   markerlist,scale.by = "size",dot.scale = 4.6
) +
  scale_colour_gradient2(low="dodgerblue",mid="white",high = "red2")+
  theme(legend.position="bottom",
        legend.text = element_text(size = 6),
        legend.title = element_text(size=6),
        legend.justification = "left",
        axis.text = element_text(size = 8),
        axis.text.x = element_text(face = "italic",angle=45,hjust = 1),
        axis.title = element_text(size = 8),
        axis.title.x = element_text(hjust = 0.5)) 


# Save annotated object
saveRDS(obj_annot, file = paste0(fig7,"/annot_ins.Rds"))

