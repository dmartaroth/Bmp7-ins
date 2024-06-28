# ## ######################################## ## #
#    IMPORT VISIUM DATA, PROCESS, SUBSET         #
# ## ######################################## ## #

# Date: Fri May 31 11:56:19 2024 ------------------
# Updated by: Daniela M. Roth

# Load packages and set directories
source(here::here("docs/packages.R"))

visium_dir <- here::here("raw-data/Visium")
dir.create(here(fig4 <-
                  here("Figure-4")), recursive = TRUE) 

# Assemble WT data --------------------------------------------------------

# Section 1: WT1
WT1_dir <- here::here(visium_dir,"WT_P3_INS/outs")
outs<- fig4
WT1 <- Load10X_Spatial(data.dir = WT1_dir, slice="WT1")

# Visualize slide
plot <- SpatialFeaturePlot(WT1, features = "nCount_Spatial")+theme(legend.position = "right")
ggsave(here("Figure-4/01a_visium_WT1_nCount_Spatial.pdf"), plot = plot,width = 5)

# Normalize
WT1<-NormalizeData(WT1, normalization.method = "RC", scale.factor = 10000)
WT1_coords<-GetTissueCoordinates(WT1)

# Get range of rows and columns
x_range <- range(WT1_coords$imagerow)
y_range <- range(WT1_coords$imagecol)
cat("row range: ",x_range[1], "to", x_range[2], "\n")
cat("col range: ",y_range[1], "to", y_range[2], "\n")


# Section 2: WT2
WT2_dir <- here::here(visium_dir,"Vx1_P3_WT/outs")
outs<- fig4
WT2 <- Load10X_Spatial(data.dir = WT2_dir, slice="WT2")

# Visualize slide
plot <- SpatialFeaturePlot(WT2, features = "nCount_Spatial")+theme(legend.position = "right")
ggsave(here("Figure-4/01b_visium_WT2_nCount_Spatial.pdf"), plot = plot,width = 5)

# Normalize
WT2<-NormalizeData(WT2, normalization.method = "RC", scale.factor = 10000)
WT2_coords<-GetTissueCoordinates(WT2)

# Get range of rows and columns
x_range <- range(WT2_coords$imagerow)
y_range <- range(WT2_coords$imagecol)
cat("row range: ",x_range[1], "to", x_range[2], "\n")
cat("col range: ",y_range[1], "to", y_range[2], "\n")


# Subset endocranial regions ----------------------------------------------

# 3 endocranial regions for slide WT1, 1 for WT2
# Endocranial part of WT2 region 2 is wiped off, and region 3 ins is off the barcodable region
# Parameters for imagerow/col and gene thresholds are individually set to 
## manually select gems interactively


## WT1 endocranial region 1 ------------------------------------------------

WT1endoR1Par<-subset(WT1, subset = wt1_imagerow>165 & wt1_imagerow<180&
                       wt1_imagecol>200&wt1_imagecol<250&
                       (Col1a1>980&Sp7>1&Col1a1 <1250&Col2a1<70))

# Visualize subset
plot <-
  SpatialFeaturePlot(
    WT1endoR1Par,
    features = "Col1a1",
    crop = FALSE,
    interactive = FALSE
  )+ theme(legend.position = "right")

ggsave(here("Figure-4/02a_visium_WT1_endo1.png"), plot = plot,width = 5)


## WT1 endocranial region 2 ------------------------------------------------

WT1endoR2Par<-subset(WT1, subset = wt1_imagerow>150 & wt1_imagerow<165&
                       wt1_imagecol>480&wt1_imagecol<520&
                       (Col2a1<90&Col1a1>750))

# Visualize subset
plot <-
  SpatialFeaturePlot(
    WT1endoR2Par,
    features = "Col1a1",
    crop = FALSE,
    interactive = FALSE
  )+ theme(legend.position = "right")

ggsave(here("Figure-4/02b_visium_WT1_endo2.png"), plot = plot,width = 5)


## WT1 endocranial region 3 ------------------------------------------------

WT1endoR3Par<-subset(WT1, subset = wt1_imagerow>390 & wt1_imagerow<400&
                       wt1_imagecol>475&wt1_imagecol<485)

# Visualize subset
plot <-
  SpatialFeaturePlot(
    WT1endoR3Par,
    features = "Col1a1",
    crop = FALSE,
    interactive = FALSE
  )+ theme(legend.position = "right")

ggsave(here("Figure-4/02c_visium_WT1_endo3.png"), plot = plot,width = 5)


## WT2 endocranial region 1 ------------------------------------------------

WT2endoR1Par<-subset(WT2, subset = wt2_imagerow>250 & wt2_imagerow<260&
                       wt2_imagecol>500&wt2_imagecol<900&
                       (Bglap<2))

# Visualize subset
plot <-
  SpatialFeaturePlot(
    WT2endoR1Par,
    features = "Col1a1",
    crop = FALSE,
    interactive = FALSE
  )+ theme(legend.position = "right")

ggsave(here("Figure-4/02d_visium_WT2_endo1.png"), plot = plot,width = 5)



## Assign identities to reach region ---------------------------------------
Idents(object = WT1endoR1Par) <-"WT1endoR1"
Idents(object = WT1endoR2Par) <-"WT1endoR2"
Idents(object = WT1endoR3Par) <-"WT1endoR3"
Idents(object = WT2endoR1Par) <-"WT2endoR1"


# Merge endocranial regions for WT1 and WT2 -------------------------------

# Merge WT1 endocranial regions together
WT1endoPar<-merge(WT1endoR1Par, y=c(WT1endoR2Par,WT1endoR3Par))

# Visualize endocranial WT1 merge
plot <- SpatialFeaturePlot(WT1endoPar, features = "Col1a1", crop = FALSE) 

ggsave(here("Figure-4/03a_visium_merged_WT1_endo.pdf"), plot = plot,width = 8, height = 4)

plot <- VlnPlot(WT1endoPar, features = c("Col1a1","Col2a1","Ifitm5", "Pecam1"), pt.size = 0.3, ncol=4,  log = FALSE) + NoLegend()
ggsave(here("Figure-4/03b_vln_visium_merged_WT1_endo.pdf"), plot = plot,width = 8, height = 4)

# Merge WT1 merge with WT2 endocranial region
WTendoPar <- merge(WT1endoPar, y=c(WT2endoR1Par))

# Visualize endocranial WT merge
plot <- SpatialFeaturePlot(WTendoPar, features = "Col1a1", crop = FALSE) 
ggsave(here("Figure-4/03c_visium_merged_WT_endo.pdf"), plot = plot,width = 9, height = 5)

plot <- VlnPlot(WTendoPar, features = c("Col1a1","Col2a1","Ifitm5", "Pecam1"), pt.size = 0.3, ncol=4,  log = FALSE) + NoLegend()
ggsave(here("Figure-4/03d_vln_visium_merged_WT_endo.pdf"), plot = plot,width = 8, height = 4)


# Subset ectocranial regions ----------------------------------------------

# Ectocranial region of WT2 R1 is cut off of barcodable region
# WT2 ectocranial region 3 off region of sequencing

## WT1 ectocranial region 1 ------------------------------------------------


WT1ectoR1Par<-subset(WT1, subset = wt1_imagerow>180 & wt1_imagerow<210&
                       wt1_imagecol>190&wt1_imagecol<260&
                       (Col1a1>980&Sp7>1&Col1a1 <1250&Col2a1<70))

# Visualize subset
plot <-
  SpatialFeaturePlot(
    WT1ectoR1Par,
    features = "Col1a1",
    crop = FALSE,
    interactive = FALSE
  )+ theme(legend.position = "right")

ggsave(here("Figure-4/04a_visium_WT1_ecto1.png"), plot = plot,width = 5)


## WT1 ectocranial region 2 ------------------------------------------------
WT1ectoR2Par<-subset(WT1, subset = wt1_imagerow>162 & wt1_imagerow<190&
                       wt1_imagecol>430&wt1_imagecol<485&
                       (Col2a1<90&Col1a1>700))
# ISSUE MUST BE SOLVED BEFORE CODE CAN BE RUN
plot <-
  SpatialFeaturePlot(
    WT1ectoR2Par,
    features = "Col1a1",
    crop = FALSE,
    interactive = FALSE
  )+ theme(legend.position = "right")

ggsave(here("Figure-4/04b_visium_WT1_ecto2.png"), plot = plot,width = 5)
# plotting not working, something to do with incorrect dimensions


## WT1 ectocranial region 3 ------------------------------------------------
WT1ectoR3Par<-subset(WT1, subset = wt1_imagerow>375 & wt1_imagerow<390&
                       wt1_imagecol>475&wt1_imagecol<485)

# Visualize subset
plot <-
  SpatialFeaturePlot(
    WT1ectoR3Par,
    features = "Col1a1",
    crop = FALSE,
    interactive = FALSE
  )+ theme(legend.position = "right")

ggsave(here("Figure-4/04c_visium_WT1_ecto3.png"), plot = plot,width = 5)


## WT1 ectocranial region 4 ------------------------------------------------
WT1ectoR4Par<-subset(WT1, subset = wt1_imagerow>55&wt1_imagerow<68&
                       wt1_imagecol>330&wt1_imagecol<360&
                       (Ctsk<37&Pecam1>0))

# Visualize subset
plot <-
  SpatialFeaturePlot(
    WT1ectoR4Par,
    features = "Col1a1",
    crop = FALSE,
    interactive = FALSE
  )+ theme(legend.position = "right")

ggsave(here("Figure-4/04d_visium_WT1_ecto4.png"), plot = plot,width = 5)


## WT2 ectocranial region 2 ------------------------------------------------
WT2ectoR2Par<-subset(WT2, subset = wt2_imagerow>500 & wt2_imagerow<530&
                       wt2_imagecol>490&wt2_imagecol<500&
                       (Bglap<2))

# Visualize subset
plot <-
  SpatialFeaturePlot(
    WT2ectoR2Par,
    features = "Col1a1",
    crop = FALSE,
    interactive = FALSE
  )+ theme(legend.position = "right")

ggsave(here("Figure-4/04e_visium_WT2_ecto2.png"), plot = plot,width = 5)


## Assign identities to each region ----------------------------------------
Idents(object = WT1ectoR1Par) <-"WT1ectoR1"
Idents(object = WT1ectoR2Par) <-"WT1ectoR2"
Idents(object = WT1ectoR3Par) <-"WT1ectoR3"
Idents(object = WT1ectoR4Par) <-"WT1ectoR4"
Idents(object = WT2ectoR2Par) <-"WT2ectoR2"


# Merge ectocranial regions -----------------------------------------------

# Merge ectocranial regions for WT1
WT1ectoPar<-merge(WT1ectoR1Par, y=c(WT1ectoR2Par,WT1ectoR3Par,WT1ectoR4Par))

# Visualize ectocranial WT1 merge
plot <- SpatialFeaturePlot(WT1ectoPar, features = "Col1a1", crop = FALSE) 

ggsave(here("Figure-4/05a_visium_merged_WT1_ecto.pdf"), plot = plot,width = 8, height = 4)

plot <- VlnPlot(WT1ectoPar, features = c("Col1a1","Col2a1","Ifitm5", "Pecam1"), pt.size = 0.3, ncol=4,  log = FALSE) + NoLegend()
ggsave(here("Figure-4/05b_vln_visium_merged_WT1_ecto.pdf"), plot = plot,width = 8, height = 4)


# Merge ectocranial WT1 and WT2 regions -------------------------------------------

WTectoPar <- merge(WT1ectoPar, y=c(WT2ectoR2Par))

# Visualize endocranial WT merge
plot <- SpatialFeaturePlot(WTectoPar, features = "Col1a1", crop = FALSE) 
ggsave(here("Figure-4/05c_visium_merged_WT_ecto.pdf"), plot = plot,width = 9, height = 5)

plot <- VlnPlot(WTendoPar, features = c("Col1a1","Col2a1","Ifitm5", "Pecam1"), pt.size = 0.3, ncol=4,  log = FALSE) + NoLegend()
ggsave(here("Figure-4/05d_vln_visium_merged_WT_ecto.pdf"), plot = plot,width = 8, height = 4)


# Merge endo- and ectocranial gems ----------------------------------------

All<-merge(WTendoPar, y=WTectoPar)

plot <- VlnPlot(All, features = c("Col1a1","Col2a1","Ifitm5", "Pecam1"), pt.size = 0.3, ncol=4,  log = FALSE) + NoLegend()
ggsave(here("Figure-4/06_vln_visium_merged_WT_all.pdf"), plot = plot,width = 8, height = 4)


## Group samples into endo or ecto -----------------------------------------

levels(All)
region <- c("ecto","ecto","ecto","ecto","endo","endo","endo","ecto","endo")
All$samples <- All@active.ident
names(region)<-levels(All)
All$region <- Idents(All)
All<-RenameIdents(All,region)
All@meta.data$region <-  All@meta.data$region
All$region <- Idents(All)
All@meta.data

# Reorder levels
levels(x = All) <- c("endo","ecto")
levels(All)


# Save merged object for further analysis ---------------------------------

saveRDS(All,file=here::here("data-output/Visium/P3_All_endovsecto.Rds"))

