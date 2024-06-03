# ## ######################################## ## #
#                 BMP7 SUBSET CELLS              #
# ## ######################################## ## #
# Date: Mon Jun 03 11:26:24 2024 ------------------
# Updated by: Daniela M. Roth

# Load annotated object
ins <- readRDS(here(fig7,"/annot_ins.Rds"))


SCpubr::do_NebulosaPlot(sample = ins,
                        features = "Bmp7",
                        pt.size =2,
                        viridis.direction = 1,
                        use_viridis = T,
                        viridis.palette = "plasma",
                        plot_cell_borders = FALSE,
                        font.size = 10,
                        plot.title.face = "italic")

# Determine top genes for Bmp7+ cells by subsetting cells

Bmp7_expression <- GetAssayData(object=ins,assay = 'RNA',layer="data")["Bmp7",]

pos_ids <- names(which(Bmp7_expression>0))
neg_ids <- names(which(Bmp7_expression==0))

pos_cells <- subset(ins,cells=pos_ids)
neg_cells <- subset(ins,cells = neg_ids)

p1 <- FeaturePlot(pos_cells,"Bmp7",pt.size = 0.5,cols = c("azure","firebrick2"))+ theme(plot.title = element_text(size = 5,face = "italic"),text = element_text(size=5),axis.text = element_text(size=4))
p2 <- FeaturePlot(neg_cells,"Bmp7",pt.size = 0.5,cols = c("azure","dodgerblue"))+ theme(plot.title = element_text(size = 5,face = "italic"),text = element_text(size=5),axis.text = element_text(size=4))

plot_grid(p1,p2,ncol=2)


markers <- FindAllMarkers(pos_cells,only.pos=FALSE,log.fc.threshold=0.25,test.use = "poisson")

top30 <- markers %>% 
  group_by(cluster) %>% 
  slice_max(n=30,order_by = avg_log2FC)

top10 <- markers %>% 
  group_by(cluster) %>% 
  slice_max(n=10,order_by = avg_log2FC)

write.csv(top30,file=here(fig7,"/top30_Bmp7cells.csv"))

genes_to_plot <- c("Tpm1","Parp8","Fbxl7","Bicc1","mt-Nd4","Rbms3","mt-Co2","mt-Co1","Zfp521","Svil")

DoHeatmap(pos_cells,slot = "data",angle=0,size=1.5,group.bar.height = 0,hjust=0.5,features=genes_to_plot)+scale_fill_gradientn(colors=c("dodgerblue","oldlace","mistyrose2","tomato2","red2"),na.value = "white")+theme(axis.text.y = element_text(face = "italic"),text = element_text(size=7),legend.position = "bottom")
          