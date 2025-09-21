## Notebook env: farnaz_spatial (R kernel)
## this notebook deconvolutes spots in visium data using the mBAT reference scRNAseq data
## to get proportion of cell types present in each spot in visium slide

library(BayesPrism)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

save_name = 'mBAT_slice_C1.rds'

mBAT <- readRDS("/home/ssobti/projects/farnaz_spatial/data/seurat_mBAT.rds")

mBAT@assays$RNA@counts -> sc_ref_mtx
t(sc_ref_mtx) -> sc_ref_mtx

cell.type.labels <- as.character(mBAT@meta.data$cluster_id)
cell.state.labels <- as.character(mBAT@meta.data$cluster_id)

spatial_seurat_obj <- Load10X_Spatial(
  "/home/ssobti/projects/farnaz_spatial/data/count-C1/outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "C1_slice",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL,
)

spatial_seurat_obj@assays$Spatial@counts -> spatial_count_mtx
t(spatial_count_mtx) -> spatial_count_mtx

as.matrix(sc_ref_mtx) -> sc_ref_mtx

as.matrix(spatial_count_mtx) -> spatial_count_mtx

plot.cor.phi (input=sc_ref_mtx,
                         input.labels=cell.state.labels,
                         title="cell state correlation",
                         #specify pdf.prefix if need to output to pdf
                         #pdf.prefix="gbm.cor.cs", 
                         cexRow=0.2, cexCol=0.2,
                         margins=c(2,2))

sc.stat <- plot.scRNA.outlier(
  input=sc_ref_mtx, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="mm", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE #return the data used for plotting. 
  #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)

spatial.stat <- plot.bulk.outlier(
  bulk.input=spatial_count_mtx,#make sure the colnames are gene symbol or ENSMEBL ID 
    sc.input=sc_ref_mtx, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="mm", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE
  #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
)

gene.list <- read.table(system.file("extdata", "genelist.mm.new.txt", package="BayesPrism"),sep="\t",header=F,stringsAsFactors=F)

genes_to_throw = gene.list$V3[gene.list$V1 %in% c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY")]

## filter out genes
sc.ref.filtered = sc_ref_mtx[, !(colnames(sc_ref_mtx) %in% genes_to_throw)]
cell_counts_above_zero <- apply(X = sc.ref.filtered, MARGIN = 2, FUN = function(x){length(which(x>0))})
sc.ref.filtered = sc.ref.filtered[, cell_counts_above_zero > 4]

library(org.Mm.eg.db)
library(AnnotationDbi)

gene_type_annotations <- mapIds(org.Mm.eg.db, keys = colnames(sc.ref.filtered),
       column = "GENETYPE", keytype = "SYMBOL")



unique(gene_type_annotations)

gene_type_df = data.frame(genes = names(gene_type_annotations), gene_type = as.character(gene_type_annotations))

gene_type_df <- gene_type_df[gene_type_df$gene_type %in% c('protein-coding', 'ncRNA', 'pseudo'),]

head(gene_type_df)

protein_coding_genes = (gene_type_df %>% filter(gene_type == 'protein-coding'))$genes
ncRNA_genes = (gene_type_df %>% filter(gene_type == 'ncRNA'))$genes
pseudo_genes = (gene_type_df %>% filter(gene_type == 'pseudo'))$genes

sc_total_exp = colSums(sc.ref.filtered)
spatial_total_exp = colSums(spatial_count_mtx)

gene_sets = list(protein_coding_genes, ncRNA_genes, pseudo_genes)
names(gene_sets) <- c('protein_coding_genes', 'ncRNA_genes', 'pseudo_genes')

names(gene_sets)

grphs = list()

for (i in 1:3){
    sc_ref <- sc_total_exp[gene_sets[[i]]]
    spatial_bulk <- spatial_total_exp[gene_sets[[i]]]
    sc_ref <- log2(sc_ref)
    spatial_bulk <- log2(spatial_bulk)
    graphing_df = data.frame(log2sc_ref = sc_ref, log2spatial_bulk = spatial_bulk)
    grphs[[i]] = ggplot(graphing_df, aes(log2sc_ref, log2spatial_bulk)) + geom_point() + geom_abline(slope = 1, color = 'red', linewidth = 1.5, linetype="dashed") + theme_classic() + ggtitle(names(gene_sets)[i])
}


options(repr.plot.width=20)
ggarrange(plotlist = grphs, nrow = 1)

sc.ref.filtered.pc <-  sc.ref.filtered[,protein_coding_genes]

dim(sc.ref.filtered)

dim(sc.ref.filtered.pc)

diff.exp.stat <- get.exp.stat(sc.dat=sc_ref_mtx[,colSums(sc_ref_mtx>0)>3],# filter genes to reduce memory use
                                          cell.type.labels=cell.type.labels,
                                          cell.state.labels=cell.state.labels,
                                          psuedo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                                          cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                                          n.cores=1 #number of threads
                                          )

sc.ref.filtered.pc.sig <- select.marker (sc.dat=sc.ref.filtered.pc,
                                                  stat=diff.exp.stat,
                                                  pval.max=0.01,
                                                  lfc.min=0.1)
#> number of markers selected for each cell type: 
#> tumor :  8686 
#> myeloid :  575 
#> pericyte :  114 
#> endothelial :  244 
#> tcell :  123 
#> oligo :  86

dim(sc.ref.filtered.pc.sig)
#> [1] 23793  7874

myPrism <- new.prism(
  reference=sc.ref.filtered.pc, 
  mixture=spatial_count_mtx,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key=NULL,
  outlier.cut=0.01,
    outlier.fraction=0.1,
)
#> number of cells in each cell state 
#> cell.state.labels
#> PJ017-tumor-6 PJ032-tumor-5     myeloid_8 PJ032-tumor-4 PJ032-tumor-3 PJ017-tumor-5         tcell PJ032-tumor-2 PJ030-tumor-5     myeloid_7 PJ035-tumor-8 PJ017-tumor-4 PJ017-tumor-3 PJ017-tumor-2 PJ017-tumor-0 PJ017-tumor-1 PJ018-tumor-5     myeloid_6     myeloid_5 PJ035-tumor-7         oligo PJ048-tumor-8 PJ032-tumor-1 PJ032-tumor-0 PJ048-tumor-7 PJ025-tumor-9 PJ035-tumor-6 PJ048-tumor-6 PJ016-tumor-6 PJ018-tumor-4     myeloid_4 PJ048-tumor-5 PJ048-tumor-4 PJ016-tumor-5 PJ025-tumor-8 PJ048-tumor-3 PJ035-tumor-5 PJ030-tumor-4 PJ018-tumor-3 PJ016-tumor-4 PJ025-tumor-7     myeloid_3 PJ035-tumor-4     myeloid_2 PJ025-tumor-6 PJ018-tumor-2 PJ030-tumor-3 PJ016-tumor-3 PJ025-tumor-5 PJ030-tumor-2 PJ018-tumor-1 PJ035-tumor-3 PJ048-tumor-2 PJ018-tumor-0 PJ048-tumor-1 PJ035-tumor-2 PJ035-tumor-1 PJ016-tumor-2 PJ030-tumor-1      pericyte   endothelial PJ035-tumor-0 PJ025-tumor-4     myeloid_1 PJ048-tumor-0     myeloid_0 PJ030-tumor-0 PJ025-tumor-3 PJ016-tumor-1 PJ016-tumor-0 PJ025-tumor-2 PJ025-tumor-1 PJ025-tumor-0 
#>            22            41            49            57            62            64            67            72            73            75            81            83            89           101           107           107           113           130           141           150           160           169           171           195           228           236           241           244           261           262           266           277           303           308           319           333           334           348           361           375           381           382           385           386           397           403           419           420           421           425           429           435           437           444           463           471           474           481           482           489           492           512           523           526           545           550           563           601           619           621           630           941           971 
#> Number of outlier genes filtered from mixture = 6 
#> Aligning reference and mixture... 
#> Nornalizing reference...

#Note that outlier.cut and outlier.fraction=0.1 filter genes in X whose expression fraction is greater than outlier.cut (Default=0.01) in more than outlier.fraction (Default=0.1) of bulk data. 
#Typically for dataset with reasonable quality control, very few genes will be filtered. 
#Removal of outlier genes will ensure that the inference will not be dominated by outliers, which sometimes may be resulted from poor QC in mapping.

bp.res <- run.prism(prism = myPrism, n.cores=5)

path = paste0('/home/ssobti/projects/farnaz_spatial/output_data/bayesprism/', save_name)
saveRDS(object = bp.res, file = path)

# extract posterior mean of cell type fraction theta
theta <- get.fraction (bp=bp.res,
            which.theta="final",
            state.or.type="type")

head(theta)

theta.cv <- bp.res@posterior.theta_f@theta.cv

head(theta.cv)




