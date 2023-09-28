# Loading packages.
library(spatialHeatmap); library(SummarizedExperiment); library(SingleCellExperiment)
source('function/bulk_dat.R')
source('function/scell_dat.R')
source('function/df_match.R')

# Obtain reproducible results.
set.seed(50)
cache.pa <- '~/.cache/shm'

# Mouse brain aSVG.
svg.mus.brain.pa <- system.file("extdata/shinyApp/data", "mus_musculus.brain.svg", package="spatialHeatmap")
svg.mus.brain <- read_svg(svg.mus.brain.pa)
svg.mus.brain <- read_svg('img/mus_musculus.brain_sp1.svg', srsc=TRUE)

# Data paths.
# Male mice strain C57/BL6, 9 weeks old.
sc.pa.mus<- '~/bigdata/single_cell/data/mouse_brain/single_cell'
# Male mice, WT 30-day old.
blk.pa.mus533 <- '~/bigdata/single_cell/data/mouse_brain/PRJNA725533/result' 

# Matching table between bulk and single cells.
# df.match.mus.brain <- df_match_mus533()

# Read data.
# Bulk tissues of mouse brain.
blk.mus.brain <- blk_dat_mus(pa=file.path(blk.pa.mus533, 'countDFeByg.xls'))
dim(blk.mus.brain); blk.mus.brain[1:2, ]

# Bulk tissues are named with aSVG features.
cvt_vecter <- spatialHeatmap:::cvt_vector
colnames(blk.mus.brain) <- unname(cvt_vecter(c('CERE', 'HIPP', 'CERE.CORTEX', 'HYPOTHA'), c('cerebellum', 'hippocampus', 'cerebral.cortex', 'hypothalamus'), colnames(blk.mus.brain)))
blk.mus.brain[1:2, ]

colnames(blk.mus.brain) %in% attribute(svg.mus.brain)[[1]]$feature

blk <- SummarizedExperiment(assays=list(counts=blk.mus.brain), colData=DataFrame(bulk=colnames(blk.mus.brain)))

blk.fil <- filter_data(data=blk, pOA=c(0.2, 5), CV=c(0.2, 100)); dim(blk.fil)

# Single cells of mouse brain.
sc.mus.brain <- sc_dat_mus_brain(sc.pa=file.path(sc.pa.mus, 'GSE147747_expr_raw_counts_table.tsv'), meta.pa=file.path(sc.pa.mus, 'GSE147747_meta_table.tsv'))
dim(sc.mus.brain); sc.mus.brain[1:3, 1:5]
blk.sc <- filter_cell(sce=sc.mus.brain, bulk=blk.fil, gen.rm=NULL, cutoff=1, p.in.cell=0.7, p.in.gen=0.5); blk.sc

bulk <- blk.sc$bulk
save(bulk, file='data/bulk')
load('data/bulk')

# Spatial enrichment.
blk.sub <- sf_var(data=bulk, feature='bulk', ft.sel=colnames(bulk), variable=NULL, var.sel=NULL, com.by='ft')


enr.res <- spatial_enrich(blk.sub, method=c('edgeR'), norm='TMM', log2.fc=1, fdr=0.05, outliers=0)

en.cort <- query_enrich(enr.res, 'cerebral.cortex')
en.hipp <- query_enrich(enr.res, 'hippocampus')
en.hypo <- query_enrich(enr.res, 'hypothalamus')

ovl_enrich(enr.res, type='up', plot='upset', font.size=5, venn.arg=list(), order.by='freq')

up.cort <- subset(rowData(en.cort), type=='up' & total==3)
dim(up.cort); up.cort[1:2, 1:8]
dn.cort <- subset(rowData(en.cort), type=='down' & total==3)
dim(dn.cort); dn.cort[1:2, 1:8]

up.hipp <- subset(rowData(en.hipp), type=='up' & total==3)
dim(up.hipp); up.hipp[1:2, 1:8]
dn.hipp <- subset(rowData(en.hipp), type=='down' & total==3)
dim(dn.hipp); dn.hipp[1:2, 1:8]

up.hypo <- subset(rowData(en.hypo), type=='up' & total==3)
dim(up.hypo); up.hypo[1:2, 1:8]
dn.hypo <- subset(rowData(en.hypo), type=='down' & total==3)
dim(dn.hypo); dn.hypo[1:2, 1:8]

# Co-visualization in manuscript.
# Obtain example data by harsh filtering, and Take overlap genes between bulk and single cells. 
# blk.mus.brain <- filter_data(data=blk.mus.brain, pOA=c(0.3, 6), CV=c(0.55, 100)); dim(blk.mus.brain)

# mus.brain <- filter_cell(sce=sc.mus.brain, bulk=blk.mus.brain, gen.rm=NULL, cutoff=1, p.in.cell=0.7, p.in.gen=0.5); mus.brain

# Joint normalization of bulk and sc data.
blk.sc.nor <- read_cache(cache.pa, 'blk.sc.nor')
if (is.null(blk.sc.nor)) {
  blk.sc.nor <- norm_cell(sce=blk.sc$cell, bulk=blk.sc$bulk, com=FALSE)
  save_cache(dir=cache.pa, overwrite=TRUE, blk.sc.nor)
}

save(blk.sc.nor, file='data/blk.sc.nor')
load('data/blk.sc.nor')

blk.aggr <- aggr_rep(data=blk.sc.nor$bulk, assay.na='logcounts', sam.factor='sample', aggr='mean')

cell <- reduce_dim(blk.sc.nor$cell)
save(cell, file='data/cell')
load('data/cell')

match2cell <- list(cerebral.cortex=c('isocort'), hypothalamus=c('hypotha'), hippocampus='hipp', cerebellum='cere')

dat.sep <- SPHM(svg=svg.mus.brain, bulk=blk.aggr, cell=cell, match=match2cell)
# Tbr1, Mef2c, Rgs4 
gg <- covis(data=dat.sep, ID=rownames(en.hypo)[c(1)], sam.factor='sample', con.factor=NULL, dimred='TSNE', cell.group='sample', tar.bulk=names(match2cell), bar.width=0.09, dim.lgd.nrow=2, dim.lgd.text.size=12, height=0.7, legend.r=0.2, legend.key.size=0.02, legend.text.size=12, legend.nrow=5, col.idp=TRUE, h=0.35, dim.axis.font.size=8, sub.title.size=0) 
ggsave(file="covis_sub.svg", plot=output(gg)$spatial_heatmap, width=10, height=8)
ggsave(file="covis_sub.jpg", plot=output(gg)$spatial_heatmap, width=10, height=8)


# Tbr1, Mef2c, Rgs4 
gg <- covis(data=dat.sep, ID=rownames(en.cort)[c(2)], sam.factor='sample', con.factor=NULL, dimred='TSNE', cell.group='sample', tar.bulk=names(match2cell), bar.width=0.09, dim.lgd.nrow=2, dim.lgd.text.size=12, height=0.7, legend.r=0.2, legend.key.size=0.02, legend.text.size=12, legend.nrow=5, col.idp=TRUE, h=0.35, dim.axis.font.size=8, sub.title.size=0) 
ggsave(file="covis_sub.svg", plot=output(gg)$spatial_heatmap, width=10, height=8)
ggsave(file="covis_sub.jpg", plot=output(gg)$spatial_heatmap, width=10, height=8)


# Slc6a11, Ahi1, Gaa, Ndn, Gprasp2 
gg <- covis(data=dat.sep, ID=rownames(en.hypo)[1], sam.factor='sample', con.factor=NULL, dimred='TSNE', cell.group='sample', tar.bulk=names(match2cell), bar.width=0.09, dim.lgd.nrow=2, dim.lgd.text.size=12, height=0.7, legend.r=0.2, legend.key.size=0.02, legend.text.size=12, legend.nrow=5, col.idp=TRUE, vari.cell=NULL, h=0.35, dim.axis.font.size=8, sub.title.size=12) 
ggsave(file="covis_clas1.svg", plot=output(gg)$spatial_heatmap, width=10, height=8)
ggsave(file="covis_clas1.jpg", plot=output(gg)$spatial_heatmap, width=10, height=8)

ids <- setdiff(rownames(blk.aggr), rownames(en.cort))
# AI413582 
gg <- covis(data=dat.sep, ID=ids[1], sam.factor='sample', con.factor=NULL, dimred='TSNE', cell.group='sample', tar.bulk=names(match2cell), bar.width=0.09, dim.lgd.nrow=2, dim.lgd.text.size=12, height=0.7, legend.r=0.15, legend.key.size=0.02, legend.text.size=12, legend.nrow=5, col.idp=TRUE, vari.cell=NULL, h=0.35, dim.axis.font.size=8, sub.title.size=12) 
ggsave(file="covis_clas2.svg", plot=output(gg)$spatial_heatmap, width=10, height=8)
ggsave(file="covis_clas2.jpg", plot=output(gg)$spatial_heatmap, width=10, height=8)

# Ddn, Cpne6, St6galnac5, Crym, Stum 
ft.trans=NULL
ft.trans='cerebral.cortex'
gg <- covis(data=dat.sep, ID=rownames(en.hipp)[8], sam.factor='sample', con.factor=NULL, dimred='TSNE', cell.group='sample', tar.bulk=names(match2cell), bar.width=0.09, dim.lgd.nrow=2, dim.lgd.text.size=12, height=0.7, legend.r=0.2, legend.key.size=0.02, legend.text.size=12, legend.nrow=5, col.idp=TRUE, vari.cell=NULL, h=0.4, dim.axis.font.size=8, ft.trans=ft.trans) 
ggsave(file="covis_clas.svg", plot=output(gg)$spatial_heatmap, width=10, height=8)

# Cell by group
cell.aggr <- aggr_rep(data=blk.sc.nor$cell, assay.na='logcounts', sam.factor='sample', aggr='mean')

match2blk <- list(isocort=c('cerebral.cortex'), hypotha=c('hypothalamus'))

dat.2blk <- SPHM(svg=svg.mus.brain, bulk=cell.aggr, cell=cell, match=match2blk)

gg <- covis(data=dat.2blk, ID=rownames(en.cort)[c(1)], sam.factor='sample', con.factor=NULL, dimred='TSNE', cell.group='sample', tar.cell=names(match2blk), bar.width=0.13, dim.lgd.nrow=1, dim.lgd.text.size=15, height=0.7, legend.r=1.5, legend.key.size=0.03, legend.text.size=15, legend.nrow=5, col.idp=FALSE, h=0.35, dim.axis.font.size=8, sub.title.size=0, bar.value.size=15) 
ggsave(file="covis_2cell.svg", plot=output(gg)$spatial_heatmap, width=10, height=8)



# Feature by group.

blk.aggr <- aggr_rep(data=blk.sc.nor$bulk, assay.na='logcounts', sam.factor='sample', aggr='mean')

match2cell <- list(cerebral.cortex=c('isocort'), hypothalamus=c('hypotha'))

dat.2cell <- SPHM(svg=svg.mus.brain, bulk=blk.aggr, cell=cell, match=match2cell)

gg <- covis(data=dat.2cell, ID=rownames(en.hypo)[c(1)], sam.factor='sample', con.factor=NULL, dimred='TSNE', cell.group='sample', tar.bulk=names(match2cell), bar.width=0.11, dim.lgd.nrow=1, dim.lgd.text.size=15, height=0.7, legend.r=1.5, legend.key.size=0.03, legend.text.size=15, legend.nrow=5, col.idp=FALSE, h=0.35, dim.axis.font.size=8, sub.title.size=0, bar.value.size=15) 
ggsave(file="covis_2cell.svg", plot=output(gg)$spatial_heatmap, width=10, height=8)

# Fixed by group
gg <- covis(data=dat.2cell, ID=rownames(en.hypo)[c(1)], sam.factor='sample', con.factor=NULL, dimred='TSNE', cell.group='sample', tar.bulk=names(match2cell), bar.width=0.11, dim.lgd.nrow=1, dim.lgd.text.size=15, height=1, legend.r=1.5, legend.key.size=0.02, legend.text.size=15, legend.nrow=5, col.idp=FALSE, h=0.35, dim.axis.font.size=8, sub.title.size=0, bar.value.size=15, profile=FALSE) 

# devtools::install_github('satijalab/seurat-data')
library(Seurat); library(SeuratData);
# InstallData("stxBrain")
set.seed(50)
brain <- LoadData("stxBrain", type = "anterior1")
SpatialFeaturePlot(brain, features = c("Hpca"))

nor.lis <- norm_srsc(cell=brain, assay='Spatial', bulk=bulk)
srt.sc <- nor.lis$cell; blk.sp <- nor.lis$bulk
blk.aggr.sp <- aggr_rep(data=blk.sp, assay.na='logcounts', sam.factor='sample', aggr='mean')

srt.sc <- RunPCA(srt.sc, assay = "SCT", verbose = FALSE)
srt.sc <- RunUMAP(srt.sc, assay = "SCT", dims = 1:5)
srt.sc <- RunTSNE(srt.sc, assay = "SCT", reduction = "pca", dims = 1:5)

srt.sc <- FindNeighbors(srt.sc, reduction = "pca", dims = 1:30)
srt.sc <- FindClusters(srt.sc, verbose = FALSE)
srt.sc$seurat_clusters <- paste0('clus', srt.sc$seurat_clusters)

svg.mus <- read_svg('img/mus_musculus.brain_sp1.svg', srsc=TRUE)
angle(svg.mus)[[1]] <- angle(svg.mus)[[1]] + 90

lis.match <- list(cerebral.cortex=c('clus3', 'clus4', 'clus7', 'clus2'))
dat.srsc <- SPHM(svg=svg.mus, bulk=blk.aggr.sp, cell=srt.sc, match=lis.match)

grp <- unique(srt.sc$seurat_clusters)
sp <- setNames(c(0, 2:14, 19:25, 32:127)[seq_along(grp)], grp)
sp[c('clus2', 'clus3', 'clus4', 'clus7')] <- 15:18 
sp <- setNames(rep(19, length(grp)), grp)

gg <- covis(data=dat.srsc, ID='Efhd2', assay.na='logcounts', dimred='TSNE', cell.group='seurat_clusters', tar.bulk=c('cerebral.cortex'), col.idp=TRUE, bar.width=0.07, dim.lgd.nrow=2, dim.lgd.text.size=8, legend.r=0.1, legend.key.size=0.013, legend.text.size=9, legend.nrow=3, h=0.6, profile=TRUE, ncol=3, vjust=5, dim.lgd.key.size=3, size.r=0.97, dim.axis.font.size=8, size.pt=1.5, shape=sp, lgd.plots.size=c(0.35, 0.25, 0.35), verbose=FALSE)
ggsave(file="covis_srsc.jpg", plot=output(gg)$spatial_heatmap, width=10, height=8)


library(MuSiC); library(Biobase); library(SingleCellExperiment); library(scater)
set.seed(50)
# Bulk data.
blk.eset <- readRDS('GSE50244bulkeset.rds')
blk.de <- exprs(blk.eset); blk.de[1:2, 1:3]
# Single-cell data.
sc.de <- readRDS('EMTABsce_healthy.rds')

# Deconvolution.
res.p <- music_prop(bulk.mtx = blk.de, sc.sce = sc.de, clusters = 'cellType', samples = 'sampleID', select.ct = c('alpha', 'beta', 'delta', 'gamma', 'acinar', 'ductal'), verbose = FALSE)
names(res.de)
prop <- res.p$Est.prop.weighted; prop[1:3, ]

library(reshape2); library(cowplot)

m.prop <- rbind(reshape2::melt(res.p$Est.prop.weighted), reshape2::melt(res.p$Est.prop.allgene))

colnames(m.prop) = c('Sub', 'CellType', 'Prop')
m.prop$CellType = factor(m.prop$CellType, levels = c('alpha', 'beta', 'delta', 'gamma', 'acinar', 'ductal'))
m.prop$Method = factor(rep(c('MuSiC', 'NNLS'), each = 89*6), levels = c('MuSiC', 'NNLS'))
m.prop$HbA1c = rep(blk.eset$hba1c, 2*6)
m.prop = m.prop[!is.na(m.prop$HbA1c), ]
m.prop$Disease = factor(c('Normal', 'T2D')[(m.prop$HbA1c > 6.5)+1], levels = c('Normal', 'T2D'))
m.prop$D = (m.prop$Disease == 'T2D')/5
m.prop = rbind(subset(m.prop, Disease == 'Normal'), subset(m.prop, Disease != 'Normal'))
m.prop[1:3, ]

# Sub1: male, normal

blk.sel <- subset(data.frame(blk.de), Sub1 >= 5)[, 'Sub1', drop=FALSE]

res.de <- lapply(X=seq_len(ncol(blk.sel)), bulk=blk.sel, FUN=function(i, bulk) {
  lis0 <- lapply(seq_len(nrow(bulk)), function(r) { 
    round(data.frame(bulk)[r, i, drop=TRUE] * prop[i, ])
  }); dat <- do.call('rbind', lis0)
  rownames(dat) <- rownames(bulk); dat
})

res.sc <- res.de[[1]]
res.sc <- SingleCellExperiment(list(counts=res.sc)); dim(res.sc) 

blk.sc.de <- filter_cell(sce=res.sc, bulk=blk.sel, gen.rm=NULL, cutoff=1, p.in.cell=0.1, p.in.gen=0.5); blk.sc.de

blk.se.de.nor <- norm_data(data=cbind(blk.sc.de$bulk, blk.sc.de$cell), norm.fun='CNF', log2.trans=TRUE)
assayNames(blk.se.de.nor) <- 'logcounts'

blk.de.nor <- blk.se.de.nor[, seq_len(ncol(blk.sc.de$bulk))]
sc.de.nor <- blk.se.de.nor[, seq_len(ncol(blk.sc.de$cell)) + ncol(blk.de.nor)]

# Deconvolution.
colnames(blk.de.nor) <- 'pancreas.islet'
sc.de.nor$bulk <- 'pancreas.islet'
sc.de.nor$cell.type <- colnames(sc.de.nor)
sc.de.nor$prop <- round(prop['Sub1', sc.de.nor$cell.type], 3)
# sc.de.nor <- reduce_dim(sc.de.nor)
sc.de.dim <- runTSNE(sc.de.nor)

islet <- read_svg('pancreas.svg')
dat.decon <- SPHM(svg=islet, bulk=blk.de.nor, cell=sc.de.dim)
cell.all <- unique(sc.de.dim$cell.type)
sp <- setNames(c(15, 17:19, 12, 8), cell.all)
sp <- setNames(rep(8, length(cell.all)), cell.all)

# AP006222.2, RP11-206L10.9
gg <- covis(data=dat.decon, ID=rownames(blk.de.nor)[4], dimred='TSNE', tar.cell=unique(sc.de.dim$cell.type), cell.group='cell.type', dim.lgd.text.size=10, dim.lgd.key.size=5, dim.lgd.nrow=4, bar.width=0.1, legend.nrow=4, h=0.4, col.idp=T, legend.r=0.2, profile=T, decon=T, size.pt=5, size.lab.pt=5, vjust.lab.pt=1.6, hjust.lab.pt=0.5, shape=sp)
ggsave(file="covis_decon.svg", plot=output(gg)$spatial_heatmap, width=10, height=8)


