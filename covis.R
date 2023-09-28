
## Example data for annotation labels. Cells of selected ages are combined.
library(spatialHeatmap); library(scRNAseq)
sce.mus <- MarquesBrainData()

# Edit colData.
cdat <- colData(sce.mus)[, -1]; colnames(cdat) <- make.names(colnames(cdat))
cdat$source_name <- make.names(cdat$source_name)
colnames(cdat)[colnames(cdat)=='source_name'] <- 'label'
# "expVar" is a reserved colname to indicate experiment variables.
colnames(cdat)[colnames(cdat)=='treatment'] <- 'expVar'
rownames(cdat) <- make.names(rownames(cdat))

cdat <- edit_tar(cdat, 'expVar', '6hr post acute stress', '6h.post.stress')
cdat <- edit_tar(cdat, 'expVar', 'No', 'control')
cdat[1:2, ]
colData(sce.mus) <- cdat

t1 <- names(table(colData(subset(sce.mus, , label=='cortex.S1'))$age))
t2 <- names(table(colData(subset(sce.mus, , label=='hippocampus.CA1'))$age))
t3 <- names(table(colData(subset(sce.mus, , label=='hypothalamus'))$age))

age.sub <- Reduce(intersect, list(t1, t2, t3))
sce.mus.sub <- subset(sce.mus, , age %in% age.sub)

# Filter cells and genes to obtain example data.
# Cells across selected ages.
sce.mus <- filter_cell(sce=sce.mus.sub, cutoff=1, p.in.cell=0.5, p.in.gen=0.5); sce.mus

sce.mus <- subset(sce.mus, , expVar=='control')
# Save the example data.
saveRDS(sce.mus, file='./cell_mouse_brain.rds')

set.seed(10)
sce.dimred.quick <- process_cell_meta(sce.mus, qc.metric=list(subsets=list(Mt=rowData(sce.mus)$featureType=='mito'), threshold=1))
colData(sce.dimred.quick)[1:3, 1:2]

sce.dimred.g1 <- subset(sce.dimred.quick, , age %in% c("p22", "p23", "p24"))
sce.dimred.g2 <- subset(sce.dimred.quick, , age %in% c("p27", "p28"))

sce.aggr.g1 <- aggr_rep(sce.dimred.g1, assay.na='logcounts', sam.factor='label', aggr='mean')
sce.aggr.g2 <- aggr_rep(sce.dimred.g2, assay.na='logcounts', sam.factor='label', aggr='mean')

# svg.mus.brain.pa <- system.file("extdata/shinyApp/data", "mus_musculus.brain.svg", package="spatialHeatmap")
svg.mus.brain.pa <- 'mus_musculus.brain.svg'
svg.mus.brain <- read_svg(svg.mus.brain.pa)

tail(attribute(svg.mus.brain)[[1]])[, 1:4]
lis.match.quick <- list(cortex.S1=c('somatosensor.cortex'), hippocampus.CA1='hippocampus.ca1', hypothalamus='hypothalamus')

dat1 <- SHM(svg=svg.mus.brain, bulk=sce.aggr.g1, cell=sce.dimred.g1, match=lis.match.quick)
dat2 <- SHM(svg=svg.mus.brain, bulk=sce.aggr.g2, cell=sce.dimred.g2, match=lis.match.quick)

shm.res.quick <- covis(data=dat1, ID=c('Apod'), dimred='TSNE', cell.group='label', tar.cell=names(lis.match.quick), assay.na='logcounts', bar.width=0.11, dim.lgd.nrow=1, height=0.7, legend.r=1.5, legend.key.size=0.02, legend.text.size=12, legend.nrow=3, profile=T) 

df <- assay(sce.aggr.g1)[, c('hypothalamus', 'cortex.S1', 'hippocampus.CA1')]
w <- which(df[, 3]>df[, 1] & df[, 3]>df[, 2])

id <- rownames(sce.aggr.g1)[w[7]]
id <- 'Phlpp1'

gg <- covis(data=dat1, ID=id, dimred='TSNE', cell.group='label', tar.cell=names(lis.match.quick), dim.lgd.nrow=1, height=0.7, legend.r=1.5, legend.key.size=0.02, legend.nrow=3, bar.width=0.11, bar.value.size=17, sub.title.size=17, legend.plot.title.size=17, legend.text.size=17, dim.lgd.text.size=7, dim.lgd.key.size=6, profile=FALSE)
library(ggplot2)
ggsave(file="covisualize_ann_cell2bulk2.svg", plot=output(gg)$spatial_heatmap, width=10, height=4)




## Example data for annotation labels. Cells across all ages are combined.
library(spatialHeatmap); library(scRNAseq)
sce.mus <- MarquesBrainData()
# Filter cells and genes to obtain example data.
sce.mus <- filter_cell(sce=sce.mus, cutoff=1, p.in.cell=0.5, p.in.gen=0.5); sce.mus

# Edit colData.
cdat <- colData(sce.mus)[, -1]; colnames(cdat) <- make.names(colnames(cdat))
cdat$source_name <- make.names(cdat$source_name)
colnames(cdat)[colnames(cdat)=='source_name'] <- 'label'
# "expVar" is a reserved colname to indicate experiment variables.
colnames(cdat)[colnames(cdat)=='treatment'] <- 'expVar'
rownames(cdat) <- make.names(rownames(cdat))

cdat <- edit_tar(cdat, 'expVar', '6hr post acute stress', '6h.post.stress')
cdat <- edit_tar(cdat, 'expVar', 'No', 'control')
cdat[1:2, ]
colData(sce.mus) <- cdat
sce.mus <- subset(sce.mus, , expVar=='control')
# Save the example data.
saveRDS(sce.mus, file='./cell_mouse_brain.rds')

set.seed(10)
sce.dimred.quick <- process_cell_meta(sce.mus, qc.metric=list(subsets=list(Mt=rowData(sce.mus)$featureType=='mito'), threshold=1))
colData(sce.dimred.quick)[1:3, 1:2]

sce.aggr.quick <- aggr_rep(sce.dimred.quick, assay.na='logcounts', sam.factor='label', aggr='mean')

svg.mus.brain.pa <- system.file("extdata/shinyApp/data", "mus_musculus.brain.svg", package="spatialHeatmap")
svg.mus.brain <- read_svg(svg.mus.brain.pa)

tail(attribute(svg.mus.brain)[[1]])[, 1:4]
lis.match.quick <- list(cortex.S1=c('cerebral.cortex'), hippocampus.CA1='hippocampus', hypothalamus='hypothalamus')

dat.quick <- SHM(svg=svg.mus.brain, bulk=sce.aggr.quick, cell=sce.dimred.quick, match=lis.match.quick)

shm.res.quick <- covis(data=dat.quick, ID=c('Apod'), dimred='TSNE', cell.group='label', tar.cell=names(lis.match.quick), assay.na='logcounts', bar.width=0.11, dim.lgd.nrow=1, height=0.7, legend.r=1.5, legend.key.size=0.02, legend.text.size=12, legend.nrow=3, profile=T) 

df <- assay(sce.aggr.quick)[, c('hypothalamus', 'cortex.S1', 'hippocampus.CA1')]
w <- which(df[, 3]>df[, 1] & df[, 3]>df[, 2])

id <- rownames(sce.aggr.quick)[w[7]]
id <- 'Cpm'

gg <- covis(data=dat.quick, ID=id, dimred='TSNE', cell.group='label',  tar.cell=names(lis.match.quick), dim.lgd.nrow=1, height=0.7, legend.r=1.5, legend.key.size=0.02, legend.nrow=3, bar.width=0.11, bar.value.size=17, sub.title.size=17, legend.plot.title.size=17, legend.text.size=17, dim.lgd.text.size=7, dim.lgd.key.size=6)
library(ggplot2)
ggsave(file="covisualize_ann_cell2bulk.svg", plot=output(gg)$spatial_heatmap, width=10, height=4)



