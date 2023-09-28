
## DataFrame with 3 rows and 3 columns
library(spatialHeatmap); library(SummarizedExperiment); library(ExpressionAtlas); library(BiocParallel)
library(biomaRt) # v2.54.0 
library(dendextend); library(enrichplot); library(ggplot2)
library(clusterProfiler) # v4.6.0

cache.pa <- '~/.cache/shm'

all.mus <- read_cache(cache.pa, 'all.mus') # Retrieve data from cache.
if (is.null(all.mus)) { # Save downloaded data to cache if it is not cached.
  all.mus <- searchAtlasExperiments(properties="kidney", species="Mus musculus")
  save_cache(dir=cache.pa, overwrite=TRUE, all.mus)
}

all.mus[7, ]
rse.mus <- read_cache(cache.pa, 'rse.mus') # Read data from cache.
if (is.null(rse.mus)) { # Save downloaded data to cache if it is not cached.
  rse.mus <- getAtlasData('E-MTAB-2801')[[1]][[1]]
  save_cache(dir=cache.pa, overwrite=TRUE, rse.mus)
}

colData(rse.mus)[1:3, ]

mus.tar <- system.file('extdata/shinyApp/data/target_mouse.txt', package='spatialHeatmap')
target.mus <- read.table(mus.tar, header=TRUE, row.names=1, sep='\t') # Importing
colData(rse.mus) <- DataFrame(target.mus) # Loading
target.mus[1:3, ]

data.sub.mus <- sf_var(data=rse.mus, feature='organism_part', ft.sel=c('brain', 'lung', 'colon', 'kidney', 'liver'), variable='strain', var.sel=c('DBA.2J', 'C57BL.6', 'CD1'), com.by='ft')

data.sub.mus.fil <- filter_data(data=data.sub.mus, pOA=c(0.5, 15), CV=c(0.8, 100), verbose=FALSE)


cvt_mart <- function(input, biomart, dataset, mirror, filters, attributes, target, use.cache=TRUE, cache.na, cache.dir='~/.cache/shm') {
  ann <- NULL
  if (use.cache==TRUE) ann <- read_cache(cache.dir, cache.na)
  if (is.null(ann)|use.cache==FALSE) { 
    ensem <- useEnsembl(biomart=biomart, dataset = dataset, mirror=mirror) 
    ann <- getBM(filters=filters, values=rownames(input), attributes=attributes, mart=ensem) 
    ann <- subset(ann, !duplicated(get(filters))) 
    # Keep all from ids.
    tgt <- ann[, target]; idx <- tgt=='' | duplicated(tgt)
    ann$target.id <- make.names(ann[, target])
    # If no target ids, from ids are maintained.
    ann$target.id[idx] <- ann[, filters][idx] 
    rownames(ann) <- ann$target.id <- make.names(ann$target.id) 
    save_cache(dir=cache.pa, overwrite=TRUE, obj=ann, na=cache.na)   
  }
  # Maintain ids not in biomaRt.
  inter <- intersect(rownames(input), ann[, filters]) 
  data.dif <- input[setdiff(rownames(input), inter)]
  # Make a rowData from ids not in biomaRt.
  rdat.dif <- ann[seq_len(nrow(data.dif)), ]
  rownames(rdat.dif) <- rownames(data.dif)
  rdat.dif[, seq_len(ncol(rdat.dif))] <- NA
  rowData(data.dif) <- rdat.dif
  # Assign target ids to ids present in biomaRt.
  input <- input[inter, ]
  ann <- subset(ann, get(filters) %in% inter) 
  input <- input[order(rownames(input)), ] 
  ann <- ann[order(ann[, filters]), ]
  all(rownames(input)==ann[, filters])
  rownames(input) <- ann$target.id
  rowData(input) <- ann
  # Maintain all ids from input.
  input <- rbind(input, data.dif); input
}

data.sub.mus.fil <- cvt_mart(input=data.sub.mus.fil, biomart='ensembl', dataset="mmusculus_gene_ensembl", mirror='www', filters='ensembl_gene_id', attributes=c("ensembl_gene_id", "uniprot_gn_symbol", "entrezgene_id", "go_id", "description"), target='uniprot_gn_symbol', use.cache=TRUE, cache.na='ann.mus', cache.dir='~/.cache/shm')
ann.mus <- rowData(data.sub.mus.fil)

ann.mus <- NULL
# ann.mus <- read_cache(cache.pa, 'ann.mus') 
if (is.null(ann.mus)) { 
  ensem <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror='www') 
  ann.mus <- getBM(filters='ensembl_gene_id', values=rownames(data.sub.mus.fil), attributes=c("ensembl_gene_id", "uniprot_gn_symbol", "entrezgene_id", "go_id", "description"), mart=ensem) 
  ann.mus <- subset(ann.mus, !duplicated(ensembl_gene_id)) 
  # Gene symbols.
  uni.sym <- ann.mus$uniprot_gn_symbol 
  idx <- uni.sym=='' | duplicated(uni.sym)
  ann.mus$symbol <- make.names(ann.mus$uniprot_gn_symbol)
  ann.mus$symbol[idx] <- ann.mus$ensembl_gene_id[idx] 
  rownames(ann.mus) <- ann.mus$symbol <- make.names(ann.mus$symbol) 
  save_cache(dir=cache.pa, overwrite=TRUE, ann.mus)   
} 

inter <- intersect(rownames(data.sub.mus.fil), ann.mus$ensembl_gene_id) 
data.dif <- data.sub.mus.fil[setdiff(rownames(data.sub.mus.fil), inter)] 
data.sub.mus.fil <- data.sub.mus.fil[inter, ]
ann.mus <- subset(ann.mus, ensembl_gene_id %in% inter) 
data.sub.mus.fil <- data.sub.mus.fil[order(rownames(data.sub.mus.fil)), ] 
ann.mus <- ann.mus[order(ann.mus$ensembl_gene_id), ] 
rownames(data.sub.mus.fil) <- ann.mus$symbol
data.sub.mus.fil <- rbind(data.sub.mus.fil, data.dif)  


# deg.lis.mus <- spatial_enrich(data.sub.mus.fil, methods=c('edgeR'), norm='TMM', log2.fc=1, fdr=0.05, aggr='mean', log2.trans.aggr=TRUE)

enr.res <- spatial_enrich(data.sub.mus.fil, method=c('edgeR'), norm='TMM', log2.fc=1, fdr=0.05, outliers=1)


en.liver <- query_enrich(enr.res, 'liver')
en.brain <- query_enrich(enr.res, 'brain')

ovl_enrich(enr.res, type='up', plot='upset', font.size=5, venn.arg=list(), order.by='freq')

up.liver <- subset(rowData(en.liver), type=='up' & total==4)
up.liver[1:5, 1:8]
up.brain <- subset(rowData(en.brain), type=='up' & total==4)
up.brain[1:5, 1:8]

dn.liver <- subset(rowData(en.liver), type=='down' & total==4)
dn.liver[1:5, 1:8]
dn.brain <- subset(rowData(en.brain), type=='down' & total==4)
dn.brain[1:5, 1:8]

# Apoa1
gene <- 'Fxn'
en.liver <- aggr_rep(en.liver, sam.factor='com.by', con.factor=NULL, aggr='mean')
gene <- cl[28:36]
en.brain <- aggr_rep(en.brain, sam.factor='com.by', con.factor=NULL, aggr='mean')

svg.mus.pa <- system.file("extdata/shinyApp/data", "mus_musculus.male.svg", package="spatialHeatmap")
svg.mus.pa <- 'mus_musculus.male.svg'
svg.mus <- read_svg(svg.mus.pa)
cord <- svg.mus[, 2][[1]]; cord[1, 4] <- 0
attribute(svg.mus)[[1]] <- cord


gene <- c('Itih1', 'Mgat3') 

dat.enrich <- SPHM(svg=svg.mus, bulk=en.liver)
shm(data=dat.enrich, ID=gene, legend.r=1, legend.nrow=3, sub.title.size=6.1, ncol=6, bar.width=0.09, lay.shm='con')


graph_line(assay(en.liver)[gene, ], lgd.pos='right', x.title = "Samples", y.title = "Assay values", lgd.guide = guides(color = guide_legend(nrow = 2, byrow = TRUE, title = NULL)))


graph_line(assay(en.brain)[rownames(up.brain)[1:50], ], lgd.pos='right', x.title = "Samples", y.title = "Assay values", lgd.guide = guides(color = FALSE))


mus.nor.aggr <- aggr_rep(enr.res$data, sam.factor='organism_part', con.factor=NULL, aggr='mean')

gene <- rownames(up.liver)[2]
gene <- rownames(up.brain)[200]
gene <- 'Fgd6' # liver, hclus, h = 3.210751
gene <- 'Serpina1b' # liver, network ds=3
gene <- 'Grik3'; h=7.149385 # brain, hclus, h=7.149385
gene <- 'Ugp2'; h=3.685629 # liver, hclus
sub.mat <- submatrix(data=mus.nor.aggr, ID=gene, p=0.15); dim(sub.mat)
png('del')
mhm.res <- matrix_hm(ID=gene, data=sub.mat, scale='row', angleCol=80, angleRow=35, cexRow=0.8, cexCol=0.8, margin=c(10, 6), static=TRUE, arg.lis1=list(offsetRow=0.01, offsetCol=0.01))
dev.off()

ply <- matrix_hm(ID=gene, data=sub.mat, scale='row', col=c('yellow', 'red'), angleCol=80, angleRow=35, cexRow=0.8, cexCol=0.8, margin=c(10, 6), static=FALSE, arg.lis1=list(offsetRow=0.01, offsetCol=0.01))

# Cut row dendrogram.
dendro.x <- mhm.res$rowDendrogram
clus <- cutree(dendro.x, h = h)[order.dendrogram(dendro.x)]
sort(unique(clus))
w <- which(names(clus)==gene) 
clus[w]
cl <- names(clus[clus==clus[w]]); cl

cut_dendro <- function(dendro, h, target) {
  clus <- cutree(dendro, h = h)[order.dendrogram(dendro)]
  w <- which(names(clus)==target)
  cl <- names(clus[clus==clus[w]]); sort(cl)
}

cl <- cut_dendro(mhm.res$rowDendrogram, h, gene); cl


df.entr <- subset(ann.mus, uniprot_gn_symbol %in% cl)
rownames(df.entr) <- seq_len(nrow(df.entr))
write.table(df.entr[, c('ensembl_gene_id', 'uniprot_gn_symbol', 'entrezgene_id')], 'cluster_hc.xlsx', col.names=T, row.names=T, sep='\t')

entriz <- unique(df.entr$entrezgene_id); entriz

adj.mod <- adj_mod(data=sub.mat)

adj.mod[['adj']][1:3, 1:3]
adj.mod[['mod']][1:3, ] 
png('del')
net <- network(ID=gene, data=sub.mat, adj.mod=adj.mod, ds=3, adj.min=0, vertex.label.cex=2, vertex.cex=7, static=TRUE)
dev.off()

df.entr.net <- subset(ann.mus, uniprot_gn_symbol %in% rownames(net))
rownames(df.entr.net) <- seq_len(nrow(df.entr.net))
write.table(df.entr.net[, c('ensembl_gene_id', 'uniprot_gn_symbol', 'entrezgene_id')], 'module.xlsx', col.names=T, row.names=T, sep='\t')

entriz <- unique(df.entr.net$entrezgene_id); entriz

kk <- enrichKEGG(gene = entriz, keyType='kegg', organism = 'mmu', pvalueCutoff = 0.05, pAdjustMethod='BH', qvalueCutoff=0.05)
df.kk <- as.data.frame(kk)
barplot(kk, showCategory=20)
browseKEGG(kk, 'mmu04146')

id.kk <- unique(unlist(strsplit(df.kk$geneID[1], '/')))
id.kk <- unique(unlist(strsplit(df.kk$geneID, '/')))
id.kk <- subset(ann.mus, entrezgene_id %in% id.kk)$symbol; id.kk


go <- enrichGO(gene = entriz, keyType = "ENTREZID", OrgDb = 'org.Mm.eg.db', ont = 'ALL', pAdjustMethod = 'BH', pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)
df.go <- as.data.frame(go)
barplot(go, showCategory=5)

g.new <- setdiff(cl, c(id.kk, unique(unlist(strsplit(df.go$geneID, '/'))))); g.new

ENSMUSG00000048388, ENSMUSG00000072769,  


output <- spatialHeatmap::output
gene <- sort(cl)[-c(11:13)]
gene <- sort(cl)
gene <- setdiff(c(g.new, gene), 'Caskin1')
gg <- shm(data=dat.enrich, ID=gene, legend.r=1, legend.nrow=3, ncol=6, lay.shm='con', bar.width=0.12, bar.value.size=10, sub.title.size=17, legend.plot.title.size=2, legend.text.size=15, line.width=0, line.color='grey30')  
ggsave(file="spatial_enrich_use.svg", plot=output(gg)$spatial_heatmap, width=10, height=8)




genes <- rownames(up.liver)[1:3]
genes <- rownames(up.brain)

# Selected: Fgd6.
hc.liver <- query_hclus(gene=genes, data=mus.nor.aggr, p=0.15, ann=ann.mus, from.id='uniprot_gn_symbol', entrez.id='entrezgene_id', clus.min=10, clus.max=50, h.inter=0.5, organism='mmu', min.count=5, OrgDb='org.Mm.eg.db')

hc.br <- query_hclus(gene=genes, data=mus.nor.aggr, p=0.15, ann=ann.mus, from.id='uniprot_gn_symbol', entrez.id='entrezgene_id', clus.min=10, clus.max=50, h.inter=0.5, organism='mmu', min.count=5, OrgDb='org.Mm.eg.db')

load('hc.br')
# Serac1, Grik3: 5, h=7.149385, Nckap1, Ank3
for (i in seq_along(hc.br)) {
  lis0 <- hc.br[[i]]
  for (j in seq_along(lis0)) {
    l0 <- lis0[[j]]
    if (!is.null(l0$kegg)) { 
      print('kegg'); print(head(l0$kegg)) 
      message(names(hc.liver)[i], ' ', length(lis0[[j]]$cluster), ' ', lis0[[j]]$height) 
    }
    if (!is.null(l0$go)) { 
      print('go'); print(head(l0$go))
      message(names(hc.liver)[i], ' ', length(lis0[[j]]$cluster), ' ', lis0[[j]]$height) 
    }
  }
}

load('hc.liver')
# Ugp2: 2, h=3.685629, Eif4ebp2: 4, h=5.185629, Rnase4: 1, h=7.71075105372941, Man1a: 3, h=4.18562906337543 
for (i in seq_along(hc.liver)) {
  lis0 <- hc.liver[[i]]
  for (j in seq_along(lis0)) {
    l0 <- lis0[[j]]
    if (!is.null(l0$kegg)) { 
      print('kegg'); print(head(l0$kegg))
      message(names(hc.liver)[i], ' ', length(lis0[[j]]$cluster), ' ', lis0[[j]]$height) 
    }
    if (!is.null(l0$go)) { 
      print('go'); print(head(l0$go))
      message(names(hc.liver)[i], ' ', length(lis0[[j]]$cluster), ' ', lis0[[j]]$height) 
    }
  }
}

load('liver.net')
# ds=2, Slc25a13, size: 20; 1 Asl 12; 3 Brawnin 9; 3 Serpina1b 9; 2 Serpina1b 9  
for (i in seq_along(liver.net)) {
  lis0 <- liver.net[[i]]
  for (j in seq_along(lis0)) {
    l0 <- lis0[[j]]
    if (!is.null(l0$kegg)) { 
      print('kegg'); print(head(l0$kegg))
      message(names(liver.net)[i], ' ', names(lis0)[j], ' ',  length(lis0[[j]]$module)) 
    }
    if (!is.null(l0$go)) { 
      print('go'); print(head(l0$go))
      message(names(liver.net)[i], ' ', names(lis0)[j], ' ',  length(lis0[[j]]$module)) 
    }
  }
}


entriz <- unique(subset(ann.mus, uniprot_gn_symbol %in% rownames(net))$entrezgene_id); entriz

keg <- en_kegg(gene=entriz, keyType='kegg', organism='mmu', pvalueCutoff=0.05, qvalueCutoff=0.05, min.count=5)
  
go <- en_go(gene=entriz, OrgDb='org.Mm.eg.db', min.count=5)

ego <- enrichGO(gene = entriz, keyType = "ENTREZID", OrgDb = 'org.Mm.eg.db', ont = "ALL", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)
head(ego)

genes <- rownames(up.brain)[200:205]

# p=0.15, ds=3
nt3 <- nt; save(nt3, file='nt3')

nt <- query_net(gene=genes, data=mus.nor.aggr, p=0.15, ds=3, ann=ann.mus, from.id='uniprot_gn_symbol', entrez.id='entrezgene_id', organism='mmu', OrgDb='org.Mm.eg.db', min.count=2)

# Selected: Fxn
load('nt3')

genes <- rownames(up.liver)
genes <- rownames(up.brain)

br.liver <- list('0'=NULL, '1'=NULL, '2'=NULL, '3'=NULL)
for (i in as.character(0:3)) {
  nt0 <- query_net(gene=genes, data=mus.nor.aggr, p=0.15, ds=i, ann=ann.mus, from.id='uniprot_gn_symbol', entrez.id='entrezgene_id', organism='mmu', OrgDb='org.Mm.eg.db', min.count=3)
 br.liver[[i]] <- nt0
}

res.br <- bplapply(as.character(0:3), BPPARAM=MulticoreParam(workers=4, RNGseed=50, log=TRUE, logdir='./'), FUN=function(i, gene, data, ann) {
  query_net(gene=gene, data=data, p=0.15, ds=i, ann=ann, from.id='uniprot_gn_symbol', entrez.id='entrezgene_id', organism='mmu', OrgDb='org.Mm.eg.db', min.count=5) }, gene=genes, data=mus.nor.aggr, ann=ann.mus)

# Functions
query_hclus <- function(gene, data, p=0.15, ann, from.id, entrez.id, clus.min=10, clus.max=50, h.inter=0.5, organism, pvalueCutoff=0.05, pAdjustMethod='BH', qvalueCutoff=0.05, min.count=5, OrgDb, ont='ALL', readable=TRUE) {
  lis.all <- NULL
  for (g in gene) {
    sub.mat <- submatrix(data=data, ID=g, p=p)
    tmp <- paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/mhm.png')
    png(tmp)
    mhm.res <- matrix_hm(ID=g, data=sub.mat, scale='row', angleCol=80, angleRow=35, cexRow=0.8, cexCol=0.8, margin=c(10, 6), static=TRUE, arg.lis1=list(offsetRow=0.01, offsetCol=0.01))
    dev.off()
    dendro <- mhm.res$rowDendrogram
    h.max <- max(get_branches_heights(dendro))
    h.min <- min(get_branches_heights(dendro))
    # Cut the tree at different heights to get clusters for each query gene.
    l0 <- query_hclus0(gene=g, ann=ann, from.id=from.id, entrez.id=entrez.id, dendro=dendro, clus.min=clus.min, clus.max=clus.max, h.min=h.min, h.max=h.max, h.inter=h.inter, organism=organism, pvalueCutoff=pvalueCutoff,  pAdjustMethod=pAdjustMethod, qvalueCutoff=qvalueCutoff, min.count=min.count, OrgDb=OrgDb, ont=ont, readable=readable)
     if (!is.null(l0)) {
       l0 <- list(l0); names(l0) <- g; lis.all <- c(lis.all, l0)
     }
  }; lis.all
}


query_hclus0 <- function(gene, ann, from.id, entrez.id, dendro, clus.min, clus.max, h.min, h.max, h.inter, organism, pvalueCutoff=pvalueCutoff, pAdjustMethod=pAdjustMethod, qvalueCutoff=qvalueCutoff, min.count, OrgDb=OrgDb, ont=ont, readable=readable) {
  # All candidate heights.
  h <- seq(h.min, h.max, by=h.inter); lis <- NULL
  for (h0 in h ) { # Cut tree at each height.
    clus <- cutree(dendro, h = h0)[order.dendrogram(dendro)]
    w <- which(names(clus)==gene) 
    clus0 <- sort(names(clus[clus==clus[w]])); len <- length(clus0)
    # Check if the cluster is within the size range.
    if (len >= clus.min & len <= clus.max) {
      # Check if the cluster already exists.
      same <- lapply(seq_along(lis), function(x) {
        c1 <- lis[[x]]$cluster
        if (length(c1)!=len) return(FALSE)
        if (all(c1 %in% clus0)) return(TRUE)
      })
      if (sum(unlist(same))==0) {
        # KEGG enrichment.
        entr.id <- unique(ann[, entrez.id][ann[, from.id] %in% clus0])
        df.keg <- en_kegg(gene=entr.id, keyType='kegg', organism=organism, pvalueCutoff=pvalueCutoff,  pAdjustMethod=pAdjustMethod, qvalueCutoff=qvalueCutoff, min.count=min.count)
        # GO enrichment.
        df.go <- en_go(gene=entr.id, keyType="ENTREZID", OrgDb=OrgDb, ont=ont, pAdjustMethod=pAdjustMethod, pvalueCutoff=pvalueCutoff, qvalueCutoff=qvalueCutoff, readable=readable, min.count=min.count)
        if (nrow(df.keg)==0 & nrow(df.go)==0) next
        lis0 <- list(height=h0, cluster=clus0)
        if (nrow(df.keg)>0) {
          if (max(df.keg$Count) >= min.count) lis0$kegg <- df.keg
        }
        if (nrow(df.go)>0) {
          if (max(df.go$Count) >= min.count) lis0$go <- df.go
        }
        if (is.null(lis0$kegg) & is.null(lis0$go)) next
        print(gene); lis <- c(lis, list(lis0))
      }
    }
  }
  if (!is.null(lis)) names(lis) <- as.character(seq_along(lis))
  lis
}

query_net <- function(gene, data, p=0.15, ds=0, ann, from.id='uniprot_gn_symbol', entrez.id='entrezgene_id', organism, pvalueCutoff = 0.05, qvalueCutoff  = 0.05, readable = TRUE, OrgDb, ont = "ALL", pAdjustMethod = "BH", min.count=5) {
  tmp <- paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/net.png')
  # Identify modules for each query gene.
  lis <- NULL; for (g in gene) { 
    sub.mat <- submatrix(data=data, ID=g, p=p)
    adj.mod <- adj_mod(data=sub.mat, ds=0:3)
    # Get the module containing the query gene.
    png(tmp); net <- network(ID=g, data=sub.mat, adj.mod=adj.mod, ds=ds, adj.min=0.9, static=TRUE); dev.off()
    # Functional enrichments for the returned module. 
    gens <- rownames(net); lis0 <- list(module=gens) 
    entr.id <- unique(subset(ann.mus, uniprot_gn_symbol %in% gens)$entrezgene_id)
    df.keg <- en_kegg(gene=entr.id, keyType='kegg', organism=organism, pvalueCutoff=pvalueCutoff,  pAdjustMethod=pAdjustMethod, qvalueCutoff=qvalueCutoff, min.count=min.count)
    df.go <- en_go(gene=entr.id, keyType="ENTREZID", OrgDb=OrgDb, ont=ont, pAdjustMethod=pAdjustMethod, pvalueCutoff=pvalueCutoff, qvalueCutoff=qvalueCutoff, readable=readable, min.count=min.count)
    if (nrow(df.keg)==0 & nrow(df.go)==0) next
    if (nrow(df.keg)>0) {
      if (max(df.keg$Count) >= min.count) lis0$kegg <- df.keg
    }
    if (nrow(df.go)>0) {
      if (max(df.go$Count) >= min.count) lis0$go <- df.go
    }
    if (is.null(lis0$kegg) & is.null(lis0$go)) next
    print(g)
    # Put enrichment results of each query gene in a nested list. 
    lis0 <- list(lis0); names(lis0) <- g; lis <- c(lis, lis0)
  }; lis
}


en_go <- function(gene, keyType="ENTREZID", OrgDb, ont='ALL', pAdjustMethod='BH', pvalueCutoff=0.05, qvalueCutoff=0.05, readable=TRUE, min.count=5) {
  go <- enrichGO(gene = gene, keyType = keyType, OrgDb = OrgDb, ont = ont, pAdjustMethod = pAdjustMethod, pvalueCutoff  = pvalueCutoff, qvalueCutoff  = qvalueCutoff, readable = readable)
  df.go <- as.data.frame(go)
  if (nrow(df.go)>0) {
    if (max(df.go$Count) >= min.count) { 
      print(head(df.go[, c('Description', 'GeneRatio', 'BgRatio', 'p.adjust')]))
    }
  }; return(df.go)
}


en_kegg <- function(gene, keyType='kegg', organism, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', min.count=5) {
  kk <- enrichKEGG(gene = gene, keyType=keyType, organism = organism, pvalueCutoff = pvalueCutoff, pAdjustMethod=pAdjustMethod, qvalueCutoff=qvalueCutoff)
  df.kk <- as.data.frame(kk)
  if (nrow(df.kk)>0) {
    if (max(df.kk$Count) >= min.count) { 
      print(head(df.kk[, c('Description', 'GeneRatio', 'BgRatio', 'p.adjust')]))
    }
  }; return(df.kk)
}

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
