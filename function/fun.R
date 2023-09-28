# Convert Ensembl ids to symbols using biomaRt.
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

# For each of given genes, subset the data matrix, perform hierarchical clustering, cut the row dendrogram tree at different heights to get clusters, and perfrom functional enrichment. The results of each gene are returned in a nested list.  
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

# Cut the row dendrogram tree at different heights to get clusters for a query gene, and perfrom functional enrichment.
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

# For each of given genes, subset the data matrix, perform network analysis, identify the cluster containing the gene under consideration, and perfrom functional enrichment. The results of each gene are returned in a nested list.  
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

# GO enrichment.
en_go <- function(gene, keyType="ENTREZID", OrgDb, ont='ALL', pAdjustMethod='BH', pvalueCutoff=0.05, qvalueCutoff=0.05, readable=TRUE, min.count=5) {
  go <- enrichGO(gene = gene, keyType = keyType, OrgDb = OrgDb, ont = ont, pAdjustMethod = pAdjustMethod, pvalueCutoff  = pvalueCutoff, qvalueCutoff  = qvalueCutoff, readable = readable)
  df.go <- as.data.frame(go)
  if (nrow(df.go)>0) {
    if (max(df.go$Count) >= min.count) { 
      print(head(df.go[, c('Description', 'GeneRatio', 'BgRatio', 'p.adjust')]))
    }
  }; return(df.go)
}

# KEGG enrichment.
en_kegg <- function(gene, keyType='kegg', organism, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', min.count=5) {
  kk <- enrichKEGG(gene = gene, keyType=keyType, organism = organism, pvalueCutoff = pvalueCutoff, pAdjustMethod=pAdjustMethod, qvalueCutoff=qvalueCutoff)
  df.kk <- as.data.frame(kk)
  if (nrow(df.kk)>0) {
    if (max(df.kk$Count) >= min.count) { 
      print(head(df.kk[, c('Description', 'GeneRatio', 'BgRatio', 'p.adjust')]))
    }
  }; return(df.kk)
}

