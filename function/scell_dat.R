# Read single cell data.
scell_dat <- function(rds, df.match, dev.zone=FALSE) {
  library(Matrix)
  rds <- readRDS(rds)
  # Single-cell metadata from processed data.
  met.sc <- rds@meta.data
  spl.sct <- rds@assays$spliced_SCT
  spl.sct.cnt <- spl.sct@counts
  # Abbreviate cell anno.
  met.sc$cell.anno <- make.names(met.sc$celltype.anno)
  cell.na <- c('lat.rt.cap', 'procam', 'endo', 'tricho', 'xylem.po.per', 'phlo.po.per', 'atricho', 'colu', 'protoxylem', 'cortex', 'metaphlo.comp.cell', 'metaxylem', 'quies.cent', 'protophlo', 'per', 'phlo', 'xylem', 'put.quies.cent', 'stem.niche')

  names(cell.na) <- c('Lateral.Root.Cap', 'Procambium', 'Endodermis', 'Trichoblast', 'Xylem.Pole.Pericycle', 'Phloem.Pole.Pericycle', 'Atrichoblast', 'Columella', 'Protoxylem', 'Cortex', 'Metaphloem...Companion.Cell', 'Metaxylem', 'Quiescent.Center', 'Protophloem', 'Pericycle', 'Phloem', 'Xylem', 'Putative.Quiescent.Center', 'Stem.Cell.Niche')
  cell.ann.uni <- unique(met.sc$cell.anno)
  id0 <- cell.ann.uni[which(!cell.ann.uni %in% names(cell.na))]
  if (length(id0)>0) { print(id0); stop('New cell names are detected!') }
  for (i in seq_along(cell.na)) {
    cell.anno <- met.sc$cell.anno
    met.sc$cell.anno[cell.anno==names(cell.na[i])] <- cell.na[i]
  }
  colnames(spl.sct.cnt) <- met.sc$cell.anno

  # Abbreviate zones.
  met.sc$zone <- make.names(met.sc$time.anno)
  zone <- c('matur', 'elong', 'meri', 'proxi.lrc', 'proxi.colu', 'dist.lrc', 'dist.colu')
  names(zone) <- make.names(c('Maturation', 'Elongation', 'Meristem', 'Proximal Lateral Root Cap', 'Proximal Columella', 'Distal Lateral Root Cap', 'Distal Columella'))
  for (i in seq_along(zone)) {
    zone.all <- met.sc$zone
    met.sc$zone[zone.all==names(zone[i])] <- zone[i]
  }
  met.sc <- met.sc[, c('orig.ident', 'celltype.anno', 'cell.anno', 'time.anno', 'zone')]
  met.sc$cell.zone <- paste0(make.names(met.sc$cell.anno), '.', met.sc$zone)

  # Transfer cell ids from processed to raw data.
  # Leave out zones.
  # cell.ids <- make.names(met.sc$cell.anno)
  # Assign zones to cells.
  cell.ids <- make.names(met.sc$cell.zone)
  names(cell.ids) <- rownames(met.sc)

  sc.all <- rds@assays$spliced_RNA@counts+rds@assays$unspliced_RNA@counts
  colnames(sc.all) <- sub('_1.*', '', colnames(sc.all))
  colnames(sc.all) <- cell.ids[colnames(sc.all)]
  sc.all <- sc.all[, !is.na(colnames(sc.all))]
  if (dev.zone==FALSE) colnames(sc.all) <- sub('\\.elong$|\\.matur$|\\.meri$', '', colnames(sc.all)) 
  # dgCMatrix: columns should not be re-ordered by sorted column names. 
  # sc.all <- sc.all[, sort(colnames(sc.all))]
  na.uni <- unique(colnames(sc.all))
  no <- na.uni[!na.uni %in% df.match$cell]
  if (length(no) > 0) stop(paste0("These cells are not in the 'df.match': ", paste0(no, collapse=', '), "!"))
  return(sc.all)
}
# Mouse kidney single cell data.
sc_dat_mus_kdn <- function(sc.pa=NULL) {
  library(data.table)
  sc.kdn <- fread(sc.pa, fill=TRUE)
  clus <- unlist(sc.kdn[1, ])[-1]
  rna <- sc.kdn$'AAACCTGAGATATGCA-1'[-1]
  sc.kdn <- as(sc.kdn[-1, -1], 'matrix')
  cell <- c('1'='endo', '2'='podo', '3'='proxi.tub', '4'='loop.hen', '5'='distal.con.tub', '6'='col.duct.prin.cell', '7'='col.duct.inter.cell', '8'='col.duct.trans.cell', '9'='novel.cell1', '10'='fibrob', '11'='macro', '12'='neutro', '13'='b.lymph', '14'='t.lymph', '15'='natural.killer', '16'='novel.cell2') 
  sum(duplicated(rna))
  rownames(sc.kdn) <- rna; colnames(sc.kdn) <- cell[clus]
  return(sc.kdn)
}
# Mouse brain.
sc_dat_mus_brain <- function(sc.pa=NULL, meta.pa=NULL) {
  library(data.table)
  sc.brain <- fread(sc.pa); sc.brain.met <- fread(meta.pa)
  # Check matching between counts and metadata.
  setkey(sc.brain, V1); setkey(sc.brain.met, V1)
  all(sc.brain$V1==sc.brain.met$V1)
  sc.brain <- sc.brain[sc.brain.met$passed_QC==TRUE]
  sc.brain.met <- sc.brain.met[passed_QC==TRUE]
  all(sc.brain$V1==sc.brain.met$V1)
  sc.brain$V1 <- make.names(sc.brain$V1)
  sc.brain.met$V1 <- make.names(sc.brain.met$V1)

  # Transpose count table.
  rna <- sc.brain$V1; sc.brain <- as(sc.brain[, -1], 'matrix')
  rownames(sc.brain) <- rna; sc.brain <- t(sc.brain)

  # Abbreviate cell names in metadata.
  sc.brain.met$ABA_parent <- make.names(sc.brain.met$ABA_parent)
  cells <- c(Cerebellum='cere', fiber.tracts='fiber.tr', Hippocampal.region='hipp', Isocortex='isocort', Olfactory.areas='olfa', Retrohippocampal.region='retrohipp', Thalamus='thalamus', ventricular.systems='ventri', Cortical.subplate='corti.sub', Hindbrain='hindbrain', Hypothalamus='hypotha', Midbrain='midbrain', Pallidum='pallidum', Striatum='striatum', Undefined.areas='undefined')

  sc.brain.met$cell <- cells[sc.brain.met$ABA_parent]
  # Assign cells names in metadata to count table.
  all(colnames(sc.brain)==sc.brain.met$V1)
  colnames(sc.brain) <- sc.brain.met$cell; return(sc.brain)
}
