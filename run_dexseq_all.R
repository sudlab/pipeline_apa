library(rtracklayer)
library(DEXSeq)
library(experimentr)
library(GenomicRanges)
library(optparse)

opts = list(make_option(c("-p","--processes"),
                        type="integer",
                        dest="processes",
                        default=1,
                        help="Number of procsses to use. Using multiple processes requires the BiocParallel library"))

args = Experiment.start(opts)
gffs <- import(args$infiles[1])

col_data = read.delim(args$infile[3], header=T, row.names = 1)
col_data$condition <- relevel(col_data$condition, "Control")
rownames(col_data) <- gsub("-","_",rownames(col_data))

counts <- read.delim(args$infiles[2])
names(counts) <- gsub(".","_",names(counts), fixed=TRUE)
counts <- counts[, rownames(col_data)]
counts <- as.matrix(counts)

cat("# design is:\n")
print(col_data)
cat('# got ', dim(counts)[2], ' samples\n')
cat('# building data set\n')

dxd <- DEXSeqDataSet(counts, col_data, design= ~ sample + exon + condition:exon,
                           featureID = mcols(gffs)$exon_id, 
                           groupID = mcols(gffs)$gene_id,
                           featureRanges=gffs)

if (args$processes > 1) {
  BPPARAM = MulticoreParam(workers=args$processes)
  dxd_results <- DEXSeq(dxd, BPPARAM=BPPARAM, quiet=F)
} else {
  dxd_results <- DEXSeq(dxd, quiet=F)
}
cat('# computing results\n')

cat('#Subset results to only introns')
dxd_results$padj <- p.adjust(dxd_results$pvalue, method="BH")
sig_exons <- rownames(subset(dxd_results, padj < 0.05))
cat('# got', length(sig_exons), ' significant exons\n')
sig_granges <- rowRanges(dxd[sig_exons,])
dxd_df <-as.data.frame(dxd_results)
cat('# outputting results')
write.table(dxd_df[,1:15],args$outfiles[1],
            quote=F, sep = "\t", row.names=FALSE)
export(sig_granges, args$outfiles[2])
save.image(args$outfiles[3])
Experiment.stop()
