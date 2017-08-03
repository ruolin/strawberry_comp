library(polyester)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)

##read expression profile
exp_prof = read.table("data/tx_exp.mat", header=T)
exp_prof.sorted = exp_prof[with(exp_prof, order(transcript_id)),]
gtf.db <- "/home/ruolin/Research/Annotations/hsapiens/gencode.v19.annotation.sqlite"
txdb <- loadDb(gtf.db)
# using the same transcript database (TxDb) as the Geuvadis analysis
# otherwise you can use a different TxDb, with the line:

# make some objects for later
# the exons of each transcript
# a data.frame linking tx to gene
ebt0 <- exonsBy(txdb, by="tx")
gene.ids <- keys(txdb, keytype="GENEID")
gene.ids <- gene.ids[gene.ids != ""]
txdf <- select(txdb, keys=gene.ids, columns=c("TXID","TXCHROM","TXNAME"), keytype="GENEID")
txdf.select <- txdf[txdf$TXNAME %in% exp_prof$transcript_id,]
txdf.select.sorted <- txdf.select[with(txdf.select, order(TXNAME)),]

# how many transcripts per gene?
tab.tx <- table(txdf.select.sorted$GENEID)
head(table(tab.tx))
#can print geneid to txid map
#txdf.select[, c("GENEID","TXNAME")]

sample_depth = c(71347640, 76531554, 56823810 , 70454210, 85193073, 68079019)

sample.tx <- txdf.select.sorted$TXID
sample.tx.names <- txdf.select.sorted$TXNAME
# subset the exons by transcript object
ebt <- ebt0[sample.tx]
head(sample.tx.names)
names(ebt) <- sample.tx.names

##simulate seq
dssl <- getSeq(Hsapiens, ebt)
dna <- unstrsplit(dssl)
library(Rsamtools)
dir.create(file.path("simulated_reads"))
writeXStringSet(dna, file="simulated_reads/transcripts.fa")

txlen = sum(width(ebt))
fpkmmat = as.matrix(exp_prof.sorted[,-c(1,2)])
rownames(fpkmmat) = exp_prof.sorted[,2]
countmat= fpkmmat %*% diag(sample_depth) * txlen / 1e9
# simulation call:
simulate_experiment_countmat(fasta = "simulated_reads/transcripts.fa", readmat = countmat, outdir='simulated_reads', fraglen = 700, readlen = 150) 


fpkm_truth = countmat  %*% diag(1/colSums(countmat)) * (1/txlen) *1e9
fpkm_truth = cbind(rownames(fpkm_truth), fpkm_truth)
colnames(fpkm_truth) = colnames(exp_prof.sorted)[-1]

count_truth = cbind(rownames(countmat), countmat)
colnames(count_truth) = colnames(exp_prof.sorted)[-1]
write.table(fpkm_truth, file="simulated_reads/fpkm_truth.mat", quote = FALSE, row.names = F, sep = "\t")
write.table(count_truth, file="simulated_reads/count_truth.mat", quote = FALSE, row.names = F, sep = "\t")

