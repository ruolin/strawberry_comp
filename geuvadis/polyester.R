library(polyester)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)

gtf.file <- "/home/ruolin/Research/Annotations/hsapiens/gencode.v19.annotation.sqlite"
txdb <- loadDb(gtf.file)
# using the same transcript database (TxDb) as the Geuvadis analysis
# otherwise you can use a different TxDb, with the line:
#gtf.file <- "/home/ruolin/Research/Annotations/hsapiens/gencode.v19.annotation.gff3"
#txdb <- makeTxDbFromGFF(gtf.file)

# make some objects for later
# the exons of each transcript
# a data.frame linking tx to gene
ebt0 <- exonsBy(txdb, by="tx")
gene.ids <- keys(txdb, keytype="GENEID")
gene.ids <- gene.ids[gene.ids != ""]
txdf <- select(txdb, keys=gene.ids, columns=c("TXID","TXCHROM","TXNAME"), keytype="GENEID")
txdf <- txdf[txdf$TXNAME %in%,]



# FASTA annotation
fasta_file = system.file('extdata', 'chr22.fa', package='polyester')
fasta = readDNAStringSet(fasta_file)

num_timepoints = 12
# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
readspertx = round(20 * width(fasta) / 100)
length(fasta)
countmat = matrix(readspertx, nrow=length(fasta), ncol=num_timepoints)
countmat
# add spikes in expression at certain timepoints to certain transcripts:
up_early = c(1,2) 
up_late = c(3,4)
countmat[up_early, 2] = 3*countmat[up_early, 2]
countmat[up_early, 3] = round(1.5*countmat[up_early, 3])
countmat[up_late, 10] = 6*countmat[up_late, 10]
countmat[up_late, 11] = round(1.2*countmat[up_late, 11])

# simulation call:
simulate_experiment_countmat(fasta_file, readmat=countmat, outdir='simulated_reads') 

