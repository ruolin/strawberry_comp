base_dir = "/home/ruolin/git/strawberry_comp/geuvadis"

setwd(base_dir)

prop_corr <- function(x,y){
  num = 2* cov( log(x+1), log(y+1))
  den = var(log(x+1)) + var(log(y+1))
  return (num/den)
}

##Read truth
#truth = "simulated_reads/fpkm_truth.mat"
#truth = "simulated_reads/count_truth.mat"
all_truth = read.table(truth, header=T)
#all_truth[all_truth$transcript_id == 'ENST00000437048.2',]
## Read Strawberry
all_strawberry = lapply(1:1,
                        function(id){
                          fname = file.path(base_dir, paste0("strawberry/", id, "/strawberry.assembled_transcripts.gtf.tmap"))
                          result = read.table(fname, header=TRUE, stringsAsFactors = FALSE)
                          result
                        })
## read StringTie
all_stringtie = lapply(1:1,
                       function(id){
                         fname = file.path(base_dir, paste0("stringtie/", id, "/stringtie.transcripts.gtf.tmap" ))
                         result = read.table(fname, header=TRUE, stringsAsFactors = FALSE)
                       })

## read Cufflinks
all_cufflinks = lapply(1:1,
                       function(id){
                         fname = file.path(base_dir, paste0("cufflinks/", id, "/cufflinks.transcripts.gtf.tmap"))
                         result = read.table(fname, header=TRUE, stringsAsFactors = FALSE)
                       })



merged_string = merge (all_stringtie[[1]], all_truth[,c(1:2)], by.x="ref_id", by.y="transcript_id")
cor(merged_string$FPKM, merged_string$ERR188021.gtf, method="spearman")
prop_corr(merged_string$FPKM, merged_string$ERR188021.gtf)

merged_straw = merge (all_strawberry[[1]], all_truth[,c(1:2)], by.x="ref_id", by.y="transcript_id")
cor(merged_straw$FPKM, merged_straw$ERR188021.gtf, method="spearman")
prop_corr(merged_straw$FPKM, merged_straw$ERR188021.gtf)

merged_cuff = merge (all_cufflinks[[1]], all_truth[,c(1:2)], by.x="ref_id", by.y="transcript_id")
cor(merged_cuff$FPKM, merged_cuff$ERR188021.gtf, method="spearman")
prop_corr(merged_cuff$FPKM, merged_cuff$ERR188021.gtf)


straw_string_merge = merge(all_strawberry[[1]], all_strawberry[[1]], by="ref_id")
cor(straw_string_merge$FPKM.x, straw_string_merge$FPKM.y)
straw_cuff_merge = merge(all_strawberry[[1]], all_cufflinks[[1]], by="ref_id")
cor(straw_cuff_merge$FPKM.x, straw_cuff_merge$FPKM.y)
string_cuff_merge = merge(all_stringtie[[1]], all_cufflinks[[1]], by="ref_id")
cor(string_cuff_merge$FPKM.x, string_cuff_merge$FPKM.y)


######
###DEBUG
cuff_tp = subset(all_cufflinks[[1]], class_code == "=")$ref_id
straw_tp = subset(all_strawberry[[1]], class_code == "=")$ref_id
string_tp = subset(all_stringtie[[1]], class_code == "=")$ref_id
setdiff(string_tp, straw_tp)
setdiff(straw_tp, string_tp)
setdiff(cuff_tp, straw_tp)
setdiff(string_tp, cuff_tp)
subset(all_strawberry[[1]], ref_id %in% setdiff(string_tp, straw_tp))

subset(all_cufflinks[[1]], ref_id %in% setdiff(straw_tp,cuff_tp))
subset(all_strawberry[[1]], ref_id %in% setdiff(cuff_tp, straw_tp))
######

x = rgamma(10000, shape=250, rate=500)
plot(density(x))
