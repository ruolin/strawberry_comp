base_dir = "/home/ruolin/git/strawberry_comp/geuvadis"
setwd(base_dir)
library(RColorBrewer)
mycol = brewer.pal(3,  name="Set1")

straw_asmb = read.table(header=F, file="assemb_straw.stats",row.names = 1)
cuff_asmb = read.table(header=F, file="assemb_cuff.stats",row.names = 1)
string_asmb = read.table(header=F, file="assemb_string.stats",row.names = 1)

num_software = 3
num_row = 6
num_sample =6
for (i in 5:6) {
  i = 6
  collect = matrix(c(rep(2,num_sample), as.numeric(straw_asmb[i,])), byrow=F, ncol=2)
  collect = rbind(collect, matrix(c(rep(3,num_sample), as.numeric(cuff_asmb[i,])), byrow=F, ncol=2) )
  collect = rbind(collect, matrix(c(rep(1,num_sample), as.numeric(string_asmb[i,])), byrow=F, ncol=2) )
  colnames(collect) = c("software", "score")
  boxplot(score ~ software, data= collect, col = mycol )
}

prop_corr <- function(x,y){
  num = 2* cov( log(x+1), log(y+1))
  den = var(log(x+1)) + var(log(y+1))
  return (num/den)
}

ARD <-function(x,y){
  if(x==y) {
    return (0)
  }else{
    num = abs(x-y)
    den = 0.5*abs(x+y)
    return (num/den)
  }
}

MARD = function(x,y){
  mean(mapply(ARD, x, y))
}


##Read truth
truth = "simulated_reads/fpkm_truth.mat"
#truth = "simulated_reads/count_truth.mat"
all_truth = read.table(truth, header=T)
#all_truth[all_truth$transcript_id == 'ENST00000437048.2',]
## Read Strawberry
all_strawberry = lapply(1:6,
                        function(id){
                          fname = file.path(base_dir, paste0("strawberry/", id, "/strawberry.assembled_transcripts.gtf.tmap"))
                          result = read.table(fname, header=TRUE, stringsAsFactors = FALSE)
                          result
                        })
## read StringTie
all_stringtie = lapply(1:6,
                       function(id){
                         fname = file.path(base_dir, paste0("stringtie/", id, "/stringtie.transcripts.gtf.tmap" ))
                         result = read.table(fname, header=TRUE, stringsAsFactors = FALSE)
                       })

## read Cufflinks
all_cufflinks = lapply(1:6,
                       function(id){
                         fname = file.path(base_dir, paste0("cufflinks/", id, "/cufflinks.transcripts.gtf.tmap"))
                         result = read.table(fname, header=TRUE, stringsAsFactors = FALSE)
                       })

all_corr <- function(all_pred, all_truth, method) {
                stopifnot(method %in% c("spearman", "proportion", "pearson"))
                result = lapply(1:6,
                               function(id){
                                 merged = merge(all_pred[[id]], all_truth[,c(1,id+1)], by.x="ref_id", by.y="transcript_id")
                                 if (method == "spearman") {
                                    cor(merged$FPKM, merged[,ncol(merged)], method="spearman")
                                 }
                                 else if (method == "pearson") {
                                   cor(merged$FPKM, merged[,ncol(merged)], method="pearson")
                                 }
                                 else {
                                   prop_corr(merged$FPKM, merged[,ncol(merged)])
                                 }
                               })
                result
}



all_MARD <- function(all_pred, all_truth) {
  result = lapply(1:6,
                  function(id){
                    merged = merge(all_pred[[id]], all_truth[,c(1,id+1)], by.x="ref_id", by.y="transcript_id")
                    MARD(merged$FPKM, merged[,ncol(merged)])
                  })
  result
}

strawberry_MARD = all_MARD(all_strawberry,all_truth)
stringtie_MARD = all_MARD(all_stringtie, all_truth)
cufflinks_MARD = all_MARD(all_cufflinks, all_truth)

strawberry_sp_corr = all_corr(all_strawberry,all_truth, method = "spearman")
stringtie_sp_corr = all_corr(all_stringtie, all_truth, method = "spearman")
cufflinks_sp_corr = all_corr(all_cufflinks, all_truth, method = "spearman")

strawberry_prop_corr = all_corr(all_strawberry,all_truth, method = "proportion")
stringtie_prop_corr = all_corr(all_stringtie, all_truth, method = "proportion")
cufflinks_prop_corr = all_corr(all_cufflinks, all_truth, method = "proportion")

pdf("geuvadis.pdf")

op = par(mfrow= c(3,1),
         mar = c(4,3,0,0.5) + 0.1,
         oma = c(0,1,1,0.5) + 0.1
)


hist(unlist(stringtie_prop_corr), col=mycol[1], 
     xlim=c(0.90,1.0), ylim=c(0,5), yaxs = "i", xaxs = "i", 
     xaxt="n", yaxt="n", xlab="", ylab="", main="", breaks=50)
title( ylab="Frequency", xlab="Proportional correlation", main="", mgp=c(2,0,0), cex.lab=1.2)
axis(1, lwd=2, tck=0, mgp=c(1,0.5,0), cex.axis=1.2)
axis(2, lwd=2, tck=0, las=2, mgp=c(1,0.5,0), cex.axis=1.2)
hist(unlist(strawberry_prop_corr), col=mycol[2], add=T, breaks=30)
hist(unlist(cufflinks_prop_corr), col=mycol[3], add=T, breaks=15)
legend('topleft',c('StringTie','Strawberry', 'Cufflinks'),
       fill = mycol, bty = 'n', cex=1.4,
       border = NA)

hist(unlist(stringtie_sp_corr), col=mycol[1], 
     xlim=c(0.9,1), ylim=c(0,5), yaxs = "i", xaxs = "i", 
     xaxt="n", yaxt="n", xlab="", ylab="", main="", breaks=50)
title( ylab="Frequency", xlab="Spearman correlation", main="", mgp=c(2,0,0), cex.lab=1.2)
axis(1, lwd=2, tck=0, mgp=c(1,0.5,0), cex.axis=1.2)
axis(2, lwd=2, tck=0, las=2, mgp=c(1,0.5,0), cex.axis=1.2)
hist(unlist(strawberry_sp_corr), col=mycol[2], add=T, breaks=30)
hist(unlist(cufflinks_sp_corr), col=mycol[3], add=T, breaks=5)


hist(unlist(stringtie_MARD), col=mycol[1], 
     xlim=c(0.1, 0.4), ylim=c(0,5), yaxs = "i", xaxs = "i", 
     xaxt="n", yaxt="n", xlab="", ylab="", main="", breaks=20)
title( ylab="Frequency", xlab="MARD", main="", mgp=c(2,0,0), cex.lab=1.2)
axis(1, lwd=2, tck=0, mgp=c(1,0.5,0), cex.axis=1.2, at=c(0,1, 0.2,0.3,0.4))
axis(2, lwd=2, tck=0, las=2, mgp=c(1,0.5,0), cex.axis=1.2)
hist(unlist(strawberry_MARD), col=mycol[2], add=T, breaks=10)
hist(unlist(cufflinks_MARD), col=mycol[3], add=T, breaks=20)

par(op)
dev.off()


##############################
##############################
##########ASSEMBLY############

##Load all data 
all_strawberry = lapply(1:6,
                        function(id){
                          fname = file.path(base_dir, paste0("strawberry/", id, "/strawberry.stats"))
                          percent = read.table(fname, skip=10, nrow=6, sep="\t", stringsAsFactors = FALSE)
                          count = read.table(fname, skip=17, nrow=2, sep="\t", stringsAsFactors = FALSE)
                          list(percent=percent, count=count)
                        })



######
###DEBUG
cuff_fp = subset(all_cufflinks[[2]], class_code == "=")$ref_id
straw_tp = subset(all_strawberry[[2]], class_code == "=")$ref_id
string_tp = subset(all_stringtie[[2]], class_code == "=")$ref_id
setdiff(string_tp, straw_tp)
setdiff(straw_tp, string_tp)

setdiff(cuff_tp, straw_tp)
setdiff(string_tp, cuff_tp)
subset(all_strawberry[[1]], ref_id %in% setdiff(string_tp, straw_tp))

subset(all_cufflinks[[1]], ref_id %in% setdiff(straw_tp,cuff_tp))
subset(all_strawberry[[1]], ref_id %in% setdiff(cuff_tp, straw_tp))
######
all_strawberry[[2]]
x = rgamma(10000, shape=250, rate=500)
plot(density(x))
