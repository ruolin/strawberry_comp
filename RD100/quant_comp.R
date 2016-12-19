#Running Time
# StringTie:   0m38.840s
# Strawberry:  1m38.146s
# Cufflinks:  10m19.780s 
library(RColorBrewer)


base_dir = "/home/ruolin/Research/strawberry_comp/RD100"

setwd(base_dir)

ARD <-function(x,y){
  if(x==y) {
    return (0)
  }else{
    num = abs(x-y)
    den = 0.5*abs(x+y)
    return (num/den)
  }
}


prop_corr <- function(x,y){
   num = 2* cov( log(x+1), log(y+1))
   den = var(log(x+1)) + var(log(y+1))
   return (num/den)
}

MARD = function(x,y){
  mean(mapply(ARD, x, y))
}

##Read ground truth
all_oracle <- lapply(1:10,
                     function(id) {
                       fname <- file.path(base_dir, paste0(id,"/RD100.control_", id, ".pro"))
                       result <- read.table(fname, header = FALSE, stringsAsFactors = FALSE) 
                       result
                     })

all_RD25_truth  = lapply(1:10,
          function(id){
            all_oracle[[id]]$true_FPKM = all_oracle[[id]]$V10*(10^3/all_oracle[[id]]$V4) * 
                                  ((10^6)/sum(all_oracle[[id]]$V10))
            all_oracle[[id]]
            
          })
rm(all_oracle)

## Read Strawberry
all_strawberry = lapply(1:10,
           function(id){
             fname = file.path(base_dir, paste0("strawberry/", id, "/strawberry.assembled_transcripts.gtf.tmap"))
             result = read.table(fname, header=TRUE, stringsAsFactors = FALSE)
             result
           })
## read StringTie
all_stringtie = lapply(1:10,
                       function(id){
                         fname = file.path(base_dir, paste0("stringtie/", id, "/stringtie.transcripts.gtf.tmap" ))
                         result = read.table(fname, header=TRUE, stringsAsFactors = FALSE)
                       })

## read Cufflinks
all_cufflinks = lapply(1:10,
                       function(id){
                         fname = file.path(base_dir, paste0("cufflinks/", id, "/cufflinks.transcripts.gtf.tmap"))
                         result = read.table(fname, header=TRUE, stringsAsFactors = FALSE)
                       })

## Correlation matrix
all_strawberry_corr = lapply(1:10,
                  function(id){
                    result = merge(all_strawberry[[id]], all_RD25_truth[[id]], by.x="ref_id", by.y="V2")
                    result
                  })

all_stringtie_corr = lapply(1:10,
                       function(id){
                         result = merge(all_stringtie[[id]], all_RD25_truth[[id]], by.x="ref_id", by.y="V2")
                       })
all_cufflinks_corr = lapply(1:10,
                            function(id){
                              result = merge(all_cufflinks[[id]], all_RD25_truth[[id]], by.x="ref_id", by.y="V2")
                            })

### proportional correlation
strawberry_prop_corr = lapply(1:10,
                       function(id){
                         result = prop_corr(all_strawberry_corr[[id]]$FPKM, all_strawberry_corr[[id]]$true_FPKM) 
                         result
                       })

stringtie_prop_corr = lapply(1:10,
                       function(id){
                         result = prop_corr(all_stringtie_corr[[id]]$FPKM, all_stringtie_corr[[id]]$true_FPKM) 
                         result
                       })

cufflinks_prop_corr = lapply(1:10,
                             function(id){
                               result = prop_corr(all_cufflinks_corr[[id]]$FPKM, all_cufflinks_corr[[id]]$true_FPKM)
                               result
                             })

### spearman correlation
strawberry_sp_corr = lapply(1:10,
                            function(id){
                              sp = cor(all_strawberry_corr[[id]]$FPKM, all_strawberry_corr[[id]]$true_FPKM, method="spearman")
                              sp          
                              })

stringtie_sp_corr = lapply(1:10,
                            function(id){
                              sp = cor(all_stringtie_corr[[id]]$FPKM, all_stringtie_corr[[id]]$true_FPKM, method="spearman")
                              sp          
                            })

cufflinks_sp_corr = lapply(1:10,
                            function(id){
                              sp = cor(all_cufflinks_corr[[id]]$FPKM, all_cufflinks_corr[[id]]$true_FPKM, method="spearman")
                              sp          
                            })


### MARD
strawberry_MARD = lapply(1:10,
                            function(id){
                              mard = MARD(all_strawberry_corr[[id]]$FPKM, all_strawberry_corr[[id]]$true_FPKM)
                              mard          
                            })

stringtie_MARD = lapply(1:10,
                         function(id){
                           mard = MARD(all_stringtie_corr[[id]]$FPKM, all_stringtie_corr[[id]]$true_FPKM)
                           mard          
                         })

cufflinks_MARD = lapply(1:10,
                         function(id){
                           mard = MARD(all_cufflinks_corr[[id]]$FPKM, all_cufflinks_corr[[id]]$true_FPKM)
                           mard          
                         })
##plot
#mycol = brewer.pal(3,  name="Pastel1")
mycol = brewer.pal(3,  name="Set1")
pdf("RD100-quant.pdf")

op = par(mfrow= c(3,1),
  mar = c(4,3,0,0.5) + 0.1,
  oma = c(0,1,1,0.5) + 0.1
  )


hist(unlist(stringtie_prop_corr), col=mycol[1], 
     xlim=c(0.8,0.95), ylim=c(0,5), yaxs = "i", xaxs = "i", 
     xaxt="n", yaxt="n", xlab="", ylab="", main="")
title( ylab="Frequency", xlab="Proportional correlation", main="", mgp=c(2,0,0), cex.lab=1.2)
axis(1, lwd=2, tck=0, mgp=c(1,0.5,0), cex.axis=1.2)
axis(2, lwd=2, tck=0, las=2, mgp=c(1,0.5,0), cex.axis=1.2)
hist(unlist(strawberry_prop_corr), col=mycol[2], add=T, breaks=6)
hist(unlist(cufflinks_prop_corr), col=mycol[3], add=T, breaks=30)
legend('topleft',c('StringTie','Strawberry', 'Cufflinks'),
       fill = mycol, bty = 'n', cex=1.4,
       border = NA)

mean(unlist(strawberry_prop_corr)) ## mean strawberry
mean(unlist(stringtie_prop_corr)) ## mean stringtie
mean(unlist(cufflinks_prop_corr)) ## mean cufflinks

hist(unlist(stringtie_sp_corr), col=mycol[1], 
     xlim=c(0.8,0.95), ylim=c(0,5), yaxs = "i", xaxs = "i", 
     xaxt="n", yaxt="n", xlab="", ylab="", main="", breaks=12)
title( ylab="Frequency", xlab="Spearman correlation", main="", mgp=c(2,0,0), cex.lab=1.2)
axis(1, lwd=2, tck=0, mgp=c(1,0.5,0), cex.axis=1.2)
axis(2, lwd=2, tck=0, las=2, mgp=c(1,0.5,0), cex.axis=1.2)
hist(unlist(strawberry_sp_corr), col=mycol[2], add=T, breaks=3)
hist(unlist(cufflinks_sp_corr), col=mycol[3], add=T, breaks=10)

mean(unlist(strawberry_sp_corr)) ## mean strawberry
mean(unlist(stringtie_sp_corr)) ## mean stringtie
mean(unlist(cufflinks_sp_corr)) ## mean cufflinks

hist(unlist(stringtie_MARD), col=mycol[1], 
     xlim=c(0.36, 0.49), ylim=c(0,5), yaxs = "i", xaxs = "i", 
     xaxt="n", yaxt="n", xlab="", ylab="", main="", breaks=12)
title( ylab="Frequency", xlab="MARD", main="", mgp=c(2,0,0), cex.lab=1.2)
axis(1, lwd=2, tck=0, mgp=c(1,0.5,0), cex.axis=1.2, at=c(0.36,0.38,0.4,0.42,0.44,0.46,0.48))
axis(2, lwd=2, tck=0, las=2, mgp=c(1,0.5,0), cex.axis=1.2)
hist(unlist(strawberry_MARD), col=mycol[2], add=T, breaks=10)
hist(unlist(cufflinks_MARD), col=mycol[3], add=T, breaks=10)

mean(unlist(strawberry_MARD)) ## mean strawberry
mean(unlist(stringtie_MARD)) ## mean stringtie
mean(unlist(cufflinks_MARD)) ## mean cufflinks

par(op)
dev.off()

