base_dir = "/home/ruolin/Research/strawberry_comp/"
library(RColorBrewer)
setwd(base_dir)
mycol = brewer.pal(3,  name="Set1")
##Load all data 
all_strawberry_RD25 = lapply(1:10,
                        function(id){
                          fname = file.path(base_dir, paste0("RD25/strawberry/", id, "/strawberry.stats"))
                          #percent = read.table(fname, skip=10, nrow=6, sep="\t", stringsAsFactors = FALSE)
                          count = read.table(fname, skip=17, nrow=2, sep=":", stringsAsFactors = FALSE)
                          list(count=count)
                        })

all_cufflinks_RD25 = lapply(1:10,
                       function(id){
                         fname = file.path(base_dir, paste0("RD25/cufflinks/", id, "/cufflinks.stats"))
                         #percent = read.table(fname, skip=10, nrow=6, sep="\t", stringsAsFactors = FALSE)
                         count = read.table(fname, skip=17, nrow=2,  sep=":", stringsAsFactors = FALSE)
                         list(count=count)
                       })

all_stringtie_RD25 = lapply(1:10,
                       function(id){
                         fname = file.path(base_dir, paste0("RD25/stringtie/", id, "/stringtie.stats"))
                         #percent = read.table(fname, skip=10, nrow=6, sep="\t", stringsAsFactors = FALSE)
                         count = read.table(fname, skip=17, nrow=2,  sep=":", stringsAsFactors = FALSE)
                         list(count=count)
                       })

all_strawberry_RD60 = lapply(1:10,
                             function(id){
                               fname = file.path(base_dir, paste0("RD60/strawberry/", id, "/strawberry.stats"))
                               #percent = read.table(fname, skip=10, nrow=6, sep="\t", stringsAsFactors = FALSE)
                               count = read.table(fname, skip=17, nrow=2, sep=":", stringsAsFactors = FALSE)
                               list(count=count)
                             })

all_cufflinks_RD60 = lapply(1:10,
                            function(id){
                              fname = file.path(base_dir, paste0("RD60/cufflinks/", id, "/cufflinks.stats"))
                              #percent = read.table(fname, skip=10, nrow=6, sep="\t", stringsAsFactors = FALSE)
                              count = read.table(fname, skip=17, nrow=2,  sep=":", stringsAsFactors = FALSE)
                              list(count=count)
                            })

all_stringtie_RD60 = lapply(1:10,
                            function(id){
                              fname = file.path(base_dir, paste0("RD60/stringtie/", id, "/stringtie.stats"))
                              #percent = read.table(fname, skip=10, nrow=6, sep="\t", stringsAsFactors = FALSE)
                              count = read.table(fname, skip=17, nrow=2,  sep=":", stringsAsFactors = FALSE)
                              list(count=count)
                            })

all_strawberry_RD100 = lapply(1:10,
                             function(id){
                               fname = file.path(base_dir, paste0("RD100/strawberry/", id, "/strawberry.stats"))
                               #percent = read.table(fname, skip=10, nrow=6, sep="\t", stringsAsFactors = FALSE)
                               count = read.table(fname, skip=17, nrow=2, sep=":", stringsAsFactors = FALSE)
                               list(count=count)
                             })

all_cufflinks_RD100 = lapply(1:10,
                            function(id){
                              fname = file.path(base_dir, paste0("RD100/cufflinks/", id, "/cufflinks.stats"))
                              #percent = read.table(fname, skip=10, nrow=6, sep="\t", stringsAsFactors = FALSE)
                              count = read.table(fname, skip=17, nrow=2,  sep=":", stringsAsFactors = FALSE)
                              list(count=count)
                            })

all_stringtie_RD100 = lapply(1:10,
                            function(id){
                              fname = file.path(base_dir, paste0("RD100/stringtie/", id, "/stringtie.stats"))
                              #percent = read.table(fname, skip=10, nrow=6, sep="\t", stringsAsFactors = FALSE)
                              count = read.table(fname, skip=17, nrow=2,  sep=":", stringsAsFactors = FALSE)
                              list(count=count)
                            })


###################
#prepare data
###################

##RD25
parse_data = function( table , software ){
 
  matching_intron_chains = data.frame(
    val = rep(0,10), software=factor(rep(software,10))  
  )
  j = 1;
  for(i in 1:10){
    matching_intron_chains$val[i] = table[[i]]$count[1,2];
  }
  list(matching_intron_chains=matching_intron_chains)
}
strawberry_summary_RD25= parse_data(all_strawberry_RD25, "strawberry")
cufflinks_summary_RD25 = parse_data(all_cufflinks_RD25, "cufflinks")
stringtie_summary_RD25 = parse_data(all_stringtie_RD25, "stringtie")



if(exists("RD25"))rm(RD25)
RD25 <- rbind(unlist(stringtie_summary_RD25$matching_intron_chain$val),
              unlist(strawberry_summary_RD25$matching_intron_chains$val), 
              unlist(cufflinks_summary_RD25$matching_intron_chains$val))
rownames(RD25) = c("stringtie", "strawberry", "cufflinks")
mean=apply(RD25, 1, mean)
se = apply(RD25, 1, function(x)
                         {sd(x)/sqrt(length(x))}
          )

RD25=cbind(RD25, mean, se, rd=rep(25,3))
rm(mean)
rm(se)


### RD60
parse_data = function( table , software ){
  
  matching_intron_chains = data.frame(
    val = rep(0,10), software=factor(rep(software,10))  
  )
  j = 1;
  for(i in 1:10){
    matching_intron_chains$val[i] = table[[i]]$count[1,2];
  }
  list(matching_intron_chains=matching_intron_chains)
}
strawberry_summary_RD60= parse_data(all_strawberry_RD60, "strawberry")
cufflinks_summary_RD60 = parse_data(all_cufflinks_RD60, "cufflinks")
stringtie_summary_RD60 = parse_data(all_stringtie_RD60, "stringtie")



if(exists("RD60"))rm(RD60)
RD60 <- rbind(unlist(stringtie_summary_RD60$matching_intron_chain$val),
              unlist(strawberry_summary_RD60$matching_intron_chains$val), 
              unlist(cufflinks_summary_RD60$matching_intron_chains$val))
rownames(RD60) = c("stringtie", "strawberry", "cufflinks")
mean=apply(RD60, 1, mean)
se = apply(RD60, 1, function(x)
{sd(x)/sqrt(length(x))}
)


RD60=cbind(RD60, mean, se, rd=rep(60,3))
rm(mean)
rm(se)


#RD100
parse_data = function( table , software ){
  
  matching_intron_chains = data.frame(
    val = rep(0,10), software=factor(rep(software,10))  
  )
  j = 1;
  for(i in 1:10){
    matching_intron_chains$val[i] = table[[i]]$count[1,2];
  }
  list(matching_intron_chains=matching_intron_chains)
}
strawberry_summary_RD100= parse_data(all_strawberry_RD100, "strawberry")
cufflinks_summary_RD100 = parse_data(all_cufflinks_RD100, "cufflinks")
stringtie_summary_RD100 = parse_data(all_stringtie_RD100, "stringtie")



if(exists("RD100"))rm(RD100)
RD100 <- rbind(unlist(stringtie_summary_RD100$matching_intron_chain$val),
              unlist(strawberry_summary_RD100$matching_intron_chains$val), 
              unlist(cufflinks_summary_RD100$matching_intron_chains$val))
rownames(RD100) = c("stringtie", "strawberry", "cufflinks")
mean=apply(RD100, 1, mean)
se = apply(RD100, 1, function(x)
{sd(x)/sqrt(length(x))}
)

RD100=cbind(RD100, mean, se, rd=rep(100,3))
rm(mean)
rm(se)

####### combine

all_data = rbind(RD25, RD60, RD100)

#### Plot

pdf("barplot.pdf")

barCenters = barplot(all_data[,"mean"], col=mycol, 
                     space=c(0,0,0,0.5,0,0,0.5,0,0),
                     ylim=c(0,7000), xaxt='n', yaxt='n')
axis(1, at=barCenters[c(2,5,8)] , label=c("RD25", "RD60", "RD100"), 
     tck=0, mgp=c(0,0.5,0), lwd=0, cex.axis=1.3)
axis(2, at=c(0,1000,2000,3000,4000,5000,6000,7000), las=2, lwd=2, cex.axis=1.3)

## Error bars
segments(barCenters, all_data[,"mean"]-all_data[,"se"]*2, barCenters,
          all_data[,"mean"] + all_data[,"se"]*2, lwd=1.5)
arrows(barCenters, all_data[,"mean"] - all_data[,"se"] * 2, barCenters,
       all_data[,"mean"] + all_data[,"se"] * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)


legend('topleft',c('StringTie','Strawberry', 'Cufflinks'),
       fill = mycol, bty = 'n', cex=1.4,
       border = NA)

(3914-2875 + 5666-4451 + 6533-5357)/3

(6533-5789 + 5666-5068 + 3914-3617)/3

1250+1350+1212

dev.off()



