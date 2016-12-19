library(RColorBrewer)
mycol = brewer.pal(3,  name="Set1")
base_dir = "/home/ruolin/Research/strawberry_comp/RD100"

setwd(base_dir)

##Load all data 
all_strawberry = lapply(1:10,
                        function(id){
                          fname = file.path(base_dir, paste0("strawberry/", id, "/strawberry.stats"))
                          percent = read.table(fname, skip=10, nrow=6, sep="\t", stringsAsFactors = FALSE)
                          count = read.table(fname, skip=17, nrow=2, sep="\t", stringsAsFactors = FALSE)
                          list(percent=percent, count=count)
                        })

all_cufflinks = lapply(1:10,
                        function(id){
                          fname = file.path(base_dir, paste0("cufflinks/", id, "/cufflinks.stats"))
                          percent = read.table(fname, skip=10, nrow=6, sep="\t", stringsAsFactors = FALSE)
                          count = read.table(fname, skip=17, nrow=2, sep="\t", stringsAsFactors = FALSE)
                          list(percent=percent, count=count)
                        })

all_stringtie = lapply(1:10,
                       function(id){
                         fname = file.path(base_dir, paste0("stringtie/", id, "/stringtie.stats"))
                         percent = read.table(fname, skip=10, nrow=6, sep="\t", stringsAsFactors = FALSE)
                         count = read.table(fname, skip=17, nrow=2, sep="\t", stringsAsFactors = FALSE)
                         list(percent=percent, count=count)
                       })

### parse the data table
parse_data = function( table , software ){
  base_level = data.frame(
    label = rep(c(1,2),10), val = rep(0,20), software=factor(rep(software,20))
  )
  exon_level = data.frame(
    label = rep(c(1,2),10), val = rep(0,20), software=factor(rep(software,20))
  )
  intron_level = data.frame(
    label = rep(c(1,2),10), val = rep(0,20), software=factor(rep(software,20))
  )
  transcript_level = data.frame(
    label = rep(c(1,2),10), val = rep(0,20), software=factor(rep(software,20))
  )
  j = 1;
  for(i in 1:10){
    base_level$val[j] = table[[i]]$percent[1,2];
    exon_level$val[j] = table[[i]]$percent[2,2];
    intron_level$val[j] = table[[i]]$percent[3,2];
    transcript_level$val[j] = table[[i]]$percent[4,2];
    j = j+1;

    base_level$val[j] = table[[i]]$percent[1,3];
    exon_level$val[j] = table[[i]]$percent[2,3];
    intron_level$val[j] = table[[i]]$percent[3,3];
    transcript_level$val[j] = table[[i]]$percent[4,3];
    j = j+1;
    
  }
  list(base_level=base_level, exon_level= exon_level, 
       intron_level = intron_level, transcript_level=transcript_level)
}
 strawberry_summary= parse_data(all_strawberry, "strawberry")
 cufflinks_summary = parse_data(all_cufflinks, "cufflinks")
 stringtie_summary = parse_data(all_stringtie, "stringtie")

#generating summaries for all the plots
all_strawberry[[2]]$percent
base_summary = rbind(strawberry_summary$base_level, cufflinks_summary$base_level, stringtie_summary$base_level)
exon_summary = rbind(strawberry_summary$exon_level, cufflinks_summary$exon_level, stringtie_summary$exon_level)
intron_summary = rbind(strawberry_summary$intron_level, cufflinks_summary$intron_level, stringtie_summary$intron_level)
transcript_summary = rbind(strawberry_summary$transcript_level, cufflinks_summary$transcript_level, stringtie_summary$transcript_level)

#str(all_summary)

#Setting up the plots
## A square 2x2 plots
pdf("RD100-assembly.pdf")
###seting up plot margin###
op = par(mfrow= c(2,2),
         oma = c(3,4,2,0.3) + 0.1,
         mar = c(3,1,0,0.3) + 0.1,
         xpd = NA)

### plot first row
plot(base_summary$val, 
     base_summary$software, 
     col=mycol[base_summary$label], 
     yaxt="n", xaxt="n",
     pch=2, 
     cex=0.8,
     axes=FALSE, 
     xlab="Nucleotides detected (%)", ylab="",
     ylim=c(0.5,3.5), xlim = c(65,100), 
     yaxs = "i", xaxs = "i",
     panel.first = abline(h=c(1,2,3), xpd = FALSE, col = "cornsilk3"))

axis(2, at=c(0.5, 1, 2, 3, 3.1), 
     labels=c("", "Strawberry", "Cufflinks", "StringTie", ""), 
     lwd=2, tck=0, las=2, mgp=c(0,0.3,0))

axis(1, at=c(65,70,75,80,85,90,95,99), lwd=2)

plot(exon_summary$val, 
     exon_summary$software, 
     col=mycol[exon_summary$label], 
     yaxt="n", xaxt="n",
     pch=2, 
     cex=0.8,
     axes=FALSE, 
     xlab="Exons detected (%)", ylab="",
     ylim=c(0.5,3.5), xlim = c(60,100), 
     yaxs = "i", xaxs = "i",
     panel.first = abline(h=c(1,2,3), xpd = FALSE, col = "cornsilk3"))

axis(2, lwd=2, at=c(0.5, 1, 2, 3, 3.1), labels=FALSE, tck=0)
axis(1, at=c(60,65,70,75,80,85,90,95,99), lwd=2)

#### Print mean
library(plyr)
print(ddply(exon_summary, ~label+software, summarise, mean=mean(val)))

##### plot legent 
legend('topright',
       c('Sensitivity','Specificity'), 
       horiz =TRUE,
       pch=2,
       col = mycol[exon_summary$label], 
       bty = 'n',
       border = NA)

######## plot second row

plot(intron_summary$val, 
     intron_summary$software, 
     col=mycol[intron_summary$label], 
     yaxt="n", xaxt="n",
     pch=2, 
     cex=0.8,
     axes=FALSE, 
     xlab="Introns detected (%)", ylab="",
     ylim=c(0.5,3.5), xlim = c(65,100), 
     yaxs = "i", xaxs = "i",
     panel.first = abline(h=c(1,2,3), xpd = FALSE, col = "cornsilk3"))

axis(2, at=c(0.5, 1, 2, 3, 3.1), 
     labels=c("", "Strawberry", "Cufflinks", "StringTie", ""), 
     lwd=2, tck=0, las=2, mgp=c(0,0.3,0))

axis(1, at=c(65,70,75,80,85,90,95,99), lwd=2)

print(ddply(intron_summary, ~label+software, summarise, mean=mean(val)))

plot(transcript_summary$val, 
     transcript_summary$software, 
     col=mycol[transcript_summary$label], 
     yaxt="n", xaxt="n",
     pch=2, 
     cex=0.8,
     axes=FALSE, 
     xlab="Transcripts detected (%)", ylab="",
     ylim=c(0.5,3.5), xlim = c(30,100), 
     yaxs = "i", xaxs = "i",
     panel.first = abline(h=c(1,2,3), xpd = FALSE, col = "cornsilk3"))

axis(2, lwd=2, at=c(0.5, 1, 2, 3, 3.1), labels=FALSE, tck=0)
axis(1, at=c(30,40,50,60,70,80,90), lwd=2)

print(ddply(transcript_summary, ~label+software, summarise, mean=mean(val)))
par(op)

dev.off()


