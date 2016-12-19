########################
### Running times:
###HepG2:  
# Strawberry: 685204ms,  11.4m, 8 threads
# Samtools view: 12.35m
# wc: 8.69m
# StringTie: 4.05m, 8 threads
# StringTie: 5.5m, 1 threads
# Cufflinks: 62.2m, 8 threads

###RD100:
# Strawberry: 46014ms, 0.77m, 8 threads
# StrintTie: 0.35m, 8 threads
# Samtools view: 0.63m
# wc: 0.53m

###RD25:
# Strawberry: 14541.6ms, 0.24m, 8 threads
# StringTie: 0.1m, 8 threads
# Samtools view: 0.064m
# wc: 0.127m
###################
setwd("/home/ruolin/git/CompareTransAbun/strawberry_comp/")

library(ggplot2)
pdf("times.pdf")
dat = data.frame(
  software = factor(c(rep("Strawberry",3), rep("StringTie",3), rep("Word Count", 3), rep("Cufflinks",3) )),
  experiment = factor(c(rep(c("2.5m", "10m", "100m"),4)), levels=c("2.5m", "10m", "100m")),
  run_times = c(0.24, 0.77, 12.35, 0.1, 0.35, 4.05, 0.127, 0.53, 8.69, 0.6, 1.91, 62.2) 
  )

ggplot(data=dat, aes(x=experiment, y=run_times, group=software, shape=software, colour=software)) +
  geom_line(aes(linetype=software), size=1) + 
  geom_point(size=3, fill="white") + 
  expand_limits(y=0) +
  scale_colour_hue(name="Software", l=30) +
  scale_shape_manual(name="Software", values = c(22,21,20,19)) +
  scale_linetype_discrete(name="Software")+
  xlab("Dataset") + ylab("Total running time in minutes") +
  theme_bw() +
  theme(axis.text.x=element_text(face="bold", size=14),
        axis.text.y=element_text(face="bold", size=14),
        axis.title.x=element_text(face="bold", size=14),
        axis.title.y=element_text(face="bold", size=14),
        legend.text=element_text(face="bold", size=14),
        legend.title=element_text(face="bold",size=14)) +
  theme(legend.position=c(.3, .6))
dev.off()
