# This script was created to generate backbone phylogenetic hypothesis of the Aedini tribe across various key literature
# The trees were arranged in a figure plate before further editing in adobe illustrator

rm(list=ls())
setwd(" ")
library(gridExtra); library(ggtree)

belkin<-read.tree(text="(sectionA,sectionB);")
(b.tree<-ggtree(belkin,lwd = 1)+geom_tiplab(size=8, offset=0.1, color="black")+ ggplot2::xlim(0,4))

rhk<-read.tree(text="((((cladeB,paracladeB),cladeA),mucidus),psorophora);")
(rhk.tree<-ggtree(rhk,lwd = 1)+geom_tiplab(size=8, offset=0.1, color="black")+ ggplot2::xlim(0,10))

molhyp<-read.tree(text="((cladeA,cladeB),psorophora);")
(molhyp.tree<-ggtree(molhyp,lwd = 1)+geom_tiplab(size=8, offset=0.1, color="black")+ ggplot2::xlim(0,6))

pdf("aedini_topologies.pdf", width=15, height=8)
plot_grid(b.tree, rhk.tree, molhyp.tree, labels=c("A", "B", "C"), ncol = 3, nrow = 1)
dev.off()
