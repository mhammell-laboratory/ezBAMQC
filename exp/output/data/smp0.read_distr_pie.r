png("/sonas-hs/bsr/hpc/data/yjin/test_BAMqc/exp/test1/figs/smp0.read_distr_pie.png",width=500,height=500,units="px")
pie(c(35107,508528),labels=c("Covered  35107 exons","Uncovered"),main="Exons",radius=0.6,clockwise=T,col=c("blue","white"))
dev.state = dev.off()
