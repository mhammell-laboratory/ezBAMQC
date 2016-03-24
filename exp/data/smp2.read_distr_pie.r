png("/sonas-hs/bsr/hpc/data/yjin/test_BAMqc/exp/test1/figs/smp2.read_distr_pie.png",width=500,height=500,units="px")
pie(c(32831,510804),labels=c("Covered  32831 exons","Uncovered"),main="Exons",radius=0.6,clockwise=T,col=c("blue","white"))
dev.state = dev.off()
