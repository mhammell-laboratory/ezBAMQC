png("/sonas-hs/bsr/hpc/data/yjin/test_BAMqc/exp/test1/figs/smp1.readlen_profile.png",width=500,height=500,units="px")
readlen_val=c(67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100)
readlen_count=c(182,173,267,238,308,340,386,354,392,415,421,399,419,456,470,469,496,483,489,488,529,499,514,561,748,813,751,827,863,912,958,1648,1875,47529)
plot(readlen_val,(readlen_count/79258),pch=20,xlab="Mapped Read Length",ylab="Proportion",col="blue")
dev.state=dev.off()
