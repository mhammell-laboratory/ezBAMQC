png("/sonas-hs/bsr/hpc/data/yjin/test_BAMqc/exp/test1/figs/smp2.readlen_profile.png",width=500,height=500,units="px")
readlen_val=c(67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100)
readlen_count=c(182,177,265,292,322,337,407,384,377,408,421,433,465,412,463,458,469,448,457,476,520,507,519,576,723,741,800,824,805,847,931,1476,1885,50582)
plot(readlen_val,(readlen_count/80603),pch=20,xlab="Mapped Read Length",ylab="Proportion",col="blue")
dev.state=dev.off()