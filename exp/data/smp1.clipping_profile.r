png("/sonas-hs/bsr/hpc/data/yjin/test_BAMqc/exp/test1/figs/smp1.clipping_profile.png",width=500,height=500,units="px")
read_pos=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99)
count=c(10080,9066,8130,7621,7132,6639,6173,5746,5321,4915,4597,4335,4076,3809,3556,3301,3060,2829,2599,2347,2127,1921,1710,1505,1283,1093,918,736,568,415,305,164,84,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,85,162,278,386,526,692,871,1040,1222,1412,1610,1804,2008,2226,2447,2677,2949,3192,3438,3680,3952,4204,4467,4757,5154,5574,5983,6410,6874,7383,7905,8793,9768)
plot(read_pos,1-(count/79258),pch=20,xlab="Position of reads",ylab="Mappability",col="blue")
dev.state=dev.off()
