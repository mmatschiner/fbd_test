

t <- read.table("~/Documents/Professional/Positions/2017_2018_Basel/projects/fbd_test/res/beast/FBD_Range1/summary/parameters.txt",header=T)
hist(t$diversification_rate, xlim=c(0,0.12), ylim=c(0,100000), main="Diversification rate", xlab="Diversification rate", breaks=seq(0,0.12,(0.12-0)/100))
hist(t$turnover, xlim=c(0,1), ylim=c(0,100000), main="Turnover", xlab="Turnover", breaks=seq(0,1,(1-0)/100))

t <- read.table("~/Documents/Professional/Positions/2017_2018_Basel/projects/fbd_test/res/beast/FBD_Range2/summary/parameters.txt",header=T)
hist(t$diversification_rate, xlim=c(0,0.12), ylim=c(0,100000), main="Diversification rate", xlab="Diversification rate", breaks=seq(0,0.12,(0.12-0)/100))
hist(t$turnover, xlim=c(0,1), ylim=c(0,100000), main="Turnover", xlab="Turnover", breaks=seq(0,1,(1-0)/100))
hist(t$sampling_proportion, xlim=c(0,0.27), ylim=c(0,100000), main="Sampling proportion", xlab="Sampling proportion", breaks=seq(0,0.27,(0.27-0)/100))

t <- read.table("~/Documents/Professional/Positions/2017_2018_Basel/projects/fbd_test/res/beast/FBD_Range3/summary/parameters.txt",header=T)
hist(t$diversification_rate, xlim=c(0,0.12), ylim=c(0,100000), main="Diversification rate", xlab="Diversification rate", breaks=seq(0,0.12,(0.12-0)/100))
hist(t$turnover, xlim=c(0,1), ylim=c(0,100000), main="Turnover", xlab="Turnover", breaks=seq(0,1,(1-0)/100))
hist(t$sampling_proportion, xlim=c(0,0.27), ylim=c(0,100000), main="Sampling proportion", xlab="Sampling proportion", breaks=seq(0,0.27,(0.27-0)/100))
