install.packages("~/repos/cpcbp_0.3.3.tar.gz")
require(cpcbp)
setwd("~/repos/EQG2014/Exercise_2.1/")
dat <- readLines("All_females_2pops.txt")
head(dat)
head.dat <- dat[1]
dat <- dat[-1]
dat <- sapply(dat, strsplit, " ")
dat <- do.call(rbind, dat)
rownames(dat) <- NULL
dat <- data.frame(dat)
dat[dat==-9999] <- NA
dat[,2:ncol(dat)] <- apply(dat[,2:ncol(dat)], 2, as.numeric)
dat[,2] <- factor(dat[,2])
colnames(dat) <- c("population", "family", strsplit(head.dat, " ")[[1]])
head(dat)
plot(dat$body, dat$sub, col=dat$population)

GS <- list(data=dat[3:ncol(dat)], f=as.numeric(dat[,1]))

phillips.cpc(GS)


## Example from the cpcrand website
install.packages("cpca")
Males Females
HumLen HumWid FemLen FemWid
mat1 <- matrix(c(1.1544, 0.9109, 1.0330, 0.7993,
                 0.9109, 2.0381, 0.7056, 1.4083,
                 1.0330, 0.7056, 1.2100, 0.7958,
                 0.7993, 1.4083, 0.7958, 2.0277), ncol=4, byrow=TRUE)
mat2 <- matrix(c(0.9617, 0.2806, 0.9841, 0.6775,
                 0.2806, 1.8475, 0.3129, 1.2960,
                 0.9841, 0.3129, 1.2804, 0.7923,
                 0.6775, 1.2960, 0.7923, 1.7819), ncol=4, byrow=TRUE)

array <- array(data=NA, dim=c(4,4,2))
array[,,1] <- mat1
array[,,2] <- mat2

n1=92
n2=47


vignette(package=cpcbp)
vignette("cpcbp")
v = c(10,5,2,1)
corr = c(0.8,0.6,0.5,0.5)
X = simdata(npts=1000,vars=v,cors=corr)
eigen(covmat(v,cor=corr))
P1 = phillips.cpc(X)
P1$cpc1.pval
P1$evecs.CPC
phillips.cpcvec(X)
pooled.cpc(X)
