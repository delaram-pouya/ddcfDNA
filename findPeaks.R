library(plyr)
x <- read.delim("CareDx.hg38.depth",header=F)
x <- na.omit(x)
y <- subset(x, V3 > 100)
y$diff <- c(1, diff(y$V2))
y$difff <- y$diff != 1
y$peakid <- cumsum(y$difff)
z <- ddply(y, .(peakid), summarise, chr=unique(V1), start=min(V2), end=max(V2), meanCoverage=mean(V3))
z <- z[,c("chr","start","end","peakid","meanCoverage")]
write.table(z, file="peaks.bed", row.names=F, col.names=T, sep="\t", quote=F)
