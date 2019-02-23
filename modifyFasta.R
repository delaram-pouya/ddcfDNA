library(data.table)
library("Biostrings")
library(stringr)
library(GenomicRanges)
library(IRanges)

setwd("/home/delaram/delaram/data/CareDx/CareDx_2019_Jan_08/bamsplit/results")
gnomad_genome <-fread("/media/pgdrive/apps/annovar/dls/gnom_genome_sub.txt", na.strings = ".", header = T, data.table = F)
colnames(gnomad_genome)[1] <- "Chr"

amplicons=fread("/home/delaram/delaram/data/CareDx/CareDx_2019_Jan_08/bamsplit/amplicons.bed",na.strings = '.',header=F,data.table = F)
colnames(amplicons)=c('chr','start','end')
amplicons$chr=gsub('chr','',amplicons$chr)
amplicons.GR = GRanges(seqnames=amplicons$chr,ranges=IRanges(start=amplicons$start,end=amplicons$end))

gnomad_genome2=subset(gnomad_genome, Ref%in%c('A','T','C','G')&Alt%in%c('A','T','C','G')&gnomAD_genome_ALL>0.35&gnomAD_genome_ALL<0.65)
colnames(gnomad_genome2)[1] <- "Chr"

gnom.GR = GRanges(seqnames=gnomad_genome2$Chr, IRanges(start=gnomad_genome2$Start,end=gnomad_genome2$End),
                  Ref=gnomad_genome2$Ref,Alt=gnomad_genome2$Alt,AF=gnomad_genome2$gnomAD_genome_ALL)

amp=as.data.frame(mergeByOverlaps(amplicons.GR,gnom.GR))
amp= amp[,c(1,2,3,7,14,15,16)]
colnames(amp)[1:4]=c('chr','amp.start','amp.end','gnom.pos')
amp$id=paste0(amp$chr,sep='_',amp$gnom.pos)
amp$amplicon= paste0(amp$chr,sep='_',amp$amp.start, sep='_', amp$amp.end)
amp$pos= amp$gnom.pos- amp$amp.start +1
dim(amp)
dim(subset(amp,!duplicated(id))) #532 unique positions in amplicons > AF cut-off is not enough
dim(subset(amp,!duplicated(amplicon)))
amp2 = subset(amp, pos>150 & pos< 350)
dim(amp2)
### 388 positions- cut-offs >> AF:(0.35, 0.65) - pos:(150,350) - unique amplicon names(one SNP in each position)
amp3 = subset(subset(amp2,!duplicated(amplicon)))
dim(amp3)
length(unique(amp3$amplicon))

hist(amp3$AF)
amp = amp3 
###################
s = readDNAStringSet("~/IDT405amplicons.fa")
seq_name = names(s)
sequence = paste(s)
df.fa<-data.frame(seq_name, sequence)
sp=as.data.frame(str_split_fixed(df.fa$seq_name, "_", 3))
sp$V2=as.numeric(as.character(sp$V2))
sp$V3=as.numeric(as.character(sp$V3))
colnames(sp)[1:3]=c('chr','start','end')
df.fa$chr=sp$chr
df.fa$start=sp$start
df.fa$end=sp$end
df.fa$chr=gsub('chr','',df.fa$chr)
df.fa$sequence=as.character(df.fa$sequence)
###############################
amp.df=merge(amp,df.fa,by.x=c('chr','amp.start','amp.end'),by.y=c('chr','start','end'))
sum(is.na(amp.df$sequence))
dim(amp.df)
x=amp$gnom.pos-amp$amp.start+1
amp.df$x= x
amp.df$pos=substr(amp.df$sequence,x,x)
amp.df$mod.str=amp.df$sequence
substr(amp.df$mod.str,x,x)<-amp.df$Alt
hist(amp.df$x)
##### QC  
amp.df[18,]
w= amp.df[18,'sequence']
substr(w, 201,201)
w2 = amp.df[18, 'mod.str']
substr(w2,201,201)

#######
amp.df$name=paste0(paste0('chr',amp.df$chr),sep='_',amp.df$amp.start,sep='_',amp.df$amp.end)
fasta=amp.df[,colnames(amp.df)%in%c('name','mod.str')]
colnames(fasta)=c('seq','name')
dim(fasta)

###### QC_2
check = paste0(paste0('chr',amplicons$chr), sep='_',amplicons$start,sep='_',amplicons$end)
sum(amp.df$name %in% check)

############
writeFasta<-function(data, filename){
  fastaLines = c()
  
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}
head(fasta)
fasta$name=paste0(fasta$name,'_pr')
writeFasta(fasta, "/home/delaram/delaram/data/CareDx/corrected388.fasta")

fasta2= subset(fasta,!duplicated(fasta$name))
writeFasta(fasta2, "/home/delaram/delaram/data/CareDx/modified2.fasta")

#c=data.frame(amp.df$Ref,amp.df$pos,amp.df$Alt)
#sapply(1:nrow(df.fa), function(i){sapply(1:nrow(amp),function(j){
#  if(df.fa$chr[i]==amp$chr[j]& df.fa$start[i]<=amp$gnom.pos[j] & df.fa$end[i]>=amp$gnom.pos[j])df.fa$pos[i]<<-amp$gnom.pos[j] })} )


