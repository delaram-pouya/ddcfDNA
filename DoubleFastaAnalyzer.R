setwd("/home/delaram/delaram/data/CareDx/CareDx_2019_Feb_12/SplitBamFiles")
library(stringr)
library(data.table)
library(GenomicRanges)
library(IRanges)
library(rlang)
library(ggplot2)
library(Rsamtools)
library(BiocManager)
library(parallel)
library(plyr)
library(UpSetR)

############################################## functions 
pileupFreq <- function(pileupres) {
  nucleotides <- levels(pileupres$nucleotide)
  res <- split(pileupres, pileupres$seqnames)
  res <- lapply(res, function (x) {split(x, x$pos)})
  res <- lapply(res, function (positionsplit) {
    nuctab <- lapply(positionsplit, function(each) {
      chr = as.character(unique(each$seqnames))
      pos = as.character(unique(each$pos))
      tablecounts <- sapply(nucleotides, function (n) {sum(each$count[each$nucleotide == n])})
      c(chr,pos, tablecounts)
    })
    nuctab <- data.frame(do.call("rbind", nuctab),stringsAsFactors=F)
    rownames(nuctab) <- NULL
    nuctab
  })
  res <- data.frame(do.call("rbind", res),stringsAsFactors=F)
  rownames(res) <- NULL
  colnames(res) <- c("seqnames","start",levels(pileupres$nucleotide))
  res[3:ncol(res)] <- apply(res[3:ncol(res)], 2, as.numeric)
  res}


clean.df <- function(df){
  df = df[,-(8:10)] 
  df$depth = df$A + df$C + df$G + df$T + df$N
  colnames(df) = c("chr" ,"position", "A" , "C" , "G" , "T", "N", "depth" )
  tmp <- t(apply(df[,3:7], 1, sort))
  df$F1 <- tmp[,5]/df$depth # Maximum frequency allele
  df$F2 <- tmp[,4]/df$depth # Second maximum frequency allele
  df$ALT.freq = (rowSums(tmp[,-5]))/df$depth 
  return(df)}

Clean<- function(df){
  df$depth = df$A + df$C + df$G + df$T
  tmp <- t(apply(df[,3:6], 1, sort))
  df$F1 <- tmp[,4]/df$depth # Maximum frequency allele
  df$F2 <- tmp[,3]/df$depth # Second maximum frequency allele
  df$ALT.freq = (rowSums(tmp[,-4]))/df$depth 
  return(df)}

####################### gnomad
gnomad_genome <-fread("/media/pgdrive/apps/annovar/dls/gnom_genome_sub.txt", na.strings = ".", header = T, data.table = F)
colnames(gnomad_genome)[1] <- "Chr"
gnomad_genome2=subset(gnomad_genome, Ref%in%c('A','T','C','G')&Alt%in%c('A','T','C','G'))
colnames(gnomad_genome2)[1]='Chr'
######################### amplicon info 

amplicons=fread("/home/delaram/amplicons.bed",na.strings = '.',header=F,data.table = F)
colnames(amplicons)=c('chr','start','end')
amplicons$chr=gsub('chr','',amplicons$chr)
amplicons.GR = GRanges(seqnames=amplicons$chr,ranges=IRanges(start=amplicons$start,end=amplicons$end))
gnom.GR = GRanges(seqnames=gnomad_genome2$Chr, IRanges(start=gnomad_genome2$Start,end=gnomad_genome2$End),
                  Ref=gnomad_genome2$Ref,Alt=gnomad_genome2$Alt,AF=gnomad_genome2$gnomAD_genome_ALL)

amp=as.data.frame(mergeByOverlaps(amplicons.GR,gnom.GR))
amp=amp[,c(1,2,3,7,14,15,16)]
colnames(amp)[1:4]=c('chr','amp.start','amp.end','pos')
amp$chr=paste0('chr',amp$chr)
ori.amp=amp
sum(duplicated(paste0(ori.amp$pos,sep='_',ori.amp$chr))) #522 pos> duplicated 
hist(ori.amp[(duplicated(paste0(ori.amp$pos,sep='_',ori.amp$chr))),]$AF) # many with very low AF 
summary(ori.amp[(duplicated(paste0(ori.amp$pos,sep='_',ori.amp$chr))),]$AF) # slight cut-off on AF > 0.01
## amp object contains all the positions not only targeted ones
amp=amp[!duplicated(paste0(amp$pos,sep='_',amp$chr))|amp$AF>0.01,] # remove duplicates on a same pos with very low AF

############################# pileup section 

p_param = PileupParam(max_depth = 50000, min_base_quality=10, min_mapq=5,
                      min_nucleotide_depth=0, min_minor_allele_depth=0,
                      distinguish_strands=TRUE, distinguish_nucleotides=TRUE,
                      ignore_query_Ns=TRUE, include_deletions=TRUE, include_insertions=FALSE,
                      left_bins=NULL, query_bins=NULL, cycle_bins=NULL)

l <- list.files(".",".bam$")
samples <- sub(".*_","", sub(".bam","",l))
bf <- mclapply(l, function(i) BamFile(i, index=paste0(i, ".bai")), mc.cores=detectCores()-2)
names(bf)<-samples
bed.pileup <- mclapply(bf, pileup, pileupParam=p_param, mc.core=detectCores()-2)
freq <- mclapply(bed.pileup, pileupFreq, mc.cores=detectCores()-2) 
names(freq) <- paste(names(bf), ".freq", sep= "")
freq <- mclapply(freq, clean.df, mc.cores = detectCores()-2)
sapply(1:length(freq), function(i)colnames(freq[[i]])[1]<<-'Contig')
tmp <- mclapply(freq, function(i) {x=data.frame(str_split_fixed(i$Contig, "_", 3));colnames(x)=c('chr','start','end');return(x)})
freq <- sapply(1:length(freq), function(i)cbind(tmp[[i]],freq[[i]]),simplify =F)
sapply(1:length(freq), function(i)freq[[i]]$T.pos<<-as.numeric(as.character(freq[[i]]$start))+as.numeric(as.character(freq[[i]]$position))-1)

## makes the table for all the positions in the bam file not only targeted positions 
m <- sapply(1:length(freq), function(i)merge(freq[[i]],amp,by.x =c('chr', 'T.pos'),by.y = c('chr', 'pos'),all.x=T,all.y=F,sort=F),simplify=F)
names(m)=samples

m <- sapply(1:length(m),function(i) m[[i]][,!colnames(m[[i]]) %in% c('amp.start','amp.end')],simplify = F)
sapply(1:length(m),function(i){x=as.data.frame(str_split_fixed(m[[i]]$Contig, "_", 4));m[[i]]$fasta<<-x$V4})
sapply(1:length(m),function(i){ m[[i]]$fasta <<-substr(as.character(m[[i]]$fasta),1,2)})
saveRDS(m,'/media/pgdrive/users/delaram/data/DoubleFastaTab2.rds')

m = readRDS('/media/pgdrive/users/delaram/data/DoubleFastaTab2.rds')
############  subset m based on AF + position

m1 = lapply(m, subset, AF>0.4 & AF<0.6) ## filters NA as well 
lapply(m1, nrow)
m1 = lapply(m1, subset, position> 150 & position< 250)
lapply(m1, nrow)
test = lapply(1:length(m1), function(i) unique.data.frame(m1[[i]][,c('chr','T.pos')]))
lapply(test, nrow) ## ~370 less but better


################################ ddply + add attribute 
m.s = lapply(m1, subset, select=c('chr','T.pos','A','T','C','G'))
#lapply(1:length(m.s), function(i) m.s[[i]]$chr<<-gsub('chr','', m.s[[i]]$chr))
lapply(1:length(m.s), function(i){m.s[[i]]$T.pos<<-as.character(m.s[[i]]$T.pos); m.s[[i]]$chr<<-as.character(m.s[[i]]$chr)})
lapply(m.s, nrow)
m.f = lapply(m.s,ddply, .(chr, T.pos), numcolwise(sum))
lapply(m.f, nrow) ## same as unique >> correct
lapply(1:length(m.f), function(i) m.f[[i]]$Sample<<-samples[i])
data=do.call(rbind, m.f)
# before merge: total #pos=2271       after merge: each sample #pos=2621 >> multiple ALt
## first filtering on AF+pos > better to add unique amplicon + depth 
data.f=merge(data,amp,by.x =c('chr', 'T.pos'),by.y = c('chr', 'pos'),all.x=T,all.y=F,sort=F)  
paste0(nrow(data.f), sep=' ', sum(table(data.f$Alt))) #### NO NA value(all pos in gnomad)
data.f= Clean(data.f)

data.f$m1=NA; data.f$m2=NA
sapply(1:nrow(data.f), function(i){
  tmp=t(apply(data.f[i,3:6],1,sort)) 
  data.f[i,]$m1<<-colnames(tmp)[4]
  ifelse(data.f[i,'ALT.freq']>0, data.f[i,]$m2<<-colnames(tmp)[3] , data.f[i,]$m2<<-NA)})

data.f$alt.C <- sapply(1:nrow(data.f),function(i){i=data.f[i,];base=c('A','C','G','T'); l =!base%in%c(i$Ref); do.call(sum,i[base[l]]) } ) / data.f$depth
data.f$E <- sapply(1:nrow(data.f),function(i){i=data.f[i,];base=c('A','C','G','T'); l =!base%in%c(i$Alt, i$Ref); do.call(sum,i[base[l]]) } ) / data.f$depth
data.f$ref_error <- !(data.f$m1==data.f$Ref & (data.f$m2==data.f$Alt|is.na(data.f$m2)) )& !(data.f$m1==data.f$Alt & (data.f$m2==data.f$Ref|is.na(data.f$m2) ) )
saveRDS(data.f,'/media/pgdrive/users/delaram/data/dFasta_finalTAB.rds')
data.f = readRDS('/media/pgdrive/users/delaram/data/dFasta_finalTAB.rds')

data.f=subset(data.f, AF>0.4 & AF<0.6 ) # merged again> remove low AF multiple ALTs
data.f$c_error <- (data.f$alt.C>0.05 & data.f$alt.C<0.4)|(data.f$alt.C>0.6 & data.f$alt.C<0.95)
data.f$label=sapply(1:nrow(data.f), function(i){
  if(data.f[i,'ref_error'] & data.f[i,'c_error'])return('both')
  else if(data.f[i,'ref_error']) return('ref-alt error')
  else if(data.f[i,'c_error']) return('CareDx error')
  else return('non')}) 

as.data.frame(table(data.f$label))
data.f.er <- data.f[ !(data.f$m1==data.f$Ref & (data.f$m2==data.f$Alt|is.na(data.f$m2))) & !(data.f$m1==data.f$Alt & (data.f$m2==data.f$Ref|is.na(data.f$m2))) ,]
data.f.cor <- data.f[ (data.f$m1==data.f$Ref & (data.f$m2==data.f$Alt| is.na(data.f$m2))) | (data.f$m1==data.f$Alt & (data.f$m2==data.f$Ref|is.na(data.f$m2))) ,]


pdf('doubleFastaPlot.pdf')
ggplot(data.f, aes(F2+1e-5, color=Sample))+geom_density()+theme_bw()+ggtitle('total data')
ggplot(data.f, aes(F2+1e-5, color=Sample))+geom_density()+scale_x_log10()+theme_bw()+ggtitle('total data(log)')
ggplot(data.f.cor, aes(F2+1e-5, color=Sample))+geom_density()+scale_x_log10()+theme_bw()+ggtitle('non.ref.error(log) subset')
ggplot(data.f.er, aes(F2+1e-5, color=Sample))+geom_density()+scale_x_log10()+theme_bw()+ggtitle('ref.error(log) subset')
ggplot(data.f, aes(y=F2+1e-5, x=Sample))+scale_y_log10()+geom_violin(aes(fill=Sample))+geom_boxplot(width=0.17)+theme_bw()
ggplot(data.f, aes(y=F2+1e-5, x=Sample))+scale_y_log10()+geom_boxplot(width=0.7,aes(fill=Sample))+theme_bw()

ggplot(data.f, aes(E+1e-5, color=Sample))+geom_density()+theme_bw()+ggtitle('total data')
ggplot(data.f, aes(E+1e-5, color=Sample))+geom_density()+scale_x_log10()+theme_bw()+ggtitle('total data(log)')


df.sub = split(data.f,data.f$Sample) ## many alternates > some with low AF
lapply(df.sub, nrow) ## number of targeted positions ~ 380

df.sub = lapply(df.sub,  function(i) i[order(i$alt.C),])
df.sub = lapply(df.sub, function(i){ i$index=1:nrow(i); i})
pdf('depth.filter.pdf')
lapply(df.sub, function(t){ggplot(t, aes(x=index,y=alt.C, color=label))+geom_point(size=0.9)+
    theme_bw()+ggtitle(t[1,'Sample'])+
    geom_text(hjust = 0, nudge_x = 0.05,check_overlap = TRUE,size=1.6,colour='black',
              aes(label=ifelse(c_error,paste0(chr,sep='_',amp.start,sep='_',amp.end),'')))
})

dev.off()
data.f=do.call(rbind,df.sub)
data.f$Sample = as.character(data.f$Sample)
ggplot(data.f, aes(x=index,y=alt.C, color=Sample))+geom_point(size=0.9)+theme_bw()+scale_color_brewer(palette="Accent")+ggtitle('total data')
ggplot(data.f[!data.f$ref_error,], aes(x=index,y=alt.C, color=Sample))+geom_point(size=0.9)+theme_bw()+scale_color_brewer(palette="Accent")+ggtitle('non.ref.error subset')#Not helpful


l = lapply(df.sub, function(i) i[i$c_error,])
l = lapply(l, function(i) sapply(1:nrow(i), function(j) paste(i[j,'chr'],i[j,'T.pos'],sep='_',collapse = NULL)))
uni = unique(unlist(l))
tab = as.data.frame(uni)
tab[,2:7] =sapply(1:length(l),function(j) {sapply(1:nrow(tab), function(i) ifelse( tab[i,'uni']%in% l[[j]],1, 0))})
colnames(tab)[1]=c('posID')
colnames(tab)[2:7]= names(l)

upset(tab,sets = names(l) , 
      matrix.color = "#990000",mainbar.y.label = 'CareDx error intersection',
      sets.bar.color = c('darkslategrey','cyan4','cyan3','cyan3','cyan2','cyan1'), 
      sets.x.label = c('(alt.C>0.05 & alt.C<0.4)|(alt.C>0.6 & alt.C<0.95)'),
      keep.order = F,nintersects=NA, point.size = 2.6,line.size = 0.7)

#################
data.f$ref_error=ifelse(data.f$ref_error,1,0)
data.f$c_error=ifelse(data.f$c_error,1,0)
tmp= data.f[,c(1,2,7,19,20,21,22)]
tmp=split.data.frame(tmp,f=tmp$Sample,drop=FALSE)

n=lapply(tmp, function(i){colnames(i)[4:7]=paste(colnames(i)[4:7],sep='_',i[1,'Sample'],collapse=NULL)})
sapply(1:length(tmp), function(i)colnames(tmp[[i]])[4:7]<<-n[[i]] )
sapply(1:length(tmp),function(i) tmp[[i]]$name<<-paste(tmp[[i]]$chr,sep='_',tmp[[i]]$T.pos,collapse=NULL))
sapply(1:length(tmp), function(i) tmp[[i]]<<-tmp[[i]][,-c(1,2,3)])
sapply(1:length(tmp),function(i)rownames(tmp[[i]])<<-tmp[[i]]$name)
name=sapply(1:length(tmp),function(i)tmp[[i]]$name)
ints=Reduce(intersect,name) # 364 positiosn present in all of the samples 
sapply(1:length(tmp), function(i) tmp[[i]]<<-tmp[[i]][,colnames(tmp[[i]])!='name'])

tmp=do.call(cbind,lapply(1:length(tmp), function(i) tmp[[i]][rownames(tmp[[i]])%in%ints,]))
tmp$ref_error.Sum=rowSums(tmp[,colnames(tmp)%in%paste0('ref_error_',names(l))])
tmp$c_error.Sum=rowSums(tmp[,colnames(tmp)%in%paste0('c_error_',names(l))])
#saveRDS(object = tmp, file = '/home/delaram/delaram/data/CareDx/DoubleFastaErrorTab.rds')
e=subset(tmp,ref_error.Sum>0)
ec=subset(tmp,c_error.Sum>0)
saveRDS(ec,file= '/home/delaram/delaram/data/CareDx/careDxErrors.rds')
ec = readRDS('/home/delaram/delaram/data/CareDx/careDxErrors.rds')
ggplot(tmp, aes(x=ref_error.Sum, y=c_error.Sum))+geom_point(alpha = 1/5,colour='blue')#geom_bin2d(bins=8)

tmp.p <- ddply(tmp, .(ref_error.Sum, c_error.Sum), summarise, num=length(ref_error.Sum))
ggplot(tmp.p, aes(ref_error.Sum, c_error.Sum))+geom_point(aes(size=num))+scale_size_continuous(trans="log2")

par(mfrow=c(1,2))
barplot(table(tmp$c_error.Sum),main = 'c_error',col='red')
barplot(table(tmp$ref_error.Sum),main='ref_error',col='blue')
par(mfrow=c(1,1))
e.tmp=tmp[,colnames(tmp)%in%paste0('E_',names(l))]
thr=1
e.tmp=subset(e.tmp,E_12707<thr&E_12708<thr&E_12709<thr&E_12710<thr&E_12711<thr&E_12712<thr)
plot(e.tmp)
alt.c.tmp=tmp[,colnames(tmp)%in%paste0('alt.C_',names(l))]
plot(alt.c.tmp,col='dark blue')


dev.off()

####### QC
nrow(amplicons)*2  ## total number of initial Contigs
lapply(1:length(bed.pileup), function(i) length(unique(bed.pileup[[i]]$seqnames))) ## number of unique contigs in bam file
lapply(1:length(bed.pileup), function(i) nrow(amplicons)*2 -length(unique(bed.pileup[[i]]$seqnames))) ##contigs with 0 aligned read
lapply(1:length(freq), function(i) length(unique(freq[[i]]$Contig))) ## same num as above
lapply(1:length(freq), function(i) nrow(freq[[i]])/length(unique(freq[[i]]$Contig))) ## ~90 positions for each contig 
lapply(1:length(freq), function(i) table(freq[[i]]$Contig) )  ## exact 
lapply(freq, function(i) sum(i$depth==0)) ## num positions with zero depth(not removed)
lapply(1:length(freq), function(i) table(round(freq[[i]]$ALT.freq,2)))  ## many ~0 ALT.freq

lapply(1:length(freq), function(i) paste0(nrow(freq[[i]]),sep=' ',nrow(m[[i]]))) # more pos after merge> many ALT positions in gnomad
sapply(1:length(m),function(i) table(as.character(m[[i]]$fasta))/nrow(m[[i]]))#ratio of aligned read to 2 ref files

lapply(1:length(m), function(i) range(as.numeric(m[[i]]$position)) ) #range of snp positions in amplicons 

test = lapply(m, subset, AF>0.4 & AF<0.6)
lapply(test, dim)
test2 = lapply(test, subset, position> 150 & position< 250)
lapply(test2, dim)
test3 = lapply(1:length(test2), function(i) unique.data.frame(test2[[i]][,c('chr','T.pos')]))
lapply(test3, dim)



######## sequencing errors: 
# subset(df1, E==1 & depth>500)
#########################


