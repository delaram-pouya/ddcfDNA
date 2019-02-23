library(data.table)
library(parallel)
library(ggplot2)
library(ggrepel)
library('UpSetR')
library(plyr)

setwd("/home/delaram/delaram/data/CareDx/CareDx_2019_Jan_08/bamsplit/results")
target <- fread("Target_SNPs.txt", na.strings = ".", header = T, data.table = F)
All <- fread("All.txt", na.strings = ".", header = T, data.table = F)
gnomad_genome <-fread("/media/pgdrive/apps/annovar/dls/gnom_genome_sub.txt", na.strings = ".", header = T, data.table = F)
colnames(gnomad_genome)[1] <- "Chr"
split=split(target ,target$Sample)
split = lapply(split, function(i) {i$chr<-as.character(i$chr);return(i)})
tmp <- do.call(rbind, split)
tmp <- unique(tmp[,c("chr","position")])
g <- subset(gnomad_genome, Start %in% tmp$position)
tmp2 <- merge(tmp, g, by.x=c("chr","position"), by.y=c("Chr","Start"))
d <- duplicated(tmp2[,1:2])
sum(d)
head(tmp2[d,])
head(tmp2[tmp2$position %in% tmp2[d,"position"]])


merged=readRDS("/home/delaram/Gnometarget.rds")
#merged=mclapply(split,function(i) merge(i,tmp2,by.x =c("chr","position"),by.y =c("chr","position"),sort = F,all.x = T),mc.cores = 2)
all=do.call(rbind,merged)
all <- subset(all, gnomAD_genome_ALL > 0.4)
all$m1=NA;all$m2=NA
sapply(1:nrow(all), function(i){tmp=t(apply(all[i,4:7],1,sort)); all[i,]$m1<<-colnames(tmp)[4]; all[i,]$m2<<-colnames(tmp)[3]})



all2 = all[all$Alt %in%c('A','T','C','G') & all$Ref %in%c('A','T','C','G'),]
all2$alt.C <- sapply(1:nrow(all2),function(i){i=all2[i,];base=c('A','C','G','T'); l =!base%in%c(i$Ref); do.call(sum,i[base[l]]) } ) / all2$depth
all2$E <- sapply(1:nrow(all2),function(i){i=all2[i,];base=c('A','C','G','T'); l =!base%in%c(i$Alt, i$Ref); do.call(sum,i[base[l]]) } ) / all2$depth
all2$ref_error <- !(all2$m1==all2$Ref & all2$m2==all2$Alt) & !(all2$m1==all2$Alt & all2$m2==all2$Ref)
all2$c_error <- (all2$alt.C>0.05 & all2$alt.C<0.35)|(all2$alt.C>0.6 & all2$alt.C<0.9)
all2$label=sapply(1:nrow(all2), function(i){
  if(all2[i,'ref_error'] & all2[i,'c_error'])return('both')
  else if(all2[i,'ref_error']) return('ref-alt error')
  else if(all2[i,'c_error']) return('CareDx error')
  else return('non')}) 
#all2$numAlt_error <- 1-(all2$F1+all2$F2) > 0.1
as.data.frame(table(all2$label))

all3 <- all2[ !(all2$m1==all2$Ref & all2$m2==all2$Alt) & !(all2$m1==all2$Alt & all2$m2==all2$Ref) ,]
all4 <- all2[ (all2$m1==all2$Ref & all2$m2==all2$Alt) | (all2$m1==all2$Alt & all2$m2==all2$Ref) ,]
all4$Sample= as.character(all4$Sample)
all3$Sample=as.character(all3$Sample)

all3$E <- sapply(1:nrow(all3),function(i){i=all3[i,];base=c('A','C','G','T'); l =!base%in%c(i$Alt, i$Ref); do.call(sum,i[base[l]]) } ) / all3$depth
all4$E <- sapply(1:nrow(all4),function(i){i=all4[i,];base=c('A','C','G','T'); l =!base%in%c(i$Alt, i$Ref); do.call(sum,i[base[l]]) } ) / all4$depth
ggplot(all3, aes(E+1e-5, color=Sample))+geom_density()+scale_x_log10()+theme_bw()

ggplot(all4, aes(E+1e-5, color=Sample))+geom_density()+theme_bw()
ggplot(all4, aes(E+1e-5, color=Sample))+geom_density()+scale_x_log10()+theme_bw()
ggplot(all4, aes(y=F2+1e-5, x=Sample))+scale_y_log10()+geom_violin(aes(fill=Sample))+geom_boxplot(width=0.17)+theme_bw()
ggplot(all4, aes(y=F2+1e-5, x=Sample))+scale_y_log10()+geom_boxplot(width=0.7,aes(fill=Sample))+theme_bw()
ggplot(all4, aes(F2+1e-5, color=Sample))+geom_density()+scale_x_log10()+theme_bw()

ggplot(all3, aes(F2+1e-5, color=Sample))+geom_density()+theme_bw()
ggplot(all3, aes(F2+1e-5, color=Sample))+geom_density()+scale_x_log10()+theme_bw()


summary(table(all4$Sample))
class(all4$Sample)
lapply(all3,function(i) median(i$F2))


######################################################
pdf('/home/delaram/delaram/plot_alter.pdf')

tmp <- target[order(target$F2),]
tmp.l <- split(tmp, f = tmp$Sample)
tmp.l = lapply(tmp.l, function(i){ i$index=1:nrow(i); i})
tmp=do.call(rbind,tmp.l)
tmp$Sample = as.character(tmp$Sample)

ggplot(tmp, aes(x=index,y=ALT.freq, color=Sample))+geom_point(size=0.9)+theme_bw()+scale_color_brewer(palette="Accent")

p=ggplot(tmp, aes(x=index,y=F2, color=Sample))+geom_point(size=0.9)+theme_bw()+scale_color_brewer(palette="Accent")
p+geom_text(aes(label=ifelse(F2>0.2&F2<0.3,as.character(chr),'')),colour='black',position=position_dodge(width =0.3),vjust=-0.4)
p+geom_label_repel(aes(label=ifelse(F2>0.2&F2<0.3,as.character(chr),'')),box.padding=0.35,point.padding = 0.5,segment.color = 'grey50')

p+geom_text_repel(aes(label=ifelse(F2>0.2&F2<0.3,as.character(chr),'')),box.padding=unit(0.45,'lines'),segment.color = 'grey50')

dev.off()
#scale_color_brewer(palette="Accent")+theme_minimal()

########################################################
pdf('/home/delaram/delaram/C.alt.plot.pdf',width=11, height = 6)
t =all2[order(all2$alt.C),]
t.l =split(t, f = t$Sample)
t.l = lapply(t.l, function(i){ i$index=1:nrow(i); i})
t=do.call(rbind,t.l)
t$Sample = as.character(t$Sample)
ggplot(t, aes(x=index,y=alt.C, color=Sample))+geom_point(size=0.9)+theme_bw()+scale_color_brewer(palette="Accent")
t.l= split(t, f = t$Sample)
lapply(t.l, function(t){ggplot(t, aes(x=index,y=alt.C, color=label))+geom_point(size=0.9)+
    theme_bw()+ggtitle(t[1,'Sample'])+
    geom_text(hjust = 0, nudge_x = 0.05,check_overlap = FALSE,size=1.6,colour='black',
              aes(label=ifelse(ref_error,'R','')))
})
l = lapply(t.l, function(i) i[(i$alt.C>0.05&i$alt.C<0.35)| (i$alt.C>0.6&i$alt.C<0.9),])
l = lapply(l, function(i) sapply(1:nrow(i), function(j) paste(i[j,'chr'],i[j,'position'],sep='_',collapse = NULL)))
uni = unique(unlist(l))
tab = as.data.frame(uni)
tab[,2:9] =sapply(1:length(l),function(j) {sapply(1:nrow(tab), function(i) ifelse( tab[i,'uni']%in% l[[j]],1, 0))})
colnames(tab)=c('posID','11148','11149','11150','11151','11152','11153','11154','11155')

upset(tab,sets = c('11148','11149','11150','11151','11152','11153','11154','11155') , 
      matrix.color = "#990000",
      sets.bar.color = c('darkslategrey','cyan4', 'cyan4','cyan3','cyan3','cyan2','cyan2' ,'cyan1'), 
      sets.x.label = c('(alt.C>0.05 & alt.C<0.35)|(alt.C>0.6 & alt.C<0.9)'),
      keep.order = F,nintersects=NA, point.size = 2.6,line.size = 0.7)

dev.off()


ggplot(t.l[[1]], aes(x=index,y=alt.C, color=label))+geom_point(size=0.9)+
  theme_bw()+ggtitle(t[1,'Sample'])+
  geom_text(hjust = 0, nudge_x = 0.05,check_overlap = FALSE,size=1.6,colour='black',
            aes(label=ifelse(ref_error,'*','')))



#(alt.C>0.05&alt.C<0.35)|(alt.C>0.6&alt.C<0.9),as.character(paste0(chr,sep="_",position)),''
#################
all2$Sample=as.character(all2$Sample)
all2.tmp=subset(all2,!Sample%in%c('11148','11155'))
all2.tmp$ref_error=ifelse(all2.tmp$ref_error,1,0)
all2.tmp$c_error=ifelse(all2.tmp$c_error,1,0)
all2.tmp= all2.tmp[,c(1,2,3,18,19,20,21)]
all2.l=split.data.frame(all2.tmp,f=all2.tmp$Sample,drop=FALSE)

n=lapply(all2.l, function(i){colnames(i)[3:7]=paste(colnames(i)[3:7],sep='_',i[1,'Sample'],collapse=NULL)})
sapply(1:length(all2.l), function(i)colnames(all2.l[[i]])[3:7]<<-n[[i]] )
sapply(1:length(all2.l),function(i) all2.l[[i]]$name<<-paste(all2.l[[i]]$chr,sep='_',all2.l[[i]]$position,collapse=NULL))
sapply(1:length(all2.l), function(i) all2.l[[i]]<<-all2.l[[i]][,-c(1,2,3)])
sapply(1:length(all2.l),function(i)rownames(all2.l[[i]])<<-all2.l[[i]]$name)
name=sapply(1:length(all2.l),function(i)all2.l[[i]]$name)
ints=Reduce(intersect,name)
sapply(1:length(all2.l), function(i) all2.l[[i]]<<-all2.l[[i]][,colnames(all2.l[[i]])!='name'])

df=do.call(cbind,lapply(1:length(all2.l), function(i) all2.l[[i]][rownames(all2.l[[i]])%in%ints,]))
df$ref_error.Sum=rowSums(df[,colnames(df)%in%paste0('ref_error_111',seq(from=49,to=54,by=1))])
df$c_error.Sum=rowSums(df[,colnames(df)%in%paste0('c_error_111',seq(from=49,to=54,by=1))])
saveRDS(object = df, file = '/home/delaram/delaram/data/CareDx/errorTable.rds')
e=subset(df,ref_error.Sum>4)
ec=subset(df,c_error.Sum>3)
pdf('oneFastaPlot.pdf')
ggplot(df, aes(x=ref_error.Sum, y=c_error.Sum))+geom_point(alpha = 1/5,colour='blue')#geom_bin2d(bins=8)

tmp <- ddply(df, .(ref_error.Sum, c_error.Sum), summarise, num=length(ref_error.Sum))
ggplot(tmp, aes(ref_error.Sum, c_error.Sum))+geom_point(aes(size=num))+scale_size_continuous(trans="log2")

par(mfrow=c(1,2))
barplot(table(df$c_error.Sum),main = 'c_error',col='red')
barplot(table(df$ref_error.Sum),main='ref_error',col='blue')
par(mfrow=c(1,1))
e.df=df[,colnames(df)%in%paste0('E_111',seq(from=49,to=54,by=1))]
thr=1
e.df=subset(e.df,E_11149<thr&E_11150<thr&E_11151<thr&E_11152<thr&E_11153<thr&E_11154<thr)
plot(e.df)
alt.c.df=df[,colnames(df)%in%paste0('alt.C_111',seq(from=49,to=54,by=1))]
plot(alt.c.df,col='dark blue')
dev.off()



