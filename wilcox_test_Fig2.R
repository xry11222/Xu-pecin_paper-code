#Upload tables of microbial genus/pathway abundance from different days, respectively
#Conducting comparisons between the groups to assess differences between them

data<-read.csv("genus.csv",header = T)
# Or data<-read.csv("pathway.csv",header = T)

n<-c(1:nrow(data))
p.w1<-rep(NA,nrow(data))
p.w2<-rep(NA,nrow(data))
p.w3<-rep(NA,nrow(data))


#The columns of 2 to 7 belong to the CON group, 8 to 13 belong to the SP group, and 14 to 19 belong to the PEC group

for ( i in n){ 
  p.w1[i]<- wilcox.test(as.numeric(data[i,2:7]),as.numeric(data[i,8:13]))$p.value;
}
p.w1<-as.numeric(p.w1)
for ( i in n){ 
  p.w2[i]<- wilcox.test(as.numeric(data[i,2:7]),as.numeric(data[i,14:19]))$p.value;
}
p.w2<-as.numeric(p.w2)
for ( i in n){ 
  p.w3[i]<- wilcox.test(as.numeric(data[i,8:13]),as.numeric(data[i,14:19]))$p.value;
}
p.w3<-as.numeric(p.w3)



res<-cbind(data,p.w1,p.w2,p.w3)


dev.off()