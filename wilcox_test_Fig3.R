#Upload tables of microbial genus/pathway abundance from different groups (CO, ABX, CON, SP and PEC), respectively
#Conducting comparisons between the groups to assess differences between them

data<-read.csv("genus.csv",header = T)
# Or data<-read.csv("pathway.csv",header = T)

n<-c(1:nrow(data))
p.w1<-rep(NA,nrow(data))
p.w2<-rep(NA,nrow(data))
p.w3<-rep(NA,nrow(data))
p.w4<-rep(NA,nrow(data))

#The columns of 2 to 5, 6 to 9, 10 to 13, 14 to 17, and 18 to 20, belong to the CO, ABX, CON, SP and PEC groups, respectively
for ( i in n){ 
  p.w1[i]<- wilcox.test(as.numeric(data[i,2:5]),as.numeric(data[i,6:9]))$p.value;
}
p.w1<-as.numeric(p.w1)
for ( i in n){ 
  p.w2[i]<- wilcox.test(as.numeric(data[i,10:13]),as.numeric(data[i,14:17]))$p.value;
}
p.w2<-as.numeric(p.w2)
for ( i in n){ 
  p.w3[i]<- wilcox.test(as.numeric(data[i,10:13]),as.numeric(data[i,18:20]))$p.value;
}
p.w3<-as.numeric(p.w3)
for ( i in n){ 
  p.w4[i]<- wilcox.test(as.numeric(data[i,14:17]),as.numeric(data[i,18:20]))$p.value;
}
p.w4<-as.numeric(p.w4)


res<-cbind(data,p.w1,p.w2,p.w3,p.w4)

dev.off()