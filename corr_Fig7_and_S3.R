
library(psych)
    
r<-corr.test(as.matrix(df),use="pairwise",method = "spearman",alpha = 0.05,adjust = "none")
r.r=r$r
r.p=r$p

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
res<-flattenCorrMatrix(r$r, r$p)


dev.off()

