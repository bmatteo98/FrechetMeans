

curvatures = c()
for (P in Ps){
  #a = c(0,runif(2, min = 0, max = 10))
  #b = c(0,runif(2, min = 0, max = 10))
  #c = c(0,runif(2, min = 0, max = 10))
  #a = c(0,sample(-10:10, size = 2, replace = TRUE))
  #b = c(0,sample(-10:10, size = 2, replace = TRUE))
  #c  = c(0,sample(-10:10, size = 2, replace = TRUE))
  #P = matrix(c(a,b, c),nrow=length(a))
  #print(P)
  Pz = rbind(rep(0,3), P)
  cv = curvature(Pz)
  curvatures = c(curvatures, cv)
  print(Pz)
  plotTRN(Pz, cv)
}
plot(1,1)
sum(curvatures==0)/length(curvatures) #0.645 # 0.6653
sum(curvatures==1)/length(curvatures) #0.282 # 0.246
sum(curvatures==-1)/length(curvatures) #0.073 #0.0887
sum(curvatures=="undefined")/length(curvatures) 
# same proportions in [-1, 0] and [0,1]


curvatures = c()
for (k in 1:1000){
  #a = c(0,runif(2, min = 0, max = 1))
  #b = c(0,runif(2, min = 0, max = 1))
  #c = c(0,runif(2, min = 0, max = 1))
  a = c(0,sample(0:10, size = 2, replace = TRUE))
  b = c(0,sample(c(0:10), size = 2, replace = TRUE))
  c  = c(0,sample(c(0:10), size = 2, replace = TRUE))
  while ((identical(a,b)) | (identical(a,c)) | (identical(b,c))){ 
    a = c(0,sample(0:10, size = 2, replace = TRUE))
    b = c(0,sample(c(0:10), size = 2, replace = TRUE))
    c  = c(0,sample(c(0:10), size = 2, replace = TRUE))
  }
  
  P = matrix(c(a,b, c),nrow=length(a))
  cv = curvature(P)
  curvatures = c(curvatures, cv)
}
#plot(1,1)
sum(curvatures==0)/length(curvatures) #0.645 # 0.6653
sum(curvatures==1)/length(curvatures) #0.282 # 0.246
sum(curvatures==-1)/length(curvatures) #0.073 #0.0887
sum(curvatures=="undefined")/length(curvatures) 
# same proportions in [-1, 0] and [0,1]





curvatures = c()
for (k in 1:10){
  a = c(0,runif(2, min = 0, max = 1))
  b = c(0,runif(2, min = 0, max = 1))
  c = c(0,runif(2, min = 0, max = 1))
  
  P = matrix(c(a,b, c),nrow=length(a))
  print(P)
  cv = curvature(P)
  curvatures = c(curvatures, cv)
  plotTRN(P, cv)
}
#plot(1,1)
sum(curvatures==0)/length(curvatures) #0.645 # 0.6653
sum(curvatures==1)/length(curvatures) #0.282 # 0.246
sum(curvatures==-1)/length(curvatures) #0.073 #0.0887
sum(curvatures=="undefined")/length(curvatures) 
# same proportions in [-1, 0] and [0,1]


a = c(0,0.2973003, 0.5541335)
b = c(0,0.4744980, 0.5445846)
c = c(0,0.9128813, 0.5161453)

a= c(0,-9,3)
b= c(0,4,3)
c = c(0,4,-8)
P = matrix(c(a,b, c),nrow=length(a))
plotTR(P)

library(ape)
library(adephylo)
# non ultrametric
Tr = rtree(3, br = runif, min = 0, max = 2)
# ultrametric
Tr_ultra=as.phylo(as.hclust(chronos(Tr, lambda=0) ))
#plot(Tr, type = "cladogram",label.offset = 0.05)
#nodelabels()
#tiplabels()
#length(Tr$edge.length)
#dists = distTips(Tr_ultra, tips = "all", method = "patristic", useC = TRUE)


projTroursPt <- function (leaves, ultra = F){
  Tr = rtree(leaves)
  if (ultra) Tr = as.phylo(as.hclust(chronos(Tr, lambda=0, quiet=TRUE) ))
  Trd = distTips(Tr, tips = "all", method = "patristic", useC = TRUE)
  pt = Trd - Trd[1]
  return( pt)
}
curvatures = c()
for (i in 1:1000){
  isTriangle = F
  while(isTriangle != TRUE){
    a = round(as.numeric(projTroursPt(5, T), 8))
    b = round(as.numeric(projTroursPt(5, T), 8))
    c = round(as.numeric(projTroursPt(5, T), 8))
    if ((identical(a,b) == F) & (identical(a,c) == FALSE) & (identical(b,c) == FALSE)) isTriangle = TRUE
  }
  P = matrix(c(a,b, c),nrow=length(a))
  curvatures = c(curvatures, curvature(P))
}
sum(curvatures==0)/length(curvatures) 
sum(curvatures==1)/length(curvatures) 
sum(curvatures==-1)/length(curvatures) 


# ultrametrics 5 leaves
# undefined 0.977
# positive 0.021
# negative 0.002



# 2-dim

for (i in 1:100){
  x= runif(2, min = 0, max = 10)
  y = runif(2, min = 0, max = 10)
  z= runif(2, min = 0, max = 10)
  P = matrix(c(x,y, z),nrow=length(x))
  cv = curvature(P)
  curvatures = c(curvatures, cv)
  #print(P)
  #plotTR2(P, cv)
  #Ps[[i]] = P
}
plot(1,1)
sum(curvatures==0)/length(curvatures) #0.645 # 0.6653
sum(curvatures==1)/length(curvatures) #0.282 # 0.246
sum(curvatures==-1)/length(curvatures) #0.073 #0.0887
sum(curvatures=="undefined")/length(curvatures) 
# same proportions in [-1, 0] and [0,1]

[1] 0
> sum(curvatures==1)/length(curvatures) #0.282 # 0.246
[1] 0.2318182
> sum(curvatures==-1)/length(curvatures) #0.073 #0.0887
[1] 0.07636364
> sum(curvatures=="undefined")/length(curvatures) 
[1] 0.6918182


[1] 0
> sum(curvatures==1)/length(curvatures) #0.282 # 0.246
[1] 0.2361538
> sum(curvatures==-1)/length(curvatures) #0.073 #0.0887
[1] 0.08230769
> sum(curvatures=="undefined")/length(curvatures) 
[1] 0.6815385

