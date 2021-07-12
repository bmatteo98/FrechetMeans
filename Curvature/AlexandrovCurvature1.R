

dtr <- function (x, y){
  x = c(0,x)
  y = c(0,y)
  return (max(x-y) - min(x-y))
}
deu <- function (x, y) return (sqrt(sum((x-y)^2)))

tropicalLine <-  function (x,y, tk) {
  if (identical(x,y)) {
    mu =  (x)
    
  }
  if ((x[1] == y[1]) | (x[2] == y[2])) {
    mu = ((1-tk)*x+y*tk)
    
  }
  
  if ((x[1] == x[2]) | (y[1] == y[2])){
    mu = ((1-tk)*x+y*tk)
    
  }
  
  if (x[1]<y[1]) {
    a = x
    b = y
  }
  if (x[1]>y[1]) {
    a = y
    b = x
  }
  if ((a[2]<b[2]) & ((a[1]-a[2])<=(b[1]-b[2]))){
    t = (b[2]-a[2])/(b[1]-a[1])
    if (tk <= t) {
      mu = c(a[1]+tk*(b[1]-a[1]), a[2]+tk*(b[1]-a[1]))
    }
    if (t < tk) {
      mu = c(a[1]+tk*(b[1]-a[1]), b[2])
      
    }
  }
  
  if ((a[2]<b[2]) & ((a[1]-a[2])>(b[1]-b[2]))){
    t = (b[1]-a[1])/(b[2]-a[2])
    if (t >= tk) {
      mu = c(a[1]+tk*(b[2]-a[2]), a[2]+tk*(b[2]-a[2]))
    }
    if (t < tk) {
      mu = c(b[1],a[2]+tk*(b[2]-a[2]))
      
    }
  }
  
  if (a[2]>b[2]){
    t = (a[2]-b[2])/(a[2]-b[2]+b[1]-a[1])
    if (t >= tk) {
      mu = c(a[1],a[2]-tk*(a[2]-b[2]+b[1]-a[1]))
    }
    if (t < tk) {
      mu = c(a[1]+b[2]-a[2]+tk*(a[2]-b[2]+b[1]-a[1]),b[2])
    }
  }
  return (mu)
}

distances <- function (a,b,c,x1,y1,c1){
  t <- seq(0,1, by= 0.01)
  distan = matrix(NA, nrow = length(t),ncol = 2)
  if (x1[1]>=y1[1]) {
    a1 = y1
    b1 = x1
  }
  else {
    a1=x1
    b1=y1
  }
  for (i in 1:length(t)){
    tk = t[i]
    muTR = tropicalLine (a,b,tk)
    muE= (1-tk)*a1+b1*tk
    distan[i,] =  c(deu(muE, c1), dtr(muTR, c))
  }
  distan = round(distan , 7)
  if ((distan[1,1] ==distan[length(t),2]) & (distan[1,2] ==distan[length(t),1])){
    distan[,1] = rev(distan[,1])
  }
  return (distan)
}

findEuclidean <- function (P){
  a = P[,1]
  b = P[,2]
  c = P[,3]
  dab = dtr(a,b)
  dbc = dtr(b,c)
  dac = dtr(c,a)
  a1 = c(0,0)
  b1 = c(dab,0)
  xc = (dac^2-dbc^2+dab^2)/(2*dab)
  if (abs(round(xc, 5)) == abs(round(dac, 5))) yc = 0
  else yc = sqrt(dac^2-xc^2)
  c1 = c(xc, yc)
  P1 <- matrix(c(a1,b1, c1),nrow=length(a1))
  return (P1)
}

curvature <- function (P){
  P1 <- findEuclidean(P)
  curvS = 0
  curvF = 0
  for (j in 1:ncol(P)){
    a = P[,j%%ncol(P)+1]
    b = P[,(j+1)%%ncol(P)+1]
    c = P[,(j+2)%%ncol(P)+1]
    a1 = P1[,j%%ncol(P1)+1]
    b1 = P1[,(j+1)%%ncol(P1)+1]
    c1 = P1[,(j+2)%%ncol(P1)+1]
    distan = round(distances (a,b,c,a1,b1,c1), 5)
    if (all(distan[,1] <= distan[,2])) curvF = curvF +1 
    else if (all(distan[,1] >= distan[,2])) curvS = curvS +1
    else (return (0))
  }
  if (curvS == 3) return (-1)
  if (curvF == 3) return (1)
  else return (0)
}


# skinny triangle
a = c(1,3)
b= c(0,0)
c= c(3,2)


# fat triangle
a = c(0,2)
b = c(1,0)
c = c(3,3)

#undefined on the paper
# but it comes positive
a = c(0,2)
b = c(2,0)
c = c(3,4)


#fat
a = c(0,4)
b = c(3,0)
c = c(5, 6)


#skinny 
a = c(0,0)
b = c(2,4)
c = c(5,1)

#undefined
a = c(0,0)
b = c(448,449)
c = c(452,256)

# could there be 0 curved regions??
a = c(0,0)
b = c(1, 0)
c = c(1,1)

#P is the tropical triangle
P<- matrix(c(a,b, c),nrow=length(a))

curvature(P)

curvatures = c()
for (i in 1:100){
  a = runif(2, min = -1, max = 0)
  b = runif(2, min = -1, max = 0)
  c = runif(2, min = -1, max = 0)
  P = matrix(c(a,b, c),nrow=length(a))
  curvatures = c(curvatures, curvature(P))
}

sum(curvatures==0)/length(curvatures) #0.645 # 0.6653
sum(curvatures==1)/length(curvatures) #0.282 # 0.246
sum(curvatures==-1)/length(curvatures) #0.073 #0.0887
# same proportions in [-1, 0] and [0,1]

library(ape)
library(adephylo)
# non ultrametric
Tr = rtree(3, br = runif, min = 0, max = 2)
# ultrametric
Tr_ultra=as.phylo(as.hclust(chronos(Tr, lambda=0) ))
plot(Tr, type = "cladogram",label.offset = 0.05)
nodelabels()
tiplabels()
length(Tr$edge.length)
dists = distTips(Tr_ultra, tips = "all", method = "patristic", useC = TRUE)


projTroursPt <- function (leaves, ultra = F){
  Tr = rtree(leaves)
  if (ultra) Tr = as.phylo(as.hclust(chronos(Tr, lambda=0) ))
  Trd = distTips(Tr, tips = "all", method = "patristic", useC = TRUE)
  pt = Trd - Trd[1]
  return( pt[-1])
}
curvatures = c()
for (i in 1:100){
  isTriangle = F
  while(isTriangle != TRUE){
    a = projTroursPt(3, F)
    b = projTroursPt(3, F)
    c = projTroursPt(3, F)
    if ((identical(a,b) == F) & (identical(a,c) == FALSE) & (identical(b,c) == FALSE)) isTriangle = TRUE
  }
  P = matrix(c(a,b, c),nrow=length(a))
  curvatures = c(curvatures, curvature(P))
}
sum(curvatures==0)/length(curvatures) # 0.675
sum(curvatures==1)/length(curvatures) # 0.211
sum(curvatures==-1)/length(curvatures) # 0.114


for (i in 1:100){
  isTriangle = F
  while(isTriangle != TRUE){
    a = projTroursPt(4, T)
    b = projTroursPt(4, T)
    c = projTroursPt(4, T)
    if ((identical(a,b) == F) & (identical(a,c) == FALSE) & (identical(b,c) == FALSE)) isTriangle = TRUE
  }
  P = matrix(c(a,b, c),nrow=length(a))
  #curvatures = c(curvatures, curvature(P))
}