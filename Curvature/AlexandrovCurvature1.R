

dtr <- function (x, y){
  x = c(0,x)
  y = c(0,y)
  return (max(x-y) - min(x-y))
}
deu <- function (x, y) return (sqrt(sum((x-y)^2)))

tropicalLine <-  function (x,y, tk) {
  if (identical(x,y)) {
    mu =  (x)
    return (mu)
  }
  if ((x[1] == y[1]) | (x[2] == y[2])) {
    mu = ((1-tk)*x+y*tk)
    return (mu)
  }
  
  if ((x[1] == x[2]) & (y[1] == y[2])){
    mu = ((1-tk)*x+y*tk)
    return (mu)
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

distances <- function (x,y,z,x1,y1,c1){
  t <- seq(0,1, length.out = 200)
  distan = matrix(NA, nrow = length(t),ncol = 2)
  if (x1[1]>=y1[1]) {
    a1 = y1
    b1 = x1
  }
  else {
    a1=x1
    b1=y1
  }
  EUsegment = c()
  for (i in 1:length(t)){
    tk = t[i]

    muTR = tropicalLine (x,y,tk)

    muE= (1-tk)*a1+b1*tk
    EUsegment = rbind(EUsegment, muE)
    distan[i,] =  c( dtr(muTR, z), deu(muE, c1))
  }
  distan = round(distan , 8)
  if ((distan[1,1] ==distan[length(t),2]) & (distan[1,2] ==distan[length(t),1])){
    distan[,1] = rev(distan[,1])
  }
  return (distan)
}

#dist2 = distances(b[-1],c[-1],a[-1],P1[,2],P1[,3],P1[,1])

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
  if (abs(round(xc, 8)) == abs(round(dac, 8))) yc = 0
  else yc = sqrt(dac^2-xc^2)
  c1 = c(xc, yc)
  P1 <- matrix(c(a1,b1, c1),nrow=length(a1))
  return (P1)
}


curvature <- function (P){
  P1 <- findEuclidean(P)
  curvS = 0
  curvF = 0
  sameDist = 0
  perm = matrix(c(1,2,3,2,3,1,1,3,2), nrow = 3, byrow = TRUE)
  for (j in 1:3){
    one = perm[j,1]
    two = perm[j,2]
    three = perm[j,3]
    a = P[,one]
    b = P[,two]
    c = P[,three]
    a1 = P1[,one]
    b1 = P1[,two]
    c1 = P1[,three]
    distan = round(distances (a,b,c,a1,b1,c1), 8)
    if (all(distan[,1] == distan[,2])) sameDist = sameDist +1
    else if (all(distan[,1] <= distan[,2])) curvS = curvS +1 
    else if (all(distan[,1] >= distan[,2])) curvF = curvF +1
    else (return ("undefined"))
  }
  if (sameDist == 3) return (0)
  else if ((sameDist+curvS) == 3) return (-1)
  else if ((sameDist+curvF) == 3) return (1)
  else return ("undefined")
}
plotTR2 <- function (P, cv){
  x = P[,1]
  y = P[,2]
  z = P[,3]
  xl = c(min(x[1], y[1], z[1]),max(x[1], y[1], z[1]) )
  yl = c(min(x[2], y[2], z[2]),max(x[2], y[2], z[2]) )
  trLine = c()
  t = seq(0,1, length.out = 200)
  for (ti in t){
    trLine = rbind(trLine, tropicalLine(x,y,ti))
  }
  par(mfrow=c(2,2), mar = rep(2, 4))
  plot(trLine[,1], trLine[,2], type = 'l', col = 'red', xlim = xl, ylim = yl, main = paste(cv, "2-dim"))
  trLine = c()
  for (ti in t){
    trLine = rbind(trLine, tropicalLine(y,z,ti))
  }
  lines(trLine[,1], trLine[,2], type = 'l', col = 'blue' )
  trLine = c()
  for (ti in t){
    trLine = rbind(trLine, tropicalLine(x,z,ti))
  }
  lines(trLine[,1], trLine[,2], type = 'l', col = 'green')
  points(x[1], x[2], pch = 16)
  points(y[1], y[2], pch = 16)
  points(z[1], z[2], pch = 16)
  P1 = findEuclidean(P)
  distan = distances (x,y,z,P1[,1],P1[,2],P1[, 3])
  yl = c(min(distan[,2], distan[,1]), max(distan[,2], distan[,1]))
  plot(distan[,2], type = 'l', col = 'blue', main = 'c from a-b', ylim = yl)
  lines( distan[,1], type = 'l', col = 'red')
  distan = distances (x,z,y,P1[,1],P1[,3],P1[,2])
  yl = c(min(distan[,2], distan[,1]), max(distan[,2], distan[,1]))
  plot(distan[,2], type = 'l', col = 'blue', main = 'b from a-c', ylim = yl)
  lines( distan[,1], type = 'l', col = 'red')
  distan = distances (y,z,x,P1[,2],P1[,3],P1[,1])
  yl = c(min(distan[,2], distan[,1]), max(distan[,2], distan[,1]))
  plot(distan[,2], type = 'l', col = 'blue', main = 'a from b-c', ylim = yl)
  lines( distan[,1], type = 'l', col = 'red')
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
P<- matrix(c(a[-1],b[-1], c[-1]),ncol=length(a))
P = matrix(c(a,b, c),nrow=length(a))
curvature(P)


P1 <- findEuclidean(P)
a1 = P1[,1]
b1= P1[, 2]
c1 = P1[, 3]
distan2 = distances (a[-1],b[-1],c[-1],P1[,1],P1[, 2],P1[, 3])
t <- seq(0,1, length.out = 200)
plot(t, distan2[,2], type = 'l')
lines(t, distan2[,1], type = 'l')

distan2 = distances (b[-1],c[-1],a[-1],P1[,2],P1[, 3],P1[, 1])
t <- seq(0,1, length.out = 200)
plot(t, distan2[,2], type = 'l')
lines(t, distan2[,1], type = 'l')

distan2 = distances (c,b,a,P1[,3],P1[, 2],P1[, 1])
#t <- seq(0,1, length.out = 200)
plot(distan2[,2], type = 'l', col = 'red')
lines(distan2[,1], type = 'l', col = 'blue')
all(distan2[,1] >= distan2[,2])

distan2 = distances (a,c,b,P1[,1],P1[, 3],P1[, 2])
t <- seq(0,1, length.out = 200)
plot( distan2[,2], type = 'l')
lines( distan2[,1], type = 'l')

plot(distan2[,1], type = 'l')
lines(distan[,2], type = 'l', col = 'red' )

c = c(8,7)
b = c(8,2)
trLine = c()
t = seq(0,1, length.out = 200)
for (ti in t){
  trLine = rbind(trLine, tropicalLine(c,b,ti))
}
plot(trLine[,1], trLine[,2], type = 'l', col = 'red')




set.seed(12346) # one type three classified wrong
set.seed(123456) # 2 0 curved
set.seed(1234567) 
curvatures = c()
Ps= list()
for (i in 1:10){
  #a = runif(2, min = 0, max = 10)
  #b = runif(2, min = 0, max = 10)
  #c = runif(2, min = 0, max = 10)
  x = sample(1:10, size = 2, replace = TRUE)
  y  = sample(c(1:10)[-x], size = 2, replace = TRUE)
  z  = sample(c(1:10)[c(-x, -y)], size = 2, replace = TRUE)
  P = matrix(c(x,y, z),nrow=length(x))
  cv = curvature(P)
  curvatures = c(curvatures, cv)
    print(P)
    plotTR2(P, cv)
  Ps[[i]] = P
}
plot(1,1)
sum(curvatures==0)/length(curvatures) #0.645 # 0.6653
sum(curvatures==1)/length(curvatures) #0.282 # 0.246
sum(curvatures==-1)/length(curvatures) #0.073 #0.0887
sum(curvatures=="undefined")/length(curvatures) 
# same proportions in [-1, 0] and [0,1]

#0.71
# 0.21
# 0.08

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