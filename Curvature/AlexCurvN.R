dtr <- function (x, y) return (max(x-y) - min(x-y))
deu <- function (x, y) return (sqrt(sum((x-y)^2)))

trSegmentN <- function (u, v){
  N = length(u) 
  lambda = v-u
  #print(lambda)
  sortLambda = sort(lambda)
  L = matrix(NA, nrow = N, ncol = N)
  y = v 
  L[1,] = y
  for (i in 2:(N-1)){
    li = sortLambda[i]
    y = pmin(u + li,v ) 
    y = y - y[1]
    L [i,] = y
  }
  L[N,] = u
  a = 0
  TrSegment = c()
  for (j in 1:(N-1)){
    y1 = L[j,]
    y2 = L[j+1,]
    b = a+dtr(y1, y2)
    t = seq(a,b,by=0.01)
    if (a == b) segment = y1
    else{
      segment = matrix(NA, nrow = length(t), ncol = N)
      for (i in 1:length(t)){
        tk = t[i]
        segment[i,] =  (1-(tk-a)/(b-a))*y1+(tk-a)*y2/(b-a)
      }
    }
    TrSegment = rbind(TrSegment, segment)
    a = b
  }
  
  return (TrSegment)
}

trSegment = trSegmentN(b,c)
#plot(trSegment[,2], trSegment[,3], type = 'l', col = 'red')
#trSegment = trSegmentN(c,b)
#plot(trSegment[,1], trSegment[,2], type = 'l')


euclSegmentN <- function(u,v,d){
  N = length(u)
  a = 0
  b = dtr(u,v)
  t = seq(a,b,length.out = d)
  EuclSegment = matrix(NA, nrow = length(t), ncol = N)
  for(i in 1:length(t)){
    tk = t[i]
    EuclSegment[i,] =  (1-(tk-a)/(b-a))*v+(tk-a)*u/(b-a)
  }
  return (EuclSegment)
}

euclSegment = euclSegmentN(c1,b1,998)
#plot(euclSegment[,1], euclSegment[,2], type = 'l')

distancesN <- function(a, b, c, a1, b1, c1){
  trSeg = trSegmentN(a,b)
  euclSeg = euclSegmentN(a1,b1, nrow(trSeg))
  #euclSeg = cbind(rep(0,998), EUsegment)
  #trSeg = cbind(rep(0,998), trLine)
  n = nrow(euclSeg) # they should be the same
  distanN = matrix(NA, nrow = n, ncol = 2)
  for (i in 1: n){
    ptTR = trSeg[i,]
    ptEU = euclSeg[i,]
    distanN[i,] = c(dtr(ptTR, c), deu(ptEU, c1))
  }
  return (distanN)
}
dist = distancesN(b,c,a,b1,c1,a1)
#dist = distancesN(P[,1],P[,2], P[,3], P1[,1], P1[,2],P1[,3] )

findEuclideanN <- function (P){
  a = P[,1]
  b = P[,2]
  c = P[,3]
  N = length(a)
  dab = dtr(a,b)
  dbc = dtr(b,c)
  dac = dtr(c,a)
  a1 = rep(0,N)
  b1 = c(0, dab, rep(0, N-2))
  xc = (dac^2-dbc^2+dab^2)/(2*dab)
  if (abs(round(xc, 8)) == abs(round(dac, 8))) yc = 0
  else yc = sqrt(dac^2-xc^2)
  c1 = c(0, xc, yc, rep(0, N-3))
  P1 <- matrix(c(a1,b1, c1),nrow=length(a1))
  return (P1)
}

#P1 = findEuclidean(P)

curvature <- function (P){
  P1 <- findEuclideanN(P)
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
    distan = round(distancesN (a,b,c,a1,b1,c1), 8)
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
plotTRN <- function (P,cv){
  a = P[,1]
  b = P[,2]
  c = P[,3]
  xl = c(min(a[2], b[2], c[2]),max(a[2], b[2], c[2]) )
  yl = c(min(a[3], b[3], c[3]),max(a[3], b[3], c[3]) )
  trSegment = trSegmentN(a,b)
  par(mfrow=c(2,2))
  plot(trSegment[,2], trSegment[,3], type = 'l', col = 'red', xlim = xl, ylim = yl, main = cv)
  trSegment = trSegmentN(c,b)
  lines(trSegment[,2], trSegment[,3], type = 'l', col = 'blue' )
  trSegment = trSegmentN(a,c)
  lines(trSegment[,2], trSegment[,3], type = 'l', col = 'green')
  points(b[2], b[3], pch = 16)
  points(c[2], c[3], pch = 16)
  points(a[2], a[3], pch = 16)
  P1 = findEuclideanN(P)
  distan = distancesN (a,b,c,P1[,1],P1[,2],P1[, 3])
  plot(distan[,2], type = 'l', col = 'blue', main = 'c from a-b')
  lines( distan[,1], type = 'l', col = 'red')
  distan = distancesN (a,c,b,P1[,1],P1[,3],P1[,2])
  plot(distan[,2], type = 'l', col = 'blue', main = 'b from a-c')
  lines( distan[,1], type = 'l', col = 'red')
  distan = distancesN (b,c,a,P1[,2],P1[,3],P1[,1])
  plot(distan[,2], type = 'l', col = 'blue', main = 'a from b-c')
  lines( distan[,1], type = 'l', col = 'red')
}


a = c(1,2,34,1,0)
b = c(0.2,4,1,0,2)
c = c(2,0,1,0.3,2.3)

# skinny triangle
a = c(0,1,3)
b= c(0,0,0)
c= c(0,3,2)


# fat triangle
a = c(0,0,2)
b = c(0,1,0)
c = c(0,3,3)

#undefined on the paper
# but it comes positive
a = c(0,0,2)
b = c(0,2,0)
c = c(0,3,4)


#fat
a = c(0,0,4)
b = c(0,3,0)
c = c(0,5, 6)


#skinny 
a = c(0,0,0)
b = c(0,2,4)
c = c(0,5,1)

#undefined
a = c(0,0,0)
b = c(0,448,449)
c = c(0,452,256)

# could there be 0 curved regions??
a = c(0,0)
b = c(1, 0)
c = c(1,1)

a = c(0,a)
b = c(0,b)
c = c(0,c)
P = matrix(c(a,b, c),nrow=length(a))
curvature(P)


P1 <- findEuclideanN(P)
a1 = P1[,1]
b1= P1[, 2]
c1 = P1[, 3]
dtr(P[,1], P[,2])
dtr(P[,1], P[,3])
dtr(P[,3], P[,2])

deu(P1[,1], P1[,2])
deu(P1[,1], P1[,3])
deu(P1[,3], P1[,2])


distan = distancesN (b,a,c,P1[,2],P1[, 1],P1[, 3])
#t = seq(0,1, length.out = nrow(distan))
plot(distan[,2], type = 'l', col = 'blue')
lines( distan[,1], type = 'l', col = 'red')

distan = distancesN (c,a,b,c1,a1,b1)
plot(distan[,2], type = 'l')
lines( distan[,1], type = 'l')

distan = distancesN (b,c,a,b1,c1,a1)

plot(distan[,2], type = 'l', col = 'red')
lines( distan[,1], type = 'l')


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















