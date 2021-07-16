

trSegment = trSegmentN(c,a)
lines(trSegment[,2], trSegment[,3], type = 'l', col = 'red')
P = matrix(c(a,b, c),nrow=length(a))

plotTRN <- function (P){
  a = P[,1]
  b = P[,2]
  c = P[,3]
  xl = c(min(a[2], b[2], c[2]),max(a[2], b[2], c[2]) )
  yl = c(min(a[3], b[3], c[3]),max(a[3], b[3], c[3]) )
  trSegment = trSegmentN(a,b)
  par(mfrow=c(2,2))
  plot(trSegment[,2], trSegment[,3], type = 'l', col = 'red', xlim = xl, ylim = yl)
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
plotTRN(P)

a = c(5,1)
c = c(9,3)
b = c(10,6)
trLine = c()
t = seq(0,1, length.out = 200)
for (ti in t){
  trLine = rbind(trLine, tropicalLine(c,b,ti))
}
plot(trLine[,1], trLine[,2], type = 'l', col = 'red')


P = matrix(c(a,b, c),nrow=length(a))
plotTR2(P)
plotTR2 <- function (P){
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
  plot(trLine[,1], trLine[,2], type = 'l', col = 'red', xlim = xl, ylim = yl)
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







plotDists <- function (P){
  a = P[,1]
  b = P[,2]
  c = P[,3]
  P1 = findEuclideanN(P)
  distan = distancesN (a,b,c,P1[,1],P1[,2],P1[, 3])
  plot(distan[,2], type = 'l', col = 'blue')
  lines( distan[,1], type = 'l', col = 'red')
}
