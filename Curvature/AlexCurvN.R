
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
  #print(L)
  t = seq(0,1,length.out = 100)
  TrSegment = matrix(NA, nrow = 1, ncol = N)
  for (j in 1:(N-1)){
    y1 = L[j,]
    y2 = L[j+1,]
    segment = matrix(NA, nrow = length(t), ncol = N)
    for (i in 1:length(t)){
      tk = t[i]
      segment[i,] =  (1-tk)*y1+tk*y2
    }
    TrSegment = rbind(TrSegment, segment)
  }
  return (TrSegment[-1,])
}

trSegment = trSegmentN(c,b)
plot(trSegment[,1], trSegment[,2], type = 'l')


euclSegmentN <- function(u,v){
  N = length(u)
  t = seq(0,1, length.out = 100*(N-1))
  EuclSegment = matrix(NA, nrow = length(t), ncol = N)
  for(i in 1:length(t)){
    tk = t[i]
    EuclSegment[i,] =  (1-tk)*v+tk*u
  }
  return (EuclSegment)
}

euclSegment = euclSegmentN(c,b)
plot(euclSegment[,1], euclSegment[,2], type = 'l')

dtr <- function (x, y) return (max(x-y) - min(x-y))
deu <- function (x, y) return (sqrt(sum((x-y)^2)))

distancesN <- function(a, b, c, a1, b1, c1){
  euclSeg = euclSegmentN(a1,b1)
  trSeg = trSegmentN(a,b)
  n = nrow(euclSeg)
  distanN = matrix(NA, nrow = n, ncol = 2)
  for (i in 1: n){
    ptTR = trSeg[i,]
    ptEU = euclSeg[i,]
    distanN[i,] = c(dtr(ptTR, c), deu(ptEU, c1))
  }
  return (distanN)
}
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
  b1 = c(dab, rep(0, N-1))
  xc = (dac^2-dbc^2+dab^2)/(2*dab)
  if (abs(round(xc, 5)) == abs(round(dac, 5))) yc = 0
  else yc = sqrt(dac^2-xc^2)
  c1 = c(xc, yc, rep(0, N-2))
  P1 <- matrix(c(a1,b1, c1),nrow=length(a1))
  return (P1)
}

#P1 = findEuclidean(P)

curvature <- function (P){
  P1 <- findEuclideanN(P)
  curvS = 0
  curvF = 0
  for (j in 1:ncol(P)){
    a = P[,j%%ncol(P)+1]
    b = P[,(j+1)%%ncol(P)+1]
    c = P[,(j+2)%%ncol(P)+1]
    a1 = P1[,j%%ncol(P1)+1]
    b1 = P1[,(j+1)%%ncol(P1)+1]
    c1 = P1[,(j+2)%%ncol(P1)+1]
    distan = round(distancesN (a,b,c,a1,b1,c1), 5)
    if (all(distan[,1] <= distan[,2])) curvS = curvS +1 
    else if (all(distan[,1] >= distan[,2])) curvF = curvF +1
    else (return (0))
  }
  if (curvS == 3) return (-1)
  if (curvF == 3) return (1)
  else return (0)
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


distan = distancesN (a,b,c,P1[,1],P1[, 2],P1[, 3])
#t = seq(0,1, length.out = nrow(distan))
plot( distan[,1], type = 'l')
lines(distan[,2], type = 'l')

distan = distancesN (a,c,b,a1,c1,b1)
plot( distan[,1], type = 'l')
lines(distan[,2], type = 'l')

distan = distancesN (c,b,a,c1,b1,a1)
plot( distan[,1], type = 'l')
lines(distan[,2], type = 'l')


