
a = c(0,0)
b = c(2,4)
c = c(5,1)


a = c(0,4)
b = c(3,0)
c = c(5, 6)

# undefined curvature triangle

a = c(0,0)
b = c(448,449)
c = c(452,256)


P<- matrix(c(a,b, c),nrow=length(a))
K = 100000
N = 10
epsilon = 0.0001

#does not converge
K = 20000
N = 100
epsilon = 0.0001


dtr <- function (x, y){
  x = c(0,x)
  y = c(0,y)
  return (max(x-y) - min(x-y))
}
deu <- function (x, y) return (sqrt(sum((x-y)^2)))
sosd <- function (d, mu, P) {
  sos = 0
  for (i in 1:ncol(P)){
    x = P[,i]
    sos = sos + d(x,mu)^2
  }
  return (sos)
}

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


FrechetStrum <- function (P, K, N, epsilon){
  mu = P[,sample(ncol(P),size=1)]
  d = c()
  for (k in 1:K){
    p = P[,sample(ncol(P),size=1)]
    muk=mu
    tk = 1/(k+1)
    mu = tropicalLine(muk, p, tk )
    d = c(d, deu(muk, mu))
    ld = length(d)
    if (ld>N){
      if (all(d[c((ld-N):ld )]<epsilon) ) {
        return (c(mu, sosd(dtr, mu, P)))
      }
    }
  }
  crit = 0
  return (c("out for iterations:", mu))
}

FrechetStrum(P, K, N, epsilon)

frMeans <- matrix(NA, nrow = 50, ncol = 3)
for (i in 1:50){
  frMeans[i,] <- FrechetStrum(P, K, N, epsilon)
}

plot(c(0, 2, 5),c(0, 4, 1), pch = 16)
lines (c(0,2), c(0,2))
lines(c(1,5), c(1,1))
lines (c(2,2), c(1,4))
points(frMeans[,1], frMeans[,2], col='magenta', type='p',pch=16, ylim=c(-1,4), xlim= c(-1,5))
points(2,1, col='green', pch = 16)

plot(c(0, 3, 5),c(4, 0, 6), pch = 16)
lines(c(0,0), c(4,0))
lines(c(0,3), c(0,0))
lines(c(3,5), c(0,2))
lines(c(5,5), c(2,6))
lines(c(2,5), c(6,6))
lines(c(2,0) ,c(6,4))
points(frMeans[,1], frMeans[,2], col='magenta', type='p',pch=16, ylim=c(-1,4), xlim= c(-1,5))
points(3,4, col='green', pch = 16)

