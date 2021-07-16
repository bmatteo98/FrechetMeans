# skinny triangle, non-positively curved

a = c(0,0)
b = c(2,4)
c = c(5,1)
#1.459623   1.361469   1.000000 627.000000  19.636286
#  1.551845    1.363357    1.000000 4338.000000   19.352629

# fat triangle, positively curved 
a = c(0,4)
b = c(3,0)
c = c(5, 6)

# undefined curvature triangle

a = c(0,0)
b = c(448,449)
c = c(452,256)

P<- matrix(c(a,b,c),nrow=length(a))
#does not converge
K = 20000
N = 100
epsilon = 0.0001

#converges
K = 1000000
N = 100
epsilon = 0.0001

# tropical distance

dtr <- function (x, y){
  x = c(0,x)
  y = c(0,y)
  return (max(x-y) - min(x-y))
}
de <- function (x, y) return (sqrt(sum((x-y)^2)))
sosd <- function (d, mu, P) {
  sos = 0
  for (i in 1:ncol(P)){
    x = P[,i]
    sos = sos + d(x,mu)^2
  }
  return (sos)
}


w <- function (k) return (1/3)
#lambda <- function(k) return (1/(k*log(k)))
lambda <- function(k) return (1/(2*k))

bacak_sequence <- function (x,y, lambda, w) {
  if (identical(x,y)) {
    return (x)
    break
  }
  tk = (2*lambda*w)/(1+2*lambda*w)
  if ((x[1] == y[1]) | (x[2] == y[2])) {
    return ((1-tk)*x+y*tk)
    break
  }
  
  if ((x[1] == x[2]) & (y[1] == y[2])){
    return ((1-tk)*x+y*tk)
    break
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



FrechetBacak<- function (P, K, N, epsilon){
  mu = P[,sample(ncol(P),size=1)]
  d = c()
  for (k in 2:K){
    p = P[,sample(ncol(P),size=1)]
    muk=mu
    lamb = lambda(k)
    ww = w(k)
    mu = bacak_sequence(muk,  p,  lamb,ww)
    d = c(d, de(muk, mu))
    ld = length(d)
    if (ld>N){
      if (all(d[c((ld-N):ld )]<epsilon) ) {
        #return (c(mu,1, k, sosd(dtr, mu, P)))
        return (c(mu, sosd(dtr, mu, P)))
        }
    }
  }
  crit = 0
  return (c(mu, crit))
  #return (c(muk,crit))
  }

FrechetBacak (P, K, N, epsilon)

frMeans <- matrix(NA, nrow = 1000, ncol = 3)
for (i in 1:1000){
  frMeans[i,] <- FrechetBacak(P, K, N, epsilon)
}
lines(frMeans[,1], frMeans[,2], col='magenta', type='p',pch=46, ylim=c(-1,6), xlim= c(-1,5))
grid(10, 8)


