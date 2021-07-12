# skinny triangle, non-positively curved

a = c(0,0)
b = c(2,4)
c = c(5,1)
#  1.575058   1.225900   1.000000 164.000000  18.969101
# 1.681355   1.354987   1.000000 539.000000  19.014077

# fat triangle, positively curved 
a = c(0,4)
b = c(3,0)
c = c(5, 6)

# undefined curvature triangle

a = c(0,0)
b = c(448,449)
c = c(452,256)

P<- matrix(c(a,b, c),nrow=length(a))
K = 20000
N = 100
epsilon = 0.001

# tropical distance

dtr <- function (x, y) return (max(x-y) - min(x-y))
de <- function (x, y) return (sqrt(sum((x-y)^2)))
sosd <- function (d, mu, P) {
  sos = 0
  mu = c(0,mu)
  for (i in 1:ncol(P)){
    x = c(0,P[,i])
    sos = sos + d(x,mu)^2
  }
  return (sos)
}

# compute inductive mean on the tropical line segment 

#set.seed(1)

inductive_mean <- function (x,y,k) {
  if (identical(x,y)) {
    return (x)
    break
  }
  if ((x[1] == y[1]) | (x[2] == y[2])) {
    return ((1-1/k)*x+y/k)
    break
  }
  
  if ((x[1] == x[2]) & (y[1] == y[2])){
    return ((1-1/k)*x+y/k)
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
  #print (c(a,b))
  if ((a[2]<b[2]) & ((a[1]-a[2])<=(b[1]-b[2]))){
    p = c(a[1]+b[2]-a[2], b[2])
    t = (b[2]-a[2])/(b[1]-a[1])
    if ((1/k) <= t) {
      #mu = c(a[1]+(1/k), a[2]+(1/k))*(b[1]-a[1])
      mu = (1-(1/k))*x+p/k
    }
    if (t < (1/k)) {
      #mu = c(a[1]+(1/k), b[2])*(b[1]-a[1])
      dk = 1/k - t
      mu = (1-dk)*p+y*dk
    }
  }
  
  if ((a[2]<b[2]) & ((a[1]-a[2])>(b[1]-b[2]))){
    p = c(b[1], a[2]+b[1]-a[1])
    t = (b[1]-a[1])/(b[2]-a[2])
    if (t >= (1/k)) {
      #mu = c(a[1]+(1/k), a[2]+(1/k))*(b[2]-a[2])
      mu = (1-(1/k))*x+p/k
      #print(c(p,x))
    }
    if (t < (1/k)) {
      #mu = c(b[1],a[2]+(1/k))*(b[2]-a[2])
      dk = 1/k - t
      mu = (1-dk)*p+y*dk
    }
  }
  
  if (a[2]>b[2]){
    p = c(a[1], b[2])
    t = (a[2]-b[2])/(a[2]-b[2]+b[1]-a[1])
    if (t >= (1/k)) {
      #mu = c(a[1],a[2]-(1/k))*(a[2]-b[2]+b[1]-a[1])
      mu = (1-(1/k))*x+p/k
    }
    if (t < (1/k)) {
      #mu = c(a[1]+b[2]-a[2]+(1/k),b[2])*(a[2]-b[2]+b[1]-a[1])
      dk = 1/k - t
      mu = (1-dk)*p+y*dk
    }
  }
  #print(p)
  return (mu)
}
#inductive_mean(a,c,2)
#points(c(0,5, 1,2.2), c(0,1, 1, 1 ), pch = 16, col = c('blue', 'blue', 'red', 'green'))
#dtr(inductive_mean(a,c,2), a)
#dtr(inductive_mean(a,c,2), c)
#dtr(a,c)
#Strum's algorithm

FrechetStrum <- function (P, K, N, epsilon){
  mu = P[,sample(ncol(P),size=1)]
  d = c()
  for (k in 1:K){
    p = P[,sample(ncol(P),size=1)]
    #print(p)
    muk=mu
    mu = inductive_mean(muk, p, k+1)
    d = c(d, de(muk, mu))
    #print(sosd(dtr, mu, P))
    ld = length(d)
    #print(mu)
    if (ld>N){
      if (all(d[c((ld-N):ld )]<epsilon) ) {
        #return (c(mu,1, k, sosd(dtr, mu, P)))
        return (c(mu, sosd(dtr, mu, P)))
        }
    }
  }
  crit = 0
  return (c(mu,crit))
  #return (d)
}

FrechetStrum(P, K, N, epsilon)

frMeans <- matrix(NA, nrow = 50, ncol = 3)
for (i in 1:50){
  frMeans[i,] <- FrechetStrum(P, K, N, epsilon)
}
plot(s[mins[,1]], s[mins[,2]], type='p',pch=16, ylim=c(-1,6), xlim= c(-1,5), col = 'green')
points(c(0, 2, 5),c(0, 4, 1), pch = 16)
lines (c(0,2), c(0,2))
lines(c(1,5), c(1,1))
lines (c(2,2), c(1,4))
plot(frMeans[,1], frMeans[,2], col='magenta', type='p',pch=16, ylim=c(-1,4), xlim= c(-1,5))

#grid(6, 7)

