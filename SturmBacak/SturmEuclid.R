
a = c(0,0)
b = c(2,4)
c = c(5,1)

P<- matrix(c(a,b,c),nrow=length(a))
K = 10000
N = 10
epsilon = 0.0001

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


inductive_mean <- function (x,y,k) {
  mu = (1-1/k)*x+y/k
  return (mu)
}

FrechetStrum <- function (P, K, N, epsilon){
  mu = P[,sample(ncol(P),size=1)]
  d = c()
  for (k in 1:K){
    p = P[,sample(ncol(P),size=1)]
    muk=mu
    mu = inductive_mean(muk, p, k+1)
    d = c(d, de(muk, mu))
    ld = length(d)
    if (ld>N){
      if (all(d[c((ld-N):ld )]<epsilon) ) {
        #return (c(mu,1, k, sosd(dtr, mu, P)))
        return (c(mu, sosd(de, mu, P)))
      }
    }
  }
  crit = 0
  return (c(mu,crit))
  #return (d)
}

FrechetStrum(P, K, N, epsilon)
plot(c(0, 2, 1),c(0, 4, 2), pch = 16, col = c('black', 'black', 'red'))

plot(c(0, 2, 5, 2.311688),c(0, 4, 1, 1.674326), pch = 16, col = c('black', 'black', 'black', 'red'))
