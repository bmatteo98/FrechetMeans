a = c(0,0)
b = c(448,449)
c = c(452,256)
K = 20000
N = 100
epsilon = 0.001
P<- matrix(c(a,b,c),nrow=length(a))

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
        return (c(mu, sosd(dtr, mu, P)))
      }
    }
  }
  crit = 0
  #return (c(mu,crit))
  return (d)
}

distances <- FrechetBacak(P, K, N, epsilon)

x = 5000:20000
plot(log(x), log(distances[x]), type = 'p', pch = 46)
lines (log(x), -log(x)+4.2, col = 'red', ylim = c(-5.6, -4),type = 'p', pch = 46)
#lines (log(x), -0.98*log(x)+4.265, col = 'blue',  type = 'p', pch = 46)
lines (log(x), -0.98*log(x)+4.2, col = 'green',  type = 'p', pch = 46)
#lines (log(x), -0.99*log(x)+4.1, col = 'green',  type = 'l')
grid(20, 10)

upper <- c()
upperx <- c()
lower <- c()
lowerx <- c()
logx <- log(x)
logd <- log(distances)[5000:20000]
benchmark <- (-0.98)*log(x)+4.2
plot(logx, logd,type = 'p', pch = 46)
lines(logx, benchmark,type = 'p', pch = 46, col = 'red')
for (i in 1:length(x)){
  if (logd[i]>benchmark[i]){
    upperx <- c(upperx, i)
    upper <- c(upper, logd[i])
  }
  else {
    lowerx <- c(lowerx, i)
    lower <- c(lower, logd[i])
  }
}
upperx <- logx[upperx]
lowerx <- logx[lowerx]
lm1 <- lm(upper~upperx)
summary(lm1)
plot(upperx, upper,type = 'p', pch = 46, ylim=c(-5.5, -4), main = 'Interpolation of the log-distances for undefined curvature triangle')
lines(upperx, 4.517 -1.003*upperx,type = 'p', pch = 46, col='red')

lm2 <- lm(lower~lowerx)
summary(lm2)

benchmark2 <-  4.1-0.9945454 *lowerx
upper2 <- c()
upperx2 <- c()
lower2 <- c()
lowerx2 <- c()
for (i in 1:length(lowerx)){
  if (lower[i]>benchmark2[i]){
    upperx2 <- c(upperx2, i)
    upper2 <- c(upper2, lower[i])
  }
  else {
    lowerx2 <- c(lowerx2, i)
    lower2 <- c(lower2, lower[i])
  }
}
upperx2 <- lowerx[upperx2]
lowerx2 <- lowerx[lowerx2]
lm3 <- lm (lower2~lowerx2)
summary(lm3)
lines(lowerx2, lower2,type = 'p', pch = 46)
lines(lowerx2, 4.173 -1.002*lowerx2,type = 'p', pch = 46, col='blue')
lm4 <- lm (upper2~upperx2)
summary(lm4)
lines(upperx2, upper2,type = 'p', pch = 46)
lines(upperx2, 4.134-0.9952*upperx2,type = 'p', pch = 46, col='green')

library(ggplot2)

xL <- c(upperx, upperx2, lowerx2)
yL <- c(upper, upper2, lower2)
group <- c(rep(1, length(upperx)), rep(2, length(upperx2)), rep(3,length(lowerx2)) )
df <- as.data.frame(cbind (xL, yL, group))
ggplot(df, aes(x = xL, y = yL, color = as.factor(group)) )+
  geom_point( ) +
  geom_smooth(method = "lm")
