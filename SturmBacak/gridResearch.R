a = c(0,0)
b = c(2,4)
c = c(5,1)

a = c(0,4)
b = c(3,0)
c = c(5, 6)

a = c(0,0)
b = c(448,449)
c = c(452,256)

a = c(0,0)
b = c(3,1)
c=c(2,5)
P<- matrix(c(a,b,c),nrow=length(a))

dtr <- function (x, y){
  x = c(0,x)
  y = c(0,y)
  return (max(x-y) - min(x-y))
}
sosd <- function (d, mu, P) {
  sos = 0
  for (i in 1:ncol(P)){
    x = P[,i]
    sos = sos + d(x,mu)^2
  }
  return (sos)
}


sosd <- function (d, mu, P) {
  sos = 0
  mu = c(0,mu)
  for (i in 1:ncol(P)){
    x = c(0,P[,i])
    sos = sos + d(x,mu)^2
  }
  return (sos)
}

sod <- function (d, mu, P) {
  sod = 0
  mu = c(0,mu)
  for (i in 1:ncol(P)){
    x = c(0,P[,i])
    sod = sod + d(x,mu)
  }
  return (sod)
}

sod(dtr, c(1,1), P)
sod(dtr, c(2,1), P)
sod(dtr, c(2,2), P)



s <- seq(0,450, by = 1)
ls = length(s)
ssq <- matrix(NA, ls, ls)

for (i in 1:ls){
  for (j in 1:ls){
    mu = c(s[i],s[j])
    ssq[i,j] = sod (dtr, mu, P)
  }
}
min(ssq)
mins <- which(ssq == min(ssq), arr.ind = T)
plot(s[mins[,1]], s[mins[,2]], type='p',pch=16, ylim=c(0,4), xlim= c(0,5), 
     col = 'orange', main='Tropical Fermat-Weber points', cex.main=2)
#points(c(0, 2, 5),c(0, 4, 1), pch = 16)

ssq <- matrix(NA, ls, ls)
for (i in 1:ls){
  for (j in 1:ls){
    mu = c(s[i],s[j])
    ssq[i,j] = sosd (dtr, mu, P)
  }
}
min(ssq)
mins <- which(ssq == min(ssq), arr.ind = T)
plot(s[mins[,1]], s[mins[,2]], type='p',pch=16, ylim=c(-3,7), xlim= c(-3,7), col = 'green')
plot(c(0, 2, 5),c(0, 4, 1), pch = 16)
lines (c(0,2), c(0,2))
lines(c(1,5), c(1,1))
lines (c(2,2), c(1,4))




plot(c(0, 2),c(0, 4),xlim=c(0,4), pch = 16)
lines (c(2,2), c(2,4))
lines(c(0,2), c(0,2))
lines(c(0,2),c(0,4))
points(1,2, col = 'red', pch = 16)

plot(c(, 3),c(4, 1), ylim = c(0,4), xlim = c(0,4),pch = 16)
lines(c(0,1), c(0, 1))
lines(c(1,2), c(1,1))

plot(c(0, 2),c(0, 1), ylim = c(0,2), pch = 16)
lines(c(0,1), c(0, 1))
lines(c(1,2), c(1,1))

plot(c(0, 3, 5),c(4, 0, 6), pch = 16)
lines(c(0,0), c(4,0))
lines(c(0,3), c(0,0))
lines(c(3,5), c(0,2))
lines(c(5,5), c(2,6))
lines(c(3,5), c(6,6))
lines(c(3,0) ,c(6,4))


plot(c(0, 448, 452), c(0, 449, 256), pch = 16)
lines(c(452, 448), c(256, 256))
lines(c(448, 448), c(256, 449))
lines(c(449, 448), c(449, 449))
lines(c(0, 449), c(0, 449))
lines(c(256, 452), c(256, 256))
points(s[mins[,1]], s[mins[,2]], type='p',pch=16,  col = 'green')


for (i in 1:nrow(mins)){
  fm <- c(s[mins[i,1]], s[mins[i,2]]) 
  print(c(dtr(P[,1], fm), dtr(P[,2], fm), dtr(P[,3], fm)))
}