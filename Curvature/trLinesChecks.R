
tropicalLine <-  function (x,y, tk) {
  if (identical(x,y)) {
    mu =  (x)
    
  }
  if ((x[1] == y[1]) | (x[2] == y[2])) {
    mu = ((1-tk)*x+y*tk)
    
  }
  
  if ((x[1] == x[2]) | (y[1] == y[2])){
    mu = ((1-tk)*x+y*tk)
    
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

a = runif(2, min = -1, max = 0)
b = runif(2, min = -1, max = 0)
c = runif(2, min = -1, max = 0)

trLine = c()
t = seq(0,1, length.out = 200)
for (ti in t){
  trLine = cbind(trLine, tropicalLine(c[-1],b[-1],ti))
}
plot(trLine[1,], trLine[2,], type = 'l')




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
plot(trSegment[,2], trSegment[,3], type = 'l', col = 'red')
points(a[2], a[3], phc = 16)


a = c(0.0000000 ,-0.4224496, -0.7683783)
b = c(0.0000000, -0.9718798, -0.1232577)
c = c(0.0000000 ,-0.5606527, -0.7082600)
dtr(c(-0.9718798,-0.5272726), c(-0.5606527, -0.7082600))
dtr(c(0,-0.9718798,-0.5272726), c(0,-0.5606527, -0.7082600))


distN = c()
for (i in 1:101){
  distN = c(distN, dtr(a, trSegment[i,]))
}
plot(distN, type = 'l')
