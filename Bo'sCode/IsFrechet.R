
#library(rje) # for powerset
#library(Dict)
library(lintools) # for is_feasable
library(limSolve)
library(lpSolve)
# distances & preparation

dtr <- function(l1, l2) max(l1-l2) - min(l1-l2) # tropical distance
ssq <- function(p,l) { # sum of squares of tropical distance
  ssq = 0
  for (x in l) {
    ssq = ssq + (dtr(x,p))^2
  }
  return(ssq)
}
NormalVec <-function (p,height) p-max(p)+height

# criterion for Frechet means

Peakpairs <- function (l1,l2){ # l1, l2 are vectors
  n <- length(l1)
  if (n != length(l2)) return ("unequal length of lists")
  d = dtr(l1, l2)
  o = c()
  for (i in 1:(n-1)){
    for (j in (i+1): n){
      if ((abs (l1[i]-l1[j]+l2[j]-l2[i] - d ) < 1e-5) | (abs (l1[i]-l1[j]+l2[j]-l2[i]) == d)){
        o = rbind(o, c(i,j))
      }
    }
  }
  return (o) 
}

L = cbind(c(1/5, 2/5, 2, 2/5, 2, 2), c(2, 2, 2, 2/5, 2/5, 1/5), c(2/5, 2/5, 2, 1/5, 2, 2))
p = c(22/15, 22/15, 2, 2/5, 2, 2)

Peakpairs(p, sp2[,3])

Peakpairs <- function (l1,l2){ # l1, l2 are vectors
  n <- length(l1)
  if (n != length(l2)) return ("unequal length of lists")
  dif = l1-l2
  maxim = max(dif)
  minim = min(dif)
  i = c()
  j = c()
  for (k in 1:length(dif)){
    if (dif[k] == maxim) i = c(i,k)
    if (dif[k] == minim) j = c(j,k)
  }
  return (as.matrix(expand.grid(i,j)))
}
Peakpairs(p, L[,1])

#Check the feasibility of a linear system
# follows from the theorem on weighted combinations of quadratic functions
#this
Peakpairs <- function (l1,l2){ # l1, l2 are vectors
  n <- length(l1)
  if (n != length(l2)) return ("unequal length of lists")
  diff = round(l1-l2, 3)
  maxim = max(diff)
  minim = min(diff)
  print(l1-l2)
  print(minim)
  print(maxim)
  i = which((diff) == maxim)
  j = which((diff) == minim)
  print(j)
  return (as.matrix(expand.grid(i,j)))
}
Peakpairs(p, L[,1])


IsFrechet <- function (p, L) {# L is a matrix, p is a vector
  m = ncol(L)
  n = length(p)
  allpeaks = Peakpairs(p, L[,1])
  allpeaks = cbind(rep(1, nrow(allpeaks)), allpeaks)
  for (i in 2:m){
    pp = Peakpairs(p, L[,i])
    pp = cbind(rep(i, nrow(pp)), pp)
    pp = pp[order(pp[,2]),]
    allpeaks = rbind(allpeaks,pp)
  }
  #allpeaks [4:6,] = rbind(c(2,1,4), c(2,1,5), c(2,1,6))
  S = matrix(0, nrow = n, ncol = nrow(allpeaks))
  for (j in 1: nrow(allpeaks)){
    peack = allpeaks[j,]
    S[peack[2], j] = p[peack[2]] - p[peack[3]] + L[ peack[3], peack[1]] - L[ peack[2], peack[1]]
    S[peack[3], j] = p[peack[3]] - p[peack[2]] + L[ peack[2], peack[1]] - L[peack[3], peack[1]]
    #print(S)
  }
  w = matrix(0, nrow = m, ncol = nrow(allpeaks))
  for (i in 1:m){
    ones = which(allpeaks[,1] == i)
    w[i,ones] = rep(1, length(ones))
  }
  pos = diag(rep(-1, nrow(allpeaks)), nrow =  nrow(allpeaks))
  b = c(rep(0,n), rep(1,m), rep(0,nrow(allpeaks) ))
  A = rbind(S,w)
  print(c(n,m,nrow(allpeaks)))
  return (A)
  #print(A)
  #showEqn(A, b)
  #return (is_feasible(A, b, neq = n+m, nleq = nrow(allpeaks) ))
}
nrows = 7
n = 6
m = 3
pos = diag(rep(1, nrows), nrow =  nrows)
A = IsFrechet(p,L)
b = c(rep(0,n), rep(1,m))
c = rep(1, nrows)
zeros = rep(0, n+m)

constr = rep(">=", n+m)
l = lp(direction="min", objective.in = c, const.mat =A, const.dir = constr,transpose.constraints = TRUE, const.rhs= zeros)
print(l$status)
print(l$solution)

library(Matrix)
rankMatrix(A)
rankMatrix(cbind(A, b))



sp2 = cbind(c(1/5, 2/5, 2, 2/5, 2, 2), c(2, 2, 2, 2/5, 2/5, 1/5), c(2/5, 2/5, 2, 1/5, 2, 2))
p = c(22/15, 22/15, 2, 2/5, 2, 2) #yes
p = c(2/5, 2/5, 2, -2/3, 14/15, 14/15) #yes
p = c(2/5, 2/5, 2, -2/3, 17/15, 14/15) #yes
p = c(2/5, 2/5, 2, 1/5, 17/15, 14/15) #yes
p = c(22/15, 22/15, 2, 2/5, 2, 2)# yes

p= c(19/15, 19/15, 9/5, 1/5, 2, 9/5)
p = c(22/15, 22/15, 2, 19/15, 2, 2)
p = c(19/15, 19/15, 9/5, 16/15, 2, 9/5)

p = c(9/5, 2, 2, 2, 2, 9/5) #no
p = c(1/3, 5/7, 0, 0, 1/2, 1)#no

l = IsFrechet(p,sp2)
showEqn(A, b)
all.equal(A, IsFrechet(p,sp2))
c( R(A), R(cbind(A,b)) ) 

p = c(2/5, 2/5, 2, -2/3, 14/15, 14/15)
IsFrechet(p, sp2)


provA <- matrix(rbinom(40,1,0.5), nrow =8, byrow = TRUE)
provA
provb <- rbinom(8,1,0.5)
is_feasible(A=provA,b=provb,neq=5,nleq=3)



sol = xranges(E = A, F = b, ispos=TRUE)
sol = ldei(E = A, F = b, G = pos, H = c)
