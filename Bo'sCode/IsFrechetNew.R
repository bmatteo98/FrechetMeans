library(limSolve)

Peakpairs <- function (l1,l2){ # l1, l2 are vectors
  n <- length(l1)
  if (n != length(l2)) return ("unequal length of lists")
  diff = round(l1-l2, 3)
  maxim = max(diff)
  minim = min(diff)
  i = which((diff) == maxim)
  j = which((diff) == minim)
  return (as.matrix(expand.grid(i,j)))
}


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
  b = c(rep(0,n), rep(1,m))
  A = rbind(S,w)
  sol = xranges(E = A, F = b, ispos=TRUE)
  return (is.numeric(sol))
  
  }


sp2 = cbind(c(1/5, 2/5, 2, 2/5, 2, 2), c(2, 2, 2, 2/5, 2/5, 1/5), c(2/5, 2/5, 2, 1/5, 2, 2))
p = c(22/15, 22/15, 2, 2/5, 2, 2) #yes
p = c(2/5, 2/5, 2, -2/3, 14/15, 14/15) #yes
p = c(2/5, 2/5, 2, -2/3, 17/15, 14/15) #yes
p = c(2/5, 2/5, 2, 1/5, 17/15, 14/15) #yes
p = c(22/15, 22/15, 2, 2/5, 2, 2)# yes
p= c(19/15, 19/15, 9/5, 1/5, 2, 9/5)# yes
p = c(22/15, 22/15, 2, 19/15, 2, 2)# yes
p = c(19/15, 19/15, 9/5, 16/15, 2, 9/5)# yes

p = c(9/5, 2, 2, 2, 2, 9/5) #no
p = c(1/3, 5/7, 0, 0, 1/2, 1)#no
p = c(1,2,3,1,2,3)#no

IsFrechet(p,sp2)


sp1 = cbind(c(9/100, 19/50, 19/50, 1, 19/50, 19/50, 1, 19/50, 1, 1), c(1, 1, 1, 1, 21/100, 57/100, 61/100, 57/100, 61/100, 61/100), c(31/50, 1, 49/50, 49/50, 1, 49/50, 49/50, 1, 1, 63/100), c(1, 1, 1, 47/100, 7/50, 4/5, 1, 4/5, 1, 1))
p = c(48/25, 48/25, 48/25, 193/100, 8/5, 43/25, 2, 43/25, 2, 193/100) #yes
p = c(171/100, 2, 2, 43/25, 139/100, 151/100, 43/25, 7/4, 43/25, 43/25) #yes
p = c(171/100, 2, 171/100, 43/25, 139/100, 151/100, 179/100, 151/100, 43/25, 43/25) #yes
p = c(171/100, 171/100, 2, 43/25, 139/100, 151/100, 43/25, 151/100, 43/25, 43/25) #yes
p = c(171/100, 2, 171/100, 43/25, 139/100, 151/100, 179/100, 7/4, 179/100, 43/25) #yes
p = c(49/25, 49/25, 49/25, 197/100, 41/25, 44/25, 197/100, 2, 197/100, 197/100) #yes

p = c(0, 2, 2, 1, 139/100, 3, 43/25, 7/4, 1/2, 43/25)
p = c(171/100, 2, 171/100, 43/25, 139/100, 1, 2, 2,  3 , 1/13)
p = c(171/100, 2, 171/100, 43/25, 2, 1, 2,171/100, 2, 171/100)
IsFrechet(p,sp1)
