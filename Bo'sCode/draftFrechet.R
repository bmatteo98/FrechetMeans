
library(rje) # for powerset
library(Dict)
library(lintools) # for is_feasable
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
        o = cbind(o, c(i,j))
      }
    }
  }
  return (o) 
}

sp2 = cbind(c(1/5, 2/5, 2, 2/5, 2, 2), c(2, 2, 2, 2/5, 2/5, 1/5), c(2/5, 2/5, 2, 1/5, 2, 2))
p = c(22/15, 22/15, 2, 2/5, 2, 2)

Peakpairs(p, sp2[,3])
# my version 
Peakpairs <- function (l1,l2){ # l1, l2 are vectors
  n <- length(l1)
  if (n != length(l2)) return ("unequal length of lists")
  maxim = max(l1-l2)
  minim = min(l1-l2)
  i = which((l1-l2) == maxim)
  j = which((l1-l2) == minim)
  return (as.matrix(expand.grid(i,j)))
}
Peakpairs(p, sp2[,1])

#Check the feasibility of a linear system
# follows from the theorem on weighted combinations of quadratic functions
#this
Peakpairs <- function (l1,l2){ # l1, l2 are vectors
  n <- length(l1)
  if (n != length(l2)) return ("unequal length of lists")
  maxim = max(l1-l2)
  minim = min(l1-l2)
  i = which((l1-l2) == maxim)
  j = which((l1-l2) == minim)
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
    allpeaks = rbind(allpeaks,pp)
    }
  S = matrix(0, nrow = n, ncol = nrow(allpeaks))
  for (j in 1: nrow(allpeaks)){
    peack = allpeaks[j,]
    S[peack[2], j] = p[peack[2]] - p[peack[3]] + L[peack[1], k] - L[peack[1], peack[2]]
    S[peack[3], j] = p[peack[3]] - p[peack[2]] + L[peack[1], peack[2]] - L[peack[1], k]
  }
}

IsFrechet(p,sp2)









IsFrechet <- function (p, L) {# L is a matrix, p is a vector
  m = ncol(L)
  n = length(p)
  allpeaks = Peakpairs(p, L[,1])
  allpeaks = cbind(rep(1, nrow(allpeaks)), allpeaks)
  for (i in 2:m){
    pp = Peakpairs(p, L[,i])
    pp = cbind(rep(i, nrow(pp)), pp)
    allpeaks = rbind(allpeaks,pp)
  }
  S = matrix(0, nrow = n, ncol = nrow(allpeaks))
  for (j in 1: nrow(allpeaks)){
    peack = allpeaks[j,]
    S[peack[2], j] = p[peack[2]] - p[peack[3]] + L[ peack[3], peack[1]] - L[ peack[2], peack[1]]
    S[peack[3], j] = p[peack[3]] - p[peack[2]] + L[ peack[2], peack[1]] - L[peack[3], peack[1]]
  }
  w = matrix(0, nrow = m, ncol = nrow(allpeaks))
  for (i in 1:m){
    ones = which(allpeaks[,1] == i)
    w[i,ones] = rep(1, length(ones))
  }
  b = c(rep(0,n), rep(1,m), rep(0, nrow(allpeaks)))
  print(c(n,m,nrow(allpeaks)))
  pos = diag(rep(-1, nrow(allpeaks)), nrow =  nrow(allpeaks))
  A = as.matrix(rbind(S,w , pos))
  return (A)
}

A = IsFrechet(p,sp2)
A = matrix(c(1.8 , 0.0 , 0.0, -1.8, -1.8, -1.8,    0, 0.0 , 1.8 , 0.0,  0.0 , 0.0 , 0.0 ,   0,0.0 , 0.0 ,
             0.0 , 0.0,  0.0 , 0.0 ,   0 ,  0.0 , 0.0 , 1.8,  1.8 , 0.0 , 0.0 ,   2,  0.0 , 0.0 , 0.0 , 0.0,  1.8 ,
             0.0  ,  0, -1.8 ,-1.8, -1.8 , 0.0 , 0.0 , 1.8 ,  -2,1.0 , 1.0 , 1.0,  0.0 , 0.0 , 0.0 ,   0 ,
             0.0 , 0.0 , 0.0 , 1.0 , 1.0,  1.0 ,   0, 0.0 , 0.0 , 0.0 , 0.0,  0.0  ,0.0  ,  1), byrow = TRUE, nrow = 9)
pos = diag(rep(-1, 7), nrow =  7)
b = c(0 ,0 ,0 ,0, 0, 0,1 ,1 ,1 ,0 ,0, 0, 0, 0 ,0 ,0 ,0)
A <- rbind(A , pos)
lpSolve(A,c(rep(0,6), rep(1,3), rep(0, 12)))
is_feasible(A,b, neq=9,nleq=7)
p = c(9/5, 2, 2, 2, 2, 9/5)
sp2 = cbind(c(1/5, 2/5, 2, 2/5, 2, 2), c(2, 2, 2, 2/5, 2/5, 1/5), c(2/5, 2/5, 2, 1/5, 2, 2))
p = c(22/15, 22/15, 2, 2/5, 2, 2)

IsFrechet(p, sp2)
A = IsFrechet(p,sp2)
b = c(rep(0,6), rep(1,3))
showEqn(A, b)

c( R(A), R(cbind(A,b)) ) 














IsFrechet <- function (p, L) { # L is a matrix, p is a vector
  m = ncol(L)
  n = length(p)
  S = c()
  gsum = c(rep(0, n))
  
  peaks = matrix(0, ncol = m, nrow = 2)
  for (i in 1:m) peacks[,i] = Peackpairs(L[,i], p) # should return only 2 indices
  
  for (i in 1:m){
    si = 0
    for (pair in peaks[,i]) {
    j = pair[1]
    k = pair[2]
    S = c(S, c(w, i, j, k)) # >=0 ????
    si = si + c(w, i, j, k) # ???
    gsum[j] = gsum[j] + (p[j]-p[k]+ L[i,k] - L[i,j]) * c(w, i, j, k)
    gsum[k] = gsum[k] + (p[k]-p[j]+ L[i,j] - L[i,k]) * c(w, i, j, k)
    }
    #S = c(S, {si = 1}) ???
  }
  # S = c(S, {gsum[i] = 0})
  is_feasible (S)
}


#Step Descent approach to search for one exact Frechet mean
Infor <- function (p1, p2) { # information of peacks and valleys
  dif = p1 - p2
  maximum = max(dif)
  minimum = min(dif)
  peacks = c()
  valleys = c()
  
  for (index in 1:length(dif)){
    if (dif[index] == maximum) peaks = c(peacks, index)
    if (dif[index] == minimum) valleys = c(valleys, index)
  }
  return (c(peacks, valleys))
}

# difference between this and peack pairs? 
Infor <- function(l1,l2){ # l1, l2 are vectors
  n <- length(l1)
  if (n != length(l2)) return ("unequal length of lists")
  i = which(max(l1-l2))
  j = which(min(l1-l2))
  return (rbind(i,j)) # returns a matrix
}


#Informations of peacks and valleys with respect to all points in the sample
Infos <- function(pt, sample) { # sample is a matrix that contains points o the columns
  n = ncol(sample)
  point = sample[,1]
  infos <- Infor (pt, point)
  for (i in 2:n){
    point = sample[,i]
    infos <- cbind(infos, Infor (pt, point))
  }
  return(infos)
}

DescendPerturb <- function (p1, p2, incr, decr) {
  dif = p1 - p2
  ifr = Infor (p1, p2)
  peacks = ifr[,1]
  valleys = ifr[,2]
  if (peacks %in% incr) {ep = 1} #??
  else if (peacks %in% decr) {ep = -1 }
  else {ep = 0}
  
  if (valleys %in% incr) {ev = 1} #??
  else if (valleys %in% decr) {ev = -1 }
  else {ev = 0}
  slope = ep - ev
  maximum = max(dif)
  minimum = min(dif)
  bounds = c()
  
  for (index in 1: length(dif)){
    if ((ep = -1) & (length(intersect(index,incr))!=0)) {bounds = c(bounds, 0.5 * (maximum - dif[index]))}
    else if ((ep = 0) & (index %in% incr)) {bounds = c(bounds,  (maximum - dif[index]))}
    else if ((ep = -1) & (index %in% decr == FALSE)) {bounds = c(bounds, (maximum - dif[index]))}
    
    if ((ev = 1) & (length(intersect(index,decr))!=0)) {bounds = c(bounds, 0.5 * (dif[index] - minimum))}
    else if ((ev = 0) & (index %in% decr)) {bounds = c(bounds,  (dif[index] - minimum))}
    else if ((ev = 1) & (index %in% incr == FALSE)) {bounds = c(bounds, (dif[index] - minimum))}
  }
  if (length(bounds == 0)) dist = 0
  else dist = min(bounds)
  return (c(slope, dist, ep, ev))
}

MinQuad <- function (f, L){ # L is a matrix
  n = ncol(L)
  f <- expression (f)
  s <- c()
  for (i in 1:n) s <- c(s, D(f, L[,i])) # L[1,]=0 ??
  return (solve (s, L)) #[1] ??
}

EquivClass(L, sample, w){ # L contains on the columns
  n = ncol(L)
  if (n <= 1) return (cbind(L[,1],w))
  else {
    for (i in 1:n-1) {
      for (j in i+1: n){
        if (L[1, i] %in% L[1,j]) {return (EquivClass(L[,-i], sample, w))}
        else if (L[1, j] %in% L[1,i]) {return (EquivClass(L[,-j], sample, w))}
        else if (length(intersect(L[1, j],L[1,i]))!=0){
          new = L[, -c(i,j)] # [L[i][1] union L[j][1], 0]
          weigths = w
          cind = sort(intersect(convert(L[i][1],L[j][1])))[1] # max ??
          if (L[2,i] > 0){
            pt = sample(L[2, i]) #????
            
            for (k in L[1,i]){
              weigths[k] = weigths[cind] - pt[cind] + pt[k]
              
            }
          }
          return (EquivClass(new, sample, weights))
          
        }
      
      }
    }
    weights = w
    if (L[2,i] > 0){
      pt = sample(L[2, i]) #????
      
      for (k in L[1,i]){
        weigths[k] = weigths[cind] - pt[cind] + pt[k]
        
      }
    }
    return (cbind(L[,1], weights))
  }
}


Cand <- function (ifs, sample, height){
  m = ncol(sample)
  if (length(ifs) !=m) return ("error, lengths do not match")
  n = nrow(sample)
  infos = c()
  for (i in 1:m){
    infos = c(infos) #[ifs[i][1], i], [ifs[i][2], i] ??
  }
  ec = EquivClass(infos, sample, rep(0, n))
  atlas = ec[2]
  if (length(ec[,1]) ==1) return (NormalVec(atlas, height))
  for (i in 1:length(ec[,1])){
    for (ind in ec[1,i]){
      atlas [ind] = c(z, i) + atlas [ind] # z???
    }
  }
  obj = 0
  for (i in 1:m){
    indp = sort(ifs[1,i])[1] #index of peacks ??
    indv = sort(ifs[2,i])[1] 
    obj = obj + (atlas[indp]-atlas[indv]+sample[i, indv]-sample[i, indp])
  }
  sol = MinQuad(obj) #expand(obj), [seq(cat(z, i), i = 1 .. nops(ec[1]))]
  NormalVec(sol, height) #map(x -> subs(sol, x), atlas), height ??
}

OneFM <- function (pt, sample, height){
  if (IsFrechet(pt, sample)){return (NormalVec(pt, height))}
  m = ncol(sample)
  n = length(pt)
  status = FALSE
  current = pt
  c = 0
  while ((status != TRUE) & (c<=10) ){
    info = Infos (current, sample)
    dtrs = c()
    for (i in 1:m){
      x = sample[,i]
      dtrs = c(dtrs, dtr(current, x))
    }
    candidates = c()
    for (i in 1:m){
      valley = info[2,i]
      pw = powerSet(1:n)
      candidates = c(candidates, valley) #valley, pw
    }
    candidates = candidates [-(1:n)]
    opti = Dict$new[current= ssq(current, sample)]
    for incr in candidates{
      newpt = current

        
    }
  }
}

FlatPerturb

closure

closures






