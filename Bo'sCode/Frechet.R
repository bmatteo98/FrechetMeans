# 2020 r code for computing Frechet means under the tropical metric 
library(limSolve) # feasibility 
library(Deriv) # compute derivative of a poly
library(rSymPy)  # solve derivative = 0
library(rje) # powerset

# criterion of Frechet means

Peakpairs <- function (l1,l2){ 
  # pairs of indices that attain the tropical distrance
  # l1, l2 are vectors
  # l1 is a potential Frechet mean
  # l2 is a point in the sample
  n <- length(l1)
  if (n != length(l2)) return ("unequal length of lists")
  diff = round(l1-l2, 3)
  maxim = max(diff)
  minim = min(diff)
  i = which((diff) == maxim)
  j = which((diff) == minim)
  # all combinations of the indices that give the mins and the maxs
  return (as.matrix(expand.grid(i,j)))
}


IsFrechet <- function (p, L) {
  # checks feasibility of the linear system
  # follows from theorem on weighted combinations of quadratic functions
  # L is a (m x n) matrix, m = dimension of the points, n =  number of points in the sample
  # L contains a sample of points on the columns
  # p is a vector, potential Frechet mean
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
  b = c(rep(0,n), rep(1,m))
  A = rbind(S,w)
  sol = xranges(E = A, F = b, ispos=TRUE)
  return (is.numeric(sol))
  
}

# normalize to a certain height (same point in the tropical projective torus)
NormalVec <-function (p,height) p-max(p)+height

# Steepest descent approach to search for one exact FM

Infor <- function (l1,l2, pv){ 
  # Information on peacks and valleys
  # l1, l2 are vectors
  n <- length(l1)
  if (n != length(l2)) return ("unequal length of lists")
  dif = round(l1-l2, 4)
  maxim = max(dif)
  minim = min(dif)
  i = c()
  j = c()
  for (k in 1:length(dif)){
    if (dif[k] == maxim) i = c(i,k)
    if (dif[k] == minim) j = c(j,k)
  }
  if (pv == 1) return (i)
  if (pv == 2) return (j)
}



DescendPerturb <- function (p1, p2, incr, decr) {
  dif = p1 - p2
  peaks = Infor(p1, p2, 1)
  valleys = Infor(p1, p2, 2)
  
  if (length(intersect(peaks,incr)) != 0){ep = 1} 
  else if (all(peaks %in% decr)) {ep = -1 }
  else {ep = 0}
  if (length(intersect(valleys,decr)) != 0) {ev = -1} 
  else if (all(valleys %in% incr)) {ev = 1 }
  else {ev = 0}
  
  slope = ep - ev
  maximum = max(dif)
  minimum = min(dif)
  bounds = c()
  
  for (index in 1: length(dif)){
    if ((ep == -1) & (index %in% incr)) {bounds = c(bounds, 0.5 * (maximum - dif[index]))}
    else if ((ep == 0) & (index %in% incr)) {bounds = c(bounds,  (maximum - dif[index]))}
    else if ((ep == -1) & ((index %in% decr) == FALSE)) {bounds = c(bounds, (maximum - dif[index]))}
    
    if ((ev == 1) & (index %in% decr)) {bounds = c(bounds, 0.5 * (dif[index] - minimum))}
    else if ((ev == 0) & (index %in% decr)) {bounds = c(bounds,  (dif[index] - minimum))}
    else if ((ev == 1) & (index %in% incr == FALSE)) {bounds = c(bounds, (dif[index] - minimum))}
  }
  if (length(bounds) == 0) {dist = 0}
  else {dist = min(bounds)}
  
  return (c(slope, dist, ep, ev))
}

Infos <- function(p1, P){
  # information on Fpeaks and valleys with respect to all points in the sample
  L = list()
  j = 0
  for (i in 1:ncol(P)){
    p2 = P[,i]
    peaks = Infor(p1,p2,1)
    valleys = Infor(p1,p2,2)
    Lp = list(peaks, i)
    Lv = list(valleys, i)
    j = j+1
    L[[j]] = Lp
    j = j+1
    L[[j]] = Lv
  }
  return(L)
}

EquivClass<- function(L, P, w){ # L contains on the columns
  if (length(L) <= 1) {
    return (list(L,w))}
  else {
    for (i in 1:(length(L)-1)) {
      for (j in (i+1): length(L)){
        if (all(unlist(L[[i]][1]) %in% unlist(L[[j]][1]))) {
          L[[i]] = NULL
          return (EquivClass(L, P, w))
        }
        else if (all(unlist(L[[j]][1]) %in% unlist(L[[i]][1])) ){
          L[[j]] = NULL
          return (EquivClass(L, P, w))
        }
        
        else if (length(intersect(unlist(L[[j]][1]),unlist(L[[i]][1])))!=0){
          new = L
          new[[j]] = NULL
          new[[i]] = NULL
          new[[length(L)+1]] = list(union(unlist(L[[j]][1]), unlist(L[[i]][1])), 0)
          ww = w
          cind = sort(intersect(unlist(L[[i]][1]),unlist(L[[j]][1])))[1] 
          if (L[[i]][2] > 0){
            pt = P[,unlist(L[[i]][2])]
            for (k in unlist(L[[i]][1])){
              ww[k] = ww[cind] - pt[cind] + pt[k]
              
            }
          }
          if (L[[j]][2] > 0){
            pt = P[,unlist(L[[j]][2])]
            for (k in unlist(L[[j]][1])){
              ww[k] = ww[cind] - pt[cind] + pt[k]
              
            }
          }
          return (EquivClass(new, P, ww))
          
        }
      }
    }
    
    ww = w
    for(i in 1:length(L)){
      
      if (L[[i]][2] > 0){
        pt = P[,unlist(L[[i]][2])]
        cind = sort(unlist(L[[i]][1]))[1] 
        
        for (k in unlist(L[[i]][1])){
          ww[k]= ww[cind] - pt [cind] + pt [k]
        }
      }
    }
    return (list(L, ww))
  }
}


MinQuad <- function (f, n){ 
  l = ""
  z = ""
  for (i in 1:n){
    zi = paste0('z', as.character(i) ) 
    zi <- Var(zi)
  }
  for (i in 1:n){
    zi = paste0('z', as.character(i) ) 
    fzi = Deriv(f , zi)
    sympy (paste0("zi=simplify(",zi, ")" ))
    z = paste0(z,  zi)
    if (i != n) z = paste0(z, ",")
    fzi = sympy( paste0("fzi=simplify(",fzi, ")" ))
    l = paste0(l,fzi,  ",")
  }
  l = paste0(l, "z1")
  sympy (paste0("l = [", l, "]"))
  sympy (paste0("z = (", z, ")"))
  solution = sympy("sol = solve(l, z)")
  return(solution)}


unstring <- function (sol, n){
  sol = unlist(strsplit(sol, ':', fixed = TRUE ))
  sol = unlist(strsplit(sol, ',', fixed = TRUE ))
  sol = unlist(strsplit(sol, '}', fixed = TRUE ))
  indices = c(1:n)*2
  sol = as.numeric(sol[indices])
  sol [length(sol)] = sol [1]
  sol [1] = 0
  return (sol)
}


Cand <- function (newpt, P, height){
  m = ncol(P)
  n = nrow(P)
  infos = Infos(newpt,P)
  if (length(infos) !=2*m) return ("error, lengths do not match")
  ec = EquivClass(infos, P, rep(0, n))
  atlas = ec[[2]]
  if (length(ec[[1]]) ==1) return (NormalVec(atlas, height))
  atlas = as.list(ec[[2]])
  for (i in 1:length(ec[[1]])){
    for (ind in ec[[1]][[i]][[1]]){
      atlas [[ind]] = paste0("z", i,  "+",  as.character(atlas [[ind]])) 
    }
  }
  obj = ""
  for (i in 1:m){
    indp = sort(infos[[2*i-1]][[1]])[1] 
    indv = sort(infos[[2*i]][[1]])[1] 
    obj = paste0(obj , "(", atlas[[indp]], "-(", atlas[[indv]], ")+", as.character(P[indv, i]), "-(", as.character(P[indp, i]), "))^2")
    if (i !=m) obj = paste0(obj, "+")
  }
  sol = MinQuad(obj,length(ec[[1]]) ) 
  sol = unstring (sol, length(ec[[1]]))
  atl = ec[[2]]
  for (i in 1:length(ec[[1]])){
    for (ind in ec[[1]][[i]][[1]]){
      atl [ind] = atl[ind] + sol[i]
    }
  }
  return(NormalVec(atl , height)) 
}

# tropical distance
dtr <- function(l1, l2) max(l1-l2) - min(l1-l2) 

#sum of squares tropical distance
ssq <- function (pt , P) {
  ssq = 0
  for (i in 1:ncol(P)){
    x = P[,i]
    ssq = ssq + (dtr(x, pt))^2
  }
  return (ssq)
}

OneFrechet <- function (pt, P, height){
  # pt is the origin (any point of the kind rep(heigth))
  # P is a (m x n) matrix, m = dimension of the points, n =  number of points in the sample
  # P contains a sample of points on the columns
  if (IsFrechet(pt, P)) return (NormalVec(pt, height))
  m = ncol(P)
  n = length(pt)
  status = FALSE
  current = pt
  c = 0
  while ((status != TRUE) & (c<=10)){
    info = Infos(current, P)
    dtrs=c()
    allcandidates = list()
    for (i in 1:m){
      dtrs = c(dtrs, dtr(current, P[,i]))
      valley = unlist(info[[i*2]][[1]])
      candidates = powerSet(c(1:n)[-valley]) 
      for (j in 1:length(candidates)){
        candidates[[j]] = sort(union(unlist(candidates[[j]]), valley))
      }
      candidates = candidates[- length(candidates)]
      allcandidates = append(candidates, allcandidates) 
    }
    allcandidates = unique(allcandidates) # allcandidates contains indices towards which we might have to move 
    opti = list(current, ssq(current, P))
    for (i in 1: length(allcandidates)){
      incr = allcandidates[[i]]
      newpt = current
      moves = apply(P, 2, function (x) DescendPerturb(current, x, incr, c()))
      slope = sum(sapply (1:m, function(j) dtrs[j]*moves[1, j])) # moves[1,] contains the slope: 1, -1, 0
      if (slope <0) { # if slope < 0 the optimal point might change otherwise remains the current
        counter = sum(abs(moves[1,]))
        diss = c(moves[2,] , - slope/counter) # moves[2,] contains the distances
        diss = unique(diss[which(diss!=0.0)]) # length of the step towards negative direction
        slide = min(diss)
        for (k in 1:n){
          if (k %in% incr) newpt[k] = newpt[k] + slide
        }
        newpt = NormalVec(newpt, height)
        newssq = ssq(newpt, P)
        if (newssq < opti[[2]]) opti = list(newpt, newssq)
      }
    }
    newpt = opti[[1]]
    c = c+1
    if (all(current == newpt)) {status = TRUE}
    else if (identical (Infos(newpt, P), Infos(current, P))){
      #return (list(Infos(newpt, P), Infos(current, P)))
      cand = Cand (Infos(newpt, P), P, height)
      if (IsFrechet(cand, P)) {
        status = TRUE
        current = cand
      }
      else{
        status = IsFrechet(newpt, P)
        current = newpt
      }
    }
    else {
      status = IsFrechet(newpt, P)
      current = newpt
    }
  }
  if (IsFrechet(current, P) == FALSE) print("Failure")
  return (current)
}

# flat directions

FlatPerturb <- function (p1, p2, incr, decr) {
  dif = p1 - p2
  peaks = Infor(p1, p2, 1)
  valleys = Infor(p1, p2, 2)
  
  if (length(intersect(peaks,incr)) != 0){ep = 1} 
  else if (all(peaks %in% decr)) {ep = -1 }
  else {ep = 0}
  if (length(intersect(valleys,decr)) != 0) {ev = -1} 
  else if (all(valleys %in% incr)) {ev = 1 }
  else {ev = 0}
  
  maximum = max(dif)
  minimum = min(dif)
  bounds = c()
  
  for (index in 1: length(dif)){
    if ((ep == -1) & (index %in% incr)) {bounds = c(bounds, 0.5 * (maximum - dif[index]))}
    else if ((ep == 0) & (index %in% incr)) {bounds = c(bounds,  (maximum - dif[index]))}
    else if ((ep == -1) & ((index %in% decr) == FALSE)) {bounds = c(bounds, (maximum - dif[index]))}
    
    if ((ev == 1) & (index %in% decr)) {bounds = c(bounds, 0.5 * (dif[index] - minimum))}
    else if ((ev == 0) & (index %in% decr)) {bounds = c(bounds,  (dif[index] - minimum))}
    else if ((ev == 1) & (index %in% incr == FALSE)) {bounds = c(bounds, (dif[index] - minimum))}
  }
  if (length(bounds) == 0) {dist = 0}
  else {dist = min(bounds)}
  #print(bounds)
  return (c(dist, ep, ev))
}
#FlatPerturb(c(0, 0, 0, 0, 2, 0.4), c(1/5, 2/5, 2, 2/5, 2, 2), c(1, 2, 3, 5, 6), c())

closure <- function(Info , pos, ind){
  #If ind belongs to incr (pos=1) or decr (pos=2), the other indices that must belong too.
  o = c()
  new = c(ind)
  while (identical(o,new) == FALSE){
    o = new
    new = c()
    for (i in 1: length(Info)){
      if (length(intersect(o, unlist(Info [[i]][[pos]]))) != 0) {
        new = union(new, unlist(Info [[i]][[3-pos]]))
      }
    }
    new = union(new, o)
  }
  return (o)
}

closures <- function (Info, pos, n){
  #The list of all possible closures for incr (pos=1) or decr (pos=2).
  cls = lapply(1:n, function(i) closure(Info, pos, i))
  o = cls
  for (cl in o){
    cl = unlist (cl)
    rem = setdiff(1:n, cl)
    cands = powerSet(rem)
    for (cand in cands){
      if (length(cand) == 0) next
      newcl = as.numeric(cl)
      newcls = list()
      for (i in 1:length(cand)){
        newcls[[i]] = cls[[cand[[i]]]]
      }
      for (s in newcls){
        newcl = union(newcl, s)
      }
      o[[length(o)+1]] = newcl
    }
  }
  o = unique (o)
  for (j in length(o):1){
    o[[j]] = sort(round(o[[j]], 2))
    if (length(o[[j]])==n) {
      o = o[-c(j)]
    }
  }
  return (unique(o))
}


AddPerturb <- function (p, incr , decr, dist){
  output = p
  for (index in 1:length(p)){
    if (index %in% incr) output [index] = p[index] + dist
    else if (index %in% decr) output [index]= p[index] - dist
  }
  return (output)
}

FlatDirs <- function (p, P){
  # The directions with slope zero at a Frechet mean p. 
  # If p is a vertex of the polytope of Frechet means, status is true.
  n = nrow(P)
  infos = apply(P, 2, function(x) list(Infor(p, x, 1), Infor(p, x, 2)))
  candidates = closures (infos, 1, n)
  if (length(candidates) == 0) return (list(list(), list()))
  dists = c()
  for (i in 1:length(candidates)){
    cand = unlist(candidates[[i]])
    di = apply (P, 2, function (y) FlatPerturb(p,y,cand, c())[1])
    dists = c(dists, min(di))
  }
  return (list(candidates, dists))
}

IsVertex <- function (p, P){
  n = length(p)
  infos = apply(P, 2, function(x) list(Infor(p, x, 1), Infor(p, x, 2)))
  cl1 = closures(infos, 1, n)
  cl2 = closures(infos, 2, n)
  common = intersect(cl1, cl2)
  return(length(common) == 0)
}



OneVertex <- function (p, P, height){
  pnormal = NormalVec(p, height)
  pt = pnormal
  n = length(p)
  while ((IsVertex(pt, P)) != TRUE){
    infos = apply(P, 2, function(x) list(Infor(pt, x, 1), Infor(pt, x, 2)))
    vcl = closures(infos , 2, n)
    dirs = c()
    for (vc in vcl){
      if (length(vc) == 1) dirs = c(dirs,vc)
    }
    dirs = sort(dirs)
    dir = c(dirs[1])
    dist = min (unlist(apply (P, 2, function (y) FlatPerturb(p,y,c(), dir)[1])))
    pt = NormalVec(AddPerturb(pt, c(), dir, dist), height)
  }
  return (pt)
}

inList <- function (vec, list){
  if (length(list)==0) return (FALSE)
  for (j in 1:length(list)){
    listj = round(unlist(list[[j]]), 5)
    vec = round(vec, 5)
    if (all( listj == vec)) return (TRUE)
  }
  return (FALSE)
}

FMPolytope <- function (p, P, height){
  v = OneVertex(p, P, height)
  vertices = list(v)
  minsum = ssq (v, P)
  pnormal = NormalVec(v, height)
  current =list(pnormal)
  nextround = list(1)
  while (length(nextround) != 0){
    nextround = list()
    for (pt in current){
      flats = FlatDirs(pt, P)
      candidates = flats[[1]]
      if (length(candidates) == 0) return(list(minsum, length(vertices), vertices))
      dists = unlist(flats[[2]])
      for (index in 1:length(candidates)){
        candidate = candidates [[index]]
        dist = dists [index]
        pnew = NormalVec(AddPerturb(pt, candidate, c(), dist), height)
        
        if ((inList(pnew, vertices) == FALSE) && (inList(pnew, current) == FALSE) && (inList(pnew, nextround) == FALSE) && (IsVertex(pnew, P))){
          nextround[[length(nextround)+1]] =  pnew
          vertices[[length(vertices)+1]] =  pnew
        }
      }
    }
    current = nextround
  }
  return (list(tropDistance = minsum, numberFrechet = length(vertices), FMs = vertices))
}


Frechet <- function (P, heigth){
  n = nrow(P)
  pt = OneFrechet(rep(heigth, n), P, heigth)
  if (IsFrechet(pt, P)==FALSE) return("Could not find one Frechet Mean")
  else (return (list (oneFrechet = pt, FMPolytope = FMPolytope(pt, P, heigth))))
}


# examples




sp2 = cbind(c(1/5, 2/5, 2, 2/5, 2, 2), c(2, 2, 2, 2/5, 2/5, 1/5), c(2/5, 2/5, 2, 1/5, 2, 2))
current= OneFrechet(c(0,0,0,0,0,0), sp2, 2)

sp1 = cbind(c(9/100, 19/50, 19/50, 1, 19/50, 19/50, 1, 19/50, 1, 1), c(1, 1, 1, 1, 21/100, 57/100, 61/100, 57/100, 61/100, 61/100), c(31/50, 1, 49/50, 49/50, 1, 49/50, 49/50, 1, 1, 63/100), c(1, 1, 1, 47/100, 7/50, 4/5, 1, 4/5, 1, 1))
current= OneFrechet(c(0,0,0,0,0,0, 0, 0, 0, 0), sp1, 2)


skinny = cbind(c(0,0,0), c(0,2,4), c(0,5,1))
OneFrechet(c(0,0,0), skinny, 2)

fat = cbind(c(0,0,4), c(0,3,0), c(0,5,6))
OneFrechet(c(0,0,0), fat, 4)

undefined = cbind(c(0,0,0), c(0,448,449), c(0,452,256))
OneFrechet(c(0,0,0), undefined, 256)


p = c(1/5, 2/5, 2, 2/5, 2, 2)
sp2 = cbind(c(1/5, 2/5, 2, 2/5, 2, 2), c(2, 2, 2, 2/5, 2/5, 1/5), c(2/5, 2/5, 2, 1/5, 2, 2))
FMPolytope(p, sp2, 2)

sp1 = cbind(c(9/100, 19/50, 19/50, 1, 19/50, 19/50, 1, 19/50, 1, 1), c(1, 1, 1, 1, 21/100, 57/100, 61/100, 57/100, 61/100, 61/100), c(31/50, 1, 49/50, 49/50, 1, 49/50, 49/50, 1, 1, 63/100), c(1, 1, 1, 47/100, 7/50, 4/5, 1, 4/5, 1, 1))
p = c(48/25, 48/25, 48/25, 193/100, 8/5, 43/25, 2, 43/25, 2, 193/100)
fmp = FMPolytope(p, sp1, 2)

  
FMsp1 = Frechet(sp1, 2)
FMsp2 =Frechet(sp2, 2)
FMskinny = Frechet(skinny, 2)
  