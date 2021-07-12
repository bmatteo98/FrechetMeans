library(rje)

Infor <- function (l1,l2, pv){ # l1, l2 are vectors
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
NormalVec <-function (p,height) p-max(p)+height

sp2 = cbind(c(1/5, 2/5, 2, 2/5, 2, 2), c(2, 2, 2, 2/5, 2/5, 1/5), c(2/5, 2/5, 2, 1/5, 2, 2))
p = c(2/5, 2/5, 2, 1/5, 14/15, 14/15)
p = c(0, 0, 0, 0, 2, 0.4)
p = c(0.5, 1, 0, 2/5, 1, 2)
infos = apply(sp2, 2, function(x) list(Infor(p, x, 1), Infor(p, x, 2)))
# because we already have a frechet mean we cannot expect to descend more
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
FlatPerturb(c(0, 0, 0, 0, 2, 0.4), c(1/5, 2/5, 2, 2/5, 2, 2), c(1, 2, 3, 5, 6), c())

closure <- function(Info , pos, ind){
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

candidates = closures(infos, 1,6)

AddPerturb <- function (p, incr , decr, dist){
  output = p
  for (index in 1:length(p)){
    if (index %in% incr) output [index] = p[index] + dist
    else if (index %in% decr) output [index]= p[index] - dist
  }
  return (output)
}

FlatDirs <- function (p, P){
  n = nrow(P)
  infos = apply(P, 2, function(x) list(Infor(p, x, 1), Infor(p, x, 2)))
  candidates = closures (infos, 1, n)
  dists = c()
  for (i in 1:length(candidates)){
    cand = unlist(candidates[[i]])
    di = apply (P, 2, function (y) FlatPerturb(p,y,cand, c())[1])
    if (i == 23) {
      print(cand)
      print(di)
    }
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
  while (IsFrechet(pt, sample) != TRUE){
    infos = apply(P, 2, function(x) list(Infor(pt, x, 1), Infor(pt, x, 2)))
    vcl = closures(infos , 2, n)
    dirs = c()
    for (vc in vcl){
      if (length(vc) == 1) dirs = c(dirs,vc)
    }
    dirs = sort(dirs)
    dir = dirs[1]
    dist = min (apply (P, 2, function (y) FlatPerturb(p,y,c, dir)[1]))
    pt = NormalVec(AddPerturb(pt, c(), dir, dist), height)
  }
  return (pt)
}
OneVertex(p, sp2, 2)

c(1/5, 2/5, 0.4, 1/5, 1.400000000, 1/5, 1/5, 1/5, 1/5, 2/5, 1/5, 2/5, 1/5, 1.400000000, 1/5, 1/5, 1/5, 1/5, 1/5, 1/5, 1/5, 1/5, 2/5, 1/5, 1/5, 1/5, 0.4, 1/5, 1/5, 1/5, 1/5, 0.4, 1/5, 8/5)
