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
p = c(0, 0, 0, 0, 2, 0.4)
p = c(2/5, 2/5, 2, 1/5, 14/15, 14/15)
sp2 = cbind(c(1/5, 2/5, 2, 2/5, 2, 2), c(2, 2, 2, 2/5, 2/5, 1/5), c(2/5, 2/5, 2, 1/5, 2, 2))
p = c(0.5, 1, 0, 2/5, 1, 2)
infos = apply(sp2, 2, function(x) list(Infor(p, x, 1), Infor(p, x, 2)))

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

#old
closures <- function (Info, pos, n){
  
  cls = lapply(1:n, function(i) closure(Info, pos, i))
  o = cls
  for (cl in o){
    cl = unlist (cl)
    cand = setdiff(1:n, cl)

      if (length(cand) == 0) next
      newcl = as.numeric(cl)
      newcls = list()
      for (i in 1:length(cand)){
      newcls[[i]] = cls[[cand[[i]]]]
      }
      for (s in newcls){
        newcl = union(newcl, s)
      }
      print(c("newcl", newcl))
      o[[length(o)+1]] = newcl
      print(c("o", o))
  }

  o = unique (o)
  for (j in length(o):1){
    if (length(o[[j]])==n) {
      o = o[-c(j)]
      }
  }
  return (o)
}



#old
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
    o[[j]] = sort(o[[j]])
    if (length(o[[j]])==n) {
      o = o[-c(j)]
    }
  }
  return (unique(o))
}

#this

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
  return (unique(o,  2))
}

closures(infos, 1,6)
o =closures(infos, 1,6)
