

Infor <- function (l1,l2, pv){ # l1, l2 are vectors
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
  if (pv == 1) return (i)
  if (pv == 2) return (j)
}


Peakpairs(c(9/5, 2, 2, 2, 2, 9/5), c(2, 2, 2, 2/5, 2/5, 1/5))
p1 = c(9/5, 2, 2, 2, 2, 9/5)
p2 = c(2, 2, 2, 2/5, 2/5, 1/5)
incr = c()
decr = c(3)
DescendPerturb <- function (p1, p2, incr, decr) {
  dif = p1 - p2
  peaks = Infor(p1, p2, 1)
  valleys = Infor(p1, p2, 2)

  if (length(intersect(peaks,incr)) != 0){ep = 1} 
  else if (all(peaks %in% decr)) {ep = -1 }
  else {ep = 0}
  if (length(intersect(valleys,decr)) != 0) {ev = -1} #??
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
  print(bounds)
  if (length(bounds) == 0) {dist = 0}
  else {dist = min(bounds)}

  return (c(slope, dist, ep, ev))
}
DescendPerturb(p1, p2, incr, decr)
