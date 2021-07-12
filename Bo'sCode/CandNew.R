library(Deriv)
library(rSymPy)  

NormalVec <-function (p,height) p-max(p)+height

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

Infos <- function(p1, P){
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
  infos = Infos(newpt,sp2)
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

p = c(9/5, 2, 2, 2, 2, 9/5)
p= c(8/5, 2, 1, 2, 2, 9/5)
p=c(2/5, 2/5, 2, -2/3, 14/15, 14/15)
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

# in these cases the solution to the system of derivetives =0 is not unique
# this should not happen : the function returns an error
p=c(1/7, 5/7, 1, 0, 1/2, 1)
p = c(1,2,3,1,2,3)

sp2 = cbind(c(1/5, 2/5, 2, 2/5, 2, 2), c(2, 2, 2, 2/5, 2/5, 1/5), c(2/5, 2/5, 2, 1/5, 2, 2))


Cand(p, sp2, 2)

