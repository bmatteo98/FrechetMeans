library(Deriv)
library(rSymPy)  
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
  for (i in 1:ncol(P)){
    p2 = P[,i]
    peaks = Infor(p1,p2,1)
    valleys = Infor(p1,p2,2)
    Li = list(peaks, valleys)

    L[[i]] = Li
  }
  return(L)
}
p = c(9/5, 2, 2, 2, 2, 9/5)
p= c(8/5, 2, 1, 2, 2, 9/5)
p=c(1/7, 5/7, 1, 0, 1/2, 1)
p=c(2/5, 2/5, 2, -2/3, 14/15, 14/15)
sp2 = cbind(c(1/5, 2/5, 2, 2/5, 2, 2), c(2, 2, 2, 2/5, 2/5, 1/5), c(2/5, 2/5, 2, 1/5, 2, 2))
L = Infos(c(0,0,0,0,0,0),sp2)
#this one
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

w = rep(0,nrow(sp2))

res = EquivClass(L, sp2, w)


Cand <- function (newpt, P, height){
  m = ncol(P)
  n = nrow(P)
  infos = Infos(newpt,sp2)
  if (length(infos) !=2*m) return ("error, lengths do not match")
  ec = EquivClass(infos, P, rep(0, n))
  atlas = as.list(ec[[2]])
  if (length(ec[[1]]) ==1) return (NormalVec(atlas, height))
  for (i in 1:length(ec[[1]])){
    for (ind in ec[[1]][[i]][[1]]){
      atlas [[ind]] = paste0("z", i,  "+",  as.character(atlas [[ind]])) # z???
    }
  }
  obj = ""
  for (i in 1:m){
    indp = sort(infos[[2*i-1]][[1]])[1] 
    indv = sort(infos[[2*i]][[1]])[1] 
    obj = paste0(obj , "(", atlas[[indp]], "-(", atlas[[indv]], ")+", as.character(P[indv, i]), "-(", as.character(P[indp, i]), "))^2")
    if (i !=m) obj = paste0(obj, "+")
  }
  sol = MinQuad(obj,length(ec[[1]]) ) #expand(obj), [seq(cat(z, i), i = 1 .. nops(ec[1]))]
  sol = unstring (sol, length(ec[[1]]))
  atl = ec[[2]]
  for (i in 1:length(ec[[1]])){
    for (ind in ec[[1]][[i]][[1]]){
      atl [ind] = atl[ind] + sol[i]
    }
  }
  return(NormalVec(atl , height)) #map(x -> subs(sol, x), atlas), height ??
}
NormalVec <-function (p,height) p-max(p)+height

 Cand(p, sp2, 2)

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

library(Deriv)
f <- function (z1,z2) ((z2 - z1 + 2)^2 + (z1 - z2 + 8/5)^2 + (z2 - z1 + 9/5)^2)
f = expression((z2 - z1 + 2)^2 + (z1 - z2 + 8/5)^2 + (z2 - z1 + 9/5)^2)
f <- "(z3+0-z4+0+2-0.4)^2+(z1+0-z2+0+2-0.2)^2+(z3+0-z4+0+2-0.4)^2"
fi = strsplit(Deriv(f, c("z1","z2")), ';')
fi = substitute(Deriv(f, c("z1","z2")), list(e1 = e1[[1]], e2 = e2[[1]]))
fi = deparse(Deriv(f, c("z1","z2"), combine = 'cbind'))
library(nlsr)
evalf(D(f , 'z2'))
Deriv(f , 'z1')
nlsSimplify(fz1)
nlsDeriv(f, 'z1')
library(rSymPy)  
z1 <- Var("z1")
z2 <- Var("z2")
z3 <- Var("z3")
z4 <- Var("z4")
obj= sympy("fiz1= simplify((z3+0-z4+0+2-0.4)^2+(z1+0-z2+0+2-0.2)^2+(z3+0-z4+0+2-0.4)^2)")
sympy("fiz2 = simplify(2 * (z2 - z1 + 2) - 2 * (z1 - z2 + 8/5) + 2 * (z2 - z1 + 9/5))")
solution = sympy("sol = solve([fiz1, -fiz1, z1], z1, z2)")

f = expression((z2 - z1 + 2)^2 + (z1 - z2 + 8/5)^2 + (z2 - z1 + 9/5)^2)
fiz1 = D(f, 'z1')
fiz2 = D(f, 'z2')
z1 <- Var("z1")
z2 <- Var("z2")
fiz1 = Var(fiz1)
fiz1 =sympy("fiz1= simplify(fiz1) ")
fiz2 =sympy("fiz2 = simplify(fiz2) ")
solution = sympy("sol = solve([fiz1, fiz2, z1], z1, z2)")

MinQuad <- function (f, n){ 
  sympy ("l= []")
  sympy ("z = []")
  for (i in 1:n){
    zi = paste0('z', as.character(i) ) 
    zi <- Var(zi)
  }
  for (i in 1:n){
  zi = paste0('z', as.character(i) ) 
  fzi = Deriv(f , zi)
  sympy (paste0("zi=simplify(",zi, ")" ))
  z = sympy("z.append(zi)")
  fzi = sympy( paste0("fzi=simplify(",fzi, ")" ))
  l = sympy ("l.append(fzi)")
  }
  return(l)
  solution = sympy("sol = solve(l, z1, z2)")
}

MinQuad(f, 4)
library(caracas)
f <- "(z2 - z1 + 2)^2 + (z1 - z2 + 8/5)^2 + (z2 - z1 + 9/5)^2"
f = "2*(z2 - z3 + 8/5)^2 + (z1 - z3 + 9/5)^2"
f = "2*(z3 - z4 + 8/5)^2 + (z1 - z2 + 9/5)^2"
fiz1 = Deriv(f , 'z1')
sympy( paste0("fiz1=simplify(",fiz1, ))
fiz1 = Deriv(f , 'z2')
fiz1 = as_sym(fiz1)
fiz1 = as_sym(fiz2)
sol = solve(fiz1)

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