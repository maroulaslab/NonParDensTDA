# IMPORTS
#require(TDA)
#require(mvtnorm)
#require(combinat)
#require(plyr)
#require(clue)
#---------------

# auxillary functions
in_region <- function(p,xlims = c(0,Inf), ylims = c(0, Inf)){
  # p: numeric representing a sample from a gaussian mixutre
  
  c1 = p[1] >= xlims[1]
  c2 = p[1] <= xlims[2]
  c3 = p[2] >= ylims[1]
  c4 = p[2] <= ylims[2]
  
  return(c1 & c2 & c3 * c4)
}

rg_samp <- function(mu,sig){
  
  I = diag(2)
  samp = mvtnorm::rmvnorm(1, mean = mu, sigma = sig*I )
  
  while(!in_region(samp)){
    samp = mvtnorm::rmvnorm(1, mean = mu, sigma = sig*I )
  }
  
  return(samp)
}

rmg_samp <- function(mat,si){
  proj.points = cbind(mat[,1],rep(0,nrow(mat)))
  comp.ind = sample(1:nrow(proj.points),1)
  samp = rg_samp(mu = proj.points[comp.ind,],sig = si)
  return(samp)
}

dist_mat <- function(upper.mat,obs.mat){
  mat = matrix(nrow = nrow(upper.mat),ncol = nrow(obs.mat))
  
  dist_sq <- function(vec){
    return(abs(sum(vec**2)))
  }
  
  for(i in 1:nrow(upper.mat)){
    for(j in 1:nrow(obs.mat)){
      mat[i,j] = dist_sq(upper.mat[i,]-obs.mat[j,])
    }
  }
  
  return(mat)
}

min_cost_injection <- function(upper.mat,obs.mat,n){
  d.mat = dist_mat(upper.mat = upper.mat,obs.mat = obs.mat)
  
  if(n == nrow(upper.mat))
  {combos = matrix(combinat::combn(1:nrow(upper.mat),m=n),ncol = 1)}
  else{
    combos = combn(1:nrow(upper.mat),m=n)}
  
  if(nrow(combos)==1){inj.costs = alply(combos,2,function(x){matrix(d.mat[x,],ncol = length(d.mat[x,]))})}
  else{inj.costs = plyr::alply(combos,2,function(x){d.mat[x,]})}
  hungarian = lapply(inj.costs,function(x){as.list(clue::solve_LSAP(x))})
  
  grab_array_entry <- function(a,l1 = as.list(1:nrow(a)),l2){
    r = mapply(function(x,y){a[x,y]},l1,l2, SIMPLIFY = FALSE)
  }
  
  optimal.injs = mapply(function(x,y){grab_array_entry(a=x,l2 = y)},inj.costs,hungarian,SIMPLIFY = FALSE)
  total.inj.costs = sapply(optimal.injs,function(x){sum(unlist(x))})
  min.inj.col = which.min(total.inj.costs)
  
  upper.rows = combos[,min.inj.col]
  obs.rows = unlist(hungarian[[min.inj.col]])
  return(matrix(cbind(upper.rows,obs.rows),ncol=2))
}

rg_eval <- function(p,mu,sig,normalize = TRUE){
  v = mvtnorm::dmvnorm(x = p, mean = mu, sigma = sig*diag(2))
  norm = ifelse(normalize,mvtnorm::pmvnorm(mean = mu,upper = c(Inf,Inf),lower = c(0,0),sigma = sig*diag(2)),1)
  return(v/norm)
}

rmg_eval <- function(po,mus,si){
  mus = cbind(mus[,1],rep(0,nrow(mus)))
  v = apply(mus,1,function(x){rg_eval(p=po,mu = x, sig = si, normalize = FALSE)})
  s = sum(v)
  return((2*s)/nrow(mus))
}

track_eval <- function(pairings,upper, lower, obs,a.probs,l.probs,si){
  
  # compute Q 
  appears = rep(0,nrow(upper))
  up.inds = 1:nrow(upper)
  there = which(up.inds %in% pairings[,1])
  appears = replace(appears,there,1)
  a = a.probs[which(as.logical(appears))]
  na = a.probs[which(!as.logical(appears))]
  na = 1-na
  na = c(1,na)
  a = Reduce('*',as.list(a))
  na = Reduce('*',as.list(na))
  Q = a*na
  
  upper.l = lapply(1:nrow(upper),function(x){upper[x,]})
  obs.l = lapply(1:nrow(obs),function(x){obs[x,]})
  upper.pres = upper.l[pairings[,1]]
  obs.pres = obs.l[pairings[,2]]
  
  u.ds = mapply(function(x,y){rg_eval(mu = x, p = y, sig = si)},upper.pres,obs.pres,SIMPLIFY = FALSE)
  u.res = Reduce('*',u.ds)
  
  obs.lower = obs.l[-pairings[,2]]
  obs.lower = do.call(rbind,obs.lower)
  
  l.ds = apply(obs.lower,1,function(x){rmg_eval(mus = lower,po = x, si = si )})
  l.res = Reduce('*',as.list(l.ds))
  
  lcard = l.probs[nrow(obs.lower)]
  
  return(Q*u.res*lcard*l.res)
  
}



# TDA Kernel class
#' Create a persistence diagram kernel density.
#'
#' A reference class to represent the persistence diagram kernel from Maroulas, Mike, and Oballe 2019.
#' It allows sampling and approximate evaluation. 
#' @field center A diagram created from the TDA package. This is the center of the kernel.
#' @field sigma A numeric. This is the bandwidth of the kernel.
#' @field lower.card.probs A vector of probabilities for to serve as the cardinality distribution of the lower portion of the diagram. The probability at index i is the probability of observing i points from the lower diagram.
#' @field upper A matrix whose rows are points in the upper portion of the kernel in birth-persistence coordinates. Created automatically on initialization.
#' @field lower  A matrix whose rows are points in the lower portion of the kernel in birth-persistence coordinates. Created automatically on initialization.
#' @field upper.probs A vector of probabilities of appearances for each point in the upper diagram. The ith entry is the probability of appearance corresponding to row i of upper. Computed automatically on intialization.
#' @export

TDAKernel <- setRefClass('TDAKernel',
                         
                         fields = list(center = 'list', sigma = 'numeric', upper = 'matrix', lower = 'matrix', upper.probs = 'numeric',lower.card.probs = 'numeric'),
                         # center: center of the kernel; should be a diagram returned from TDA package
                         # sigma: sigma parameter from Mike and Maroulas
                         # upper: matrix whose rows are points in the upper diagram
                         # lower: matrix whose rows are points in the lower diagram
                         #    points in upper and lower are given in birth persistence coordinates
                         
                         methods = list(
                           initialize = function(center,sigma,lower.card.probs, fast.eval = TRUE){
                             "Creates the TDAKernel object. Requires input of center, sigma, and lower.card.probs. User can also designate fast.eval = FALSE to assume all upper.probs are 1."
                             .self$center = center
                             .self$sigma = sigma
                             .self$lower.card.probs = lower.card.probs
                             D = center$diagram
                             TD = cbind(D[,2], D[,3]-D[,2])
                             .self$upper = matrix(TD[which(TD[,2] > sigma),],ncol = 2)
                             if(nrow(.self$upper) == 0){print('empty upper')}
                             .self$lower = matrix(TD[which(TD[,2] <= sigma),],ncol = 2)
                             
                             I = diag(1,2)
                             if(fast.eval){.self$upper.probs = rep(1,nrow(.self$upper))}
                             else{.self$upper.probs = apply(.self$upper,1,function(x){pmvnorm(mean = x,upper = c(Inf,Inf),lower = c(-Inf,0),sigma = sigma*I)})}
                             
                           },
                           
                           sample = function(n){
                             "Draw n samples from the kernel."
                             draw <- function(){
                               appears = sapply(.self$upper.probs,function(x){rbinom(1,1,x)})
                               appears = as.logical(appears)
                               appears.upper = .self$upper[which(appears),]
                               appears.upper = matrix(appears.upper,ncol=2)
                               card = sample.int(length(.self$lower.card.probs),size=1,prob = .self$lower.card.probs)
                               if(nrow(appears.upper) > 0){
                                 u.samp = lapply(1:nrow(appears.upper),function(x){rg_samp(mu = appears.upper[x,],sig = .self$sigma)})
                                 u.samp = do.call(rbind,u.samp)
                               }
                               else{u.samp = matrix(nrow = 0, ncol = 2)}
                               l.samp = lapply(1:card,function(x){rmg_samp(mat = .self$lower,si = .self$sigma)})
                               l.samp = do.call(rbind,l.samp)
                               return(rbind(u.samp,l.samp))
                             }
                             
                             return(lapply(1:n,function(x){draw()}))
                           },
                           
                           eval = function(D){
                             "Approximate evaluation of a diagram D (a diagram from the TDA package) in the kernel. Evaluation is the sum over highest probability matchings for all possible upper diagram cardinalities."
                             # D: diagram from TDA package
                             
                             if(is.list(D)){
                               D = D$diagram
                               D = cbind(D[,2],D[,3]-D[,2])
                             }
                             n.Du = nrow(.self$upper)
                             from.upper = 1:n.Du
                             
                             matchings = lapply(from.upper,function(x){min_cost_injection(upper.mat = .self$upper,obs.mat = D,n = x)})
                             results = lapply(matchings,function(x){track_eval(pairings = x, upper = .self$upper,lower = .self$lower, obs = D,si = .self$sigma, a.probs = .self$upper.probs, l.probs = .self$lower.card.probs)})
                             
                             R = Reduce('+', results)
                             return(R)
                             
                             
                           }
                         ))