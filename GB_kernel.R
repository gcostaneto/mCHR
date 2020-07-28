# X a list of molecular markers

GB_kernel <-function(X){
  cat(paste0('------------------------------------','\n'))
  cat(paste0('Starting at ',Sys.time(),'\n'))
  
  WWG <- function(M, p){
    w <- scale(x = M, center = T, scale = F)
    
    S <- ((M==2)*1) * - rep(2*(1-p)^2, each=N) + ((M==1)*1) * rep(2*p*(1-p), each=N) + ((M==0)*1) * (-rep(2*p^2, each=N))
    
    WWl <- w %*% t(w)
    Ga <- WWl/(sum(diag(WWl))/N) + diag(1e-6, nrow(WWl))
    
    SSl <- S %*% t(S)
    Gd <- SSl/(sum(diag(SSl))/N)
    
    return(list(Ga=Ga,Gd=Gd,M=M,S=S))
  }
  Ga = list()
  Gd = list()
  
  if(!is.list(X)){
    M <- X
    N <- nrow(M)
    m <- ncol(M)
    p <- colMeans(M)/2
    return(WWG(M=M,p=p))
  }
  
  for(i in 1:length(X)){
    cat(paste0('------------------------------------','\n'))
    cat(paste0('Computing Kernel ',i,'\n'))
    M <- X[[i]]
    N <- nrow(M)
    m <- ncol(M)
    p <- colMeans(M)/2
    out = WWG(M=M,p=p)
    Ga[[i]]=out$Ga
    Gd[[i]]=out$Gd
  }
  cat(paste0('------------------------------------','\n'))
  cat(paste0('Organizing outputs','\n'))
  
  names(Ga) = paste0('Ga_',names(X))
  names(Gd) = paste0('Gd_',names(X))
  cat(paste0('------------------------------------','\n'))
  cat(paste0('Done at ',Sys.time(),'\n'))
  return(list(Ga = Ga, Gd = Gd))
}
