
SolSistParalelo=function(n,m,H1,H2,Amod,rc,rb,iter,cluster){
  require(dplyr)
  p <- detectCores()
  if(iter==0){
    cluster <- makeCluster(p)
    clusterEvalQ(cluster,{library(Matrix); NULL})
  }

  # Asignamos las columnas lo más equitativamente posible entre los procesadores.
  np <- floor(seq(1,n+m,length.out=(p+1)))
  # Calcular Ai * Hi_inv * Ait = Ci
  fun1 <- function(...){
    ...$A%*%solve(...$Hs,t(...$A))
  }
  # Calcular Ai * Hi_inv * r_i
  fun2 <- function(...){
    ...$A%*%solve(...$Hs,...$r)
  }
  fun3 <- function(...){
    solve(...$Hs,...$r) - solve(t(...$Hs),t(...$A)%*%...$dl)
  }
  ejecutar.tarea <- function(fun, L){
    fun(L)
  }
  
  H <- rBind(H1,H2)
  AHpart <- lapply(1:p, function(t){
    is_first <- as.numeric(t==1)
    list(A = Amod[, (np[t]+(1-is_first)):np[t+1]],
         Hs = .sparseDiagonal((np[t+1]-np[t]+is_first),H@x[(np[t]+(1-is_first)):np[t+1]]),
         r = rc[(np[t]+(1-is_first)):np[t+1]])
  })
  

  C <- Cholesky(-1*Reduce('+', clusterApply(cluster, AHpart, function(x) ejecutar.tarea(fun1,x))))

  suma1 <- -1*Reduce('+', clusterApply(cluster, AHpart, function(x) ejecutar.tarea(fun2,x)))
  
  dlambda <- solve(C,rb+suma1,'LDLt')
  
  dzeta <- do.call(rBind, clusterApply(cluster, AHpart, function(x) ejecutar.tarea(fun3,c(x, dl=dlambda))))
  
  list(dzeta[1:n], dzeta[(n+1):(n+m)], dlambda, cluster)

}


