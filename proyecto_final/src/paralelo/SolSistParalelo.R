
SolSistParalelo=function(n,m,H1,H2,Amod,rc,rb,iter,cluster,p=-1){
  
  # Si no se especifica, se hace un nodo por core
  if(p == -1 | abs(as.integer(p)) != p) p <- detectCores()
  
  # En la primera iteración se crea el cluster
  if(iter==0){
    cluster <- makeCluster(p)
    clusterEvalQ(cluster,{library(Matrix); NULL})
  }

  # Asignar las columnas lo más equitativamente posible entre los procesadores.
  np <- floor(seq(1,n+m,length.out=(p+1)))
  
  #######
  # Funciones básicas. Reciben una lista con los parámetros necesarios
  #######
  
  # Calcular en un nodo: Ai * Hi_inv * Ait = Ci
  fun1 <- function(...){
    ...$A%*%solve(...$Hs,t(...$A))
  }
  # Calcular en un nodo: LAi * rci = (Ai * Hi_inv) * r_i
  fun2 <- function(...){
    ...$A%*%solve(...$Hs,...$r)
  }
  # Calcular en un nodo: H_inv * rci - Lai_t * dlambda = H_inv * rci - H_inv_t * At * dlambda
  fun3 <- function(...){
    solve(...$Hs,...$r) - solve(t(...$Hs),t(...$A)%*%...$dl)
  }
  # Ejecutar una de las funciones 1 a 3. L debe ser una lista que contenga los argumentos necesarios
  # para llamar la función fun
  ejecutar.tarea <- function(fun, L){
    fun(L)
  }
  
  #######
  # Estructuras de datos para procesamiento en paralelo
  #######
  
  # Generar H en formato sparse en master
  H <- rBind(H1,H2)
  
  # Una lista de la misma longitud que el número de núcleos. En cada elemento tiene una lista
  # con la información que utilizará cada nodo. Así evitamos mandar todas las matrices a todos los nodos.
  AHpart <- lapply(1:p, function(t){
    is_first <- as.numeric(t==1)
    list(A = Amod[, (np[t]+(1-is_first)):np[t+1]],
         Hs = .sparseDiagonal((np[t+1]-np[t]+is_first),H@x[(np[t]+(1-is_first)):np[t+1]]),
         r = rc[(np[t]+(1-is_first)):np[t+1]])
  })
  
  #######
  # Cálculos
  #######
  
  # Calcula C en el master
  C <- Cholesky(-1*Reduce('+', clusterApply(cluster, AHpart, function(x) ejecutar.tarea(fun1,x))))
  
  # Calcula suma de Lai * rci. Cada producto en un nodo; la suma en el master
  suma1 <- Reduce('+', clusterApply(cluster, AHpart, function(x) ejecutar.tarea(fun2,x)))
  
  # Calcula dlambda = H_inv * (rb - suma) = (L * D * Lt)_inv * (rb - suma). Todo en el master
  dlambda <- solve(C,rb - suma1,'LDLt')
  
  # Calcula dzeta
  dzeta <- do.call(rBind, clusterApply(cluster, AHpart, function(x) ejecutar.tarea(fun3,c(x, dl=dlambda))))
  
  list(dy = dzeta[1:n], dx = dzeta[(n+1):(n+m)], dlambda = dlambda, cluster=cluster)

}


