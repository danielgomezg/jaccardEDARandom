silueta <- function(datos, individuo, numCluster, matrizPearson){
  promSilueta <- NULL
  resultadoCluster <- NULL
  distancia <- NULL
  items <- NULL
  sumatoria <- NULL
  distanciaOtroCLuster <- NULL
  
  #cantidad de elementos que tiene un cluster
  for (z in 1:numCluster) {
    items <- c(items, length(individuo[[2]][[z]]))
  }
  
  for (i in 1:numCluster) {
    #cantidad de elementos que tiene un cluster
    #items <- length(individuo[[2]][[i]])
    
    
    #si items tiene valor 1 el calculo de homogeneidad es 1
    if(items[i] == 1){
      promSilueta <- c(promSilueta, 1)
      
    }else{
      for (j in 1:items[i]) {
        
        #obtener dsitancia promedio del item a los demas items del cluster
        for(k in  1:items[i]){
          if(j != k){
            distancia <- c(distancia, matrizPearson[individuo[[2]][[i]][[j]], individuo[[2]][[i]][[k]]])
          }
        }
        
        distanciaCluster <- mean(distancia) 
        #resultadoCluster <- c(resultadoCluster, resultado)
        distancia <- NULL
        
        #Obetiendo bi
        for (w in 1:numCluster) {
          if(w != i){
            for (y in 1:items[w]) {
              distancia <- c(distancia, matrizPearson[individuo[[2]][[i]][[j]], individuo[[2]][[w]][[y]]])
            }
            distanciaOtroCLuster <- c(distanciaOtroCLuster, mean(distancia))
            distancia <- NULL
          }
        }
        
        clusterMasCercano <- max(distanciaOtroCLuster)
        resultado <- (clusterMasCercano - distanciaCluster) / max(clusterMasCercano, distanciaCluster)
        sumatoria <- c(sumatoria, resultado)
      }
      #se guarda el resultado del cluster
      promSilueta <- c(promSilueta, mean(sumatoria))
      sumatoria <- NULL
    }
  }
  return(mean(promSilueta))
}

homogeneidad <- function(datos, individuo, numCluster, matrizPearson){
  #calcular la matriz de distancia
  promHomogeneidad <- NULL
  resultadoCluster <- NULL
  distancia <- NULL
  
  for (i in 1:numCluster) {
    #cantidad de elementos que tiene un cluster
    items <- length(individuo[[2]][[i]])
    
    #si items tiene valor 1 el calculo de homogeneidad es 1
    if(items == 1){
      promHomogeneidad <- c(promHomogeneidad, 1)
    }
    
    else{
      for (j in 1:items) {
        
        #obtener dsitancia promedio del item a los demas items del cluster
        for(k in  1:items){
          if(j != k){
            distancia <- c(distancia, matrizPearson[individuo[[2]][[i]][[j]], individuo[[2]][[i]][[k]]])
          }
        }
        
        resultado <- mean(distancia) / ((items * (items - 1)) / 2)
        resultadoCluster <- c(resultadoCluster, resultado)
        distancia <- NULL
      }
      #se guarda el resultado del cluster
      promHomogeneidad <- c(promHomogeneidad, sum(resultadoCluster))
      resultadoCluster <- NULL
    }
  }
  
  return(sum(promHomogeneidad)/numCluster)
}


mejor <- function(fitness, cantIndividuos){
  maxi <- max(fitness)
  for(i in 1:cantIndividuos){
    if(fitness[i] == maxi){
      resultado <- c(fitness[i], i)
      return(resultado)
    }
  }
  return(0)
}


vectorJaccard <- function(a, b, c, d){
  VectorFinal <- 1:1000
  for (i in 1:length(a)) {
    posicion <- a[i]
    VectorFinal[posicion] <- 1
  }
  for (i in 1:length(b)) {
    posicion <- b[i]
    VectorFinal[posicion] <- 2
  }
  for (i in 1:length(c)) {
    posicion <- c[i]
    VectorFinal[posicion] <- 3
  }
  for (i in 1:length(d)) {
    posicion <- d[i]
    VectorFinal[posicion] <- 4
  }
  return(VectorFinal)
}


