source("funcionesEDA.R")
#library("optparse")

#Funcion principal del algoritmo Evolucion diferencial
#DEA <- function(datos, tamanoP, fMutacion, CR){
DEAFinal <- function(datos, tamanoP){
  #datos es el dataset, numCluster cantidad de cluste que se van a crear, numGeneraciones cantidad de generaciones que realizara el algoritmo
  #tamanoP es el tamano de la poblacion, fMutacion es el factor de mutacion y CR es factor de cruce 
  tiempo <- proc.time()
  numCluster <- 4
  numGeneraciones <- 200
  #fMutacion <- round(runif(1, 0.3, 0.9), 1)
  #EDA
  fMutacion <- 0.5012
  CR <- 0.8808 # [0.5, 0.9]
  #fMutacion <- round(runif(1, 0.3, 0.9), 4)
  #CR <- round(runif(1, 0.4, 0.9), 4)
  #fMutacionInicial <- fMutacion
  #CRInicial <- CR
  mejorIndividuoNuevo <- 0
  
  maximo <- round(0.05 * numGeneraciones, 0)
  #maximo <- 10
  numeroSinMejora <- 0
  numeroConMejora <- 0
  mutacion <- NULL
  crMax <- NULL
  crMin <- NULL
  contador <- 1
  contadorCR <- 1
  
  #Grafico
  #matrizBox <- numeric(tamanoP * (numGeneraciones + 2)) 
  #dim(matrizBox) <- c(tamanoP, (numGeneraciones + 2)) 
  matrizBox2 <- numeric(tamanoP * (numGeneraciones + 1)) 
  dim(matrizBox2) <- c(tamanoP, (numGeneraciones + 1)) 
  
  g <- 1
  bestSolucion <- NULL
  pos <- NULL
  
  #obtener cantidad de muestras y genes
  dimensionesDatos <- dim(datos)
  features <- dimensionesDatos[1]
  cantDatos <- dimensionesDatos[2]
  
  #Crear lista de la poblacion
  pNew <- list()
  pNew <- vector("list", length = tamanoP)
  
  #Crear matriz de distancia
  matrizPearson <- matPearson(datos)
  
  
  #Contador para obtener evaluaciones
  contadorEva <- 0
  
  #crear poblacion inicial
  p <- setup(datos, numCluster, features, tamanoP, cantDatos, g)
  
  #Grafico de la poblacion incial 
  #matrizBox[,1] <- 0.9
  #matrizBox[,2] <- individuosRend(datos, p, numCluster, tamanoP, matrizPearson)
  #matrizBox2[,1] <- matrizBox[,2]
  
  #EDA
  #Rendimientos de los individuos
  RendimientoIndividuos <- individuosRend(datos, p, numCluster, tamanoP, matrizPearson)
  mejorIndividuoActual <- bestIndividuo2(RendimientoIndividuos)
  
  #matrizBox2[,1] <- individuosRend(datos, p, numCluster, tamanoP, matrizPearson)
  matrizBox2[,1] <- RendimientoIndividuos
  
  #mejorIndividuoActual <- bestIndividuo(datos, p, features, numCluster, tamanoP, matrizPearson)
  #mejorIndividuoActual <- mejor(matrizBox2[,1], tamanoP)
  print(paste0(( c( g, mejorIndividuoActual[1], fMutacion, "con: ", numeroConMejora, "sin: ", numeroSinMejora, maximo) )))
  
  promedio <- NULL

  while (g <= numGeneraciones) {
    
    #cr
    if((0.1 * numGeneraciones) < g){
      #fitnessPoblacion <- individuosRend(datos, p, numCluster, tamanoP, matrizPearson)
      #PromFitness <- sum(fitnessPoblacion) / length(fitnessPoblacion)
      PromFitness <- sum(RendimientoIndividuos) / length(RendimientoIndividuos)
    }
    
    #mutacion
    if(g %% maximo == 1 && g != 1){
      mejorIndividuoActual <- mejorIndividuoNuevo
      #mejorIndividuoActual <- bestIndividuo2(RendimientoIndividuos)
      #mejorIndividuoActual <- bestIndividuo(datos, p, features, numCluster, tamanoP, matrizPearson)
      print(paste0((c("Individuo Actual", mejorIndividuoActual[1], g, fMutacion) )))
    }
    
    #Grafico de CR
    limiteMaxCR <- 0
    limiteMinCR <- 0
    
    for(i in 1:tamanoP){
      #Se obtienen los individuos aleatorios 
      individuosRandom <- sample(1:tamanoP, 6, replace = FALSE)
      
      #CR
      if((0.1 * numGeneraciones) < g){
        #if(fitnessPoblacion[i] >= PromFitness){
        if(RendimientoIndividuos[i] >= PromFitness){
          
          #experimento cr
          if(contadorCR >= 4 && g > (0.55 * g)){
            CR <- 0.5 + ((0.9 - 0.5) * (g / numGeneraciones))
          }else{
          CR <- 0.5 - ((0.5 - 0.1) * (g / numGeneraciones))    }
          
          #Grafico de CR
          if(limiteMinCR == 0){
            crMin <- c(crMin, CR)
            limiteMinCR <- 1
            }
        }
        else{
          CR <- 0.5 + ((0.9 - 0.5) * (g / numGeneraciones))
          
          #Grafico de CR
          if(limiteMaxCR == 0){
            crMax <- c(crMax, CR)
            limiteMaxCR <- 1
          }
        }
      }else{
        
        #Grafico de CR
        if(limiteMinCR == 0){
          crMin <- c(crMin, CR)
          crMax <- c(crMax, CR)
          limiteMinCR <- 1
        }
      }
      if(i != individuosRandom[1] && i != individuosRandom[2] && i != individuosRandom[3]){
        #se realiza el crossover y se obtiene el mejoor individuo
        iPrueba <- crossover(p[[i]], p[[individuosRandom[1]]], p[[individuosRandom[2]]], p[[individuosRandom[3]]], datos, fMutacion, CR, numCluster, features, cantDatos, g)
        iPruebaRend <- silueta(datos, iPrueba, numCluster, matrizPearson)
        
        contadorEva <- contadorEva + 1
        
        #Selecciono el mejor individuo
        if(RendimientoIndividuos[i] > iPruebaRend){
          pNew[[i]] <- p[[i]]
        }else{
          pNew[[i]] <- iPrueba
          RendimientoIndividuos[i] <- iPruebaRend
        }
        #iNew <- reemplazo(datos, p[[i]], iPrueba, features, numCluster, matrizPearson)
        #pNew[[i]] <- iNew
      }else{
        #se realiza el crossover y se obtiene el mejoor individuo
        iPrueba <- crossover(p[[i]], p[[individuosRandom[4]]], p[[individuosRandom[5]]], p[[individuosRandom[6]]], datos, fMutacion, CR, numCluster, features, cantDatos, g)
        iPruebaRend <- silueta(datos, iPrueba, numCluster, matrizPearson)
        
        contadorEva <- contadorEva + 1
        
        #Selecciono el mejor individuo
        if(RendimientoIndividuos[i] > iPruebaRend){
          pNew[[i]] <- p[[i]]
        }else{
          pNew[[i]] <- iPrueba
          RendimientoIndividuos[i] <- iPruebaRend
        }
        #iNew <- reemplazo(datos, p[[i]], iPrueba, features, numCluster, matrizPearson)
        #pNew[[i]] <- iNew
      }
    }
    
    #Grafica de los individuos de la  generacion nueva
    #matrizBox[,g + 2] <- individuosRend(datos, pNew, numCluster, tamanoP, matrizPearson)
    #matrizBox2[,g + 1] <- matrizBox[,g + 2]
    #matrizBox2[,g + 1] <- individuosRend(datos, pNew, numCluster, tamanoP, matrizPearson)
    matrizBox2[,g + 1] <- RendimientoIndividuos
    
    
    #EDA 
    #fMutacion
    mutacion <- c(mutacion, fMutacion)
    if(g %% maximo == 0){
      mejorIndividuoNuevo <- bestIndividuo2(RendimientoIndividuos)
      #mejorIndividuoNuevo <- bestIndividuo(datos, pNew, features, numCluster, tamanoP, matrizPearson)
      #mejorIndividuoNuevo <- mejor(matrizBox2[, g + 1], tamanoP)
      print(paste0((c("Individuo Nuevo", mejorIndividuoNuevo[1], g, fMutacion) )))
    }
    
    if(contador == maximo){
      contador <- 1
      if(mejorIndividuoNuevo[1] > mejorIndividuoActual[1]){
        fMutacion <- 0.5 - (0.4 * (g / numGeneraciones))
        contadorCR <- 1
      }
      else{
        fMutacion <- 0.5 + (0.4 * (g / numGeneraciones))
        
        #experimento cr
        if(g > (0.55 * numGeneraciones)){
            contadorCR <- contadorCR + 1   
          }
      }
    }
    else{
      contador <- contador + 1
    }
    
    p <- pNew
    pNew <- list()
    pNew <- vector("list", length = tamanoP)
    g <- g + 1
    
  }
  
  generaciones <- 1:numGeneraciones
  generacionesCRMax <- 1:length(crMax)
  
  par(mar = c(5, 4, 4, 4) + 0.25)
  boxplot(x = matrizBox2, y = matrizBox2, xlab = "Generación ", ylab = "Valor del fitness")
  par(new = TRUE)
  plot(generaciones , crMin , type = "b", col = "green", lwd = 1.3, lty = 3, pch = 2, axes = FALSE, bty = "n", xlab = "", ylab = "")
  lines(generaciones ,mutacion , type = "b", col = "red", lwd = 1.3, lty = 3, pch = 20)
  lines(generacionesCRMax ,crMax , type = "b", col = "green", lwd = 1.3, lty = 3, pch = 2)
  L = legend(x = 'bottomleft', legend = c(0.5012,0.8808), col=2:3, lty=c(3,3), ncol=1, bty='n', x.intersp=0.5, pch=c(20,2), inset=0.02)
  
  # use position data of previous legend to draw legend with invisble lines and points but with labels and box. x.intersp controls distance between lines and labels
  legend(x = L$rect$left, y = L$rect$top, legend = c("Factor de mutación", "Cruce"), col=rep(NA,2), lty=c(1,1), ncol=1, x.intersp = 4, bg = NA)
  
  axis(4)
  mtext("Valor Cruce y Mutación", side = 4, line = 3, col = 1)
  
  print(paste0((proc.time() - tiempo)))
  
  print(paste0(c("Numero de evaluaciones: ", contadorEva)))
  
  #mejorSolucion <- bestIndividuo(datos, p, features, numCluster, tamanoP, matrizPearson)
  mejorSolucion <- bestIndividuo2(RendimientoIndividuos)
  print(paste0(c(mejorSolucion[1], mejorSolucion[2])))
  return(p)
}


#option_list = list(
#  make_option(c("--seed"), type = "integer"),
#  make_option(c("--fMutacion"), type = "integer"),
#  make_option(c("--CR"), type = "double"),
#  make_option(c("--tamanoP"), type = "double"),
#  make_option(c("--tries"), type = "integer"),
#  make_option(c("--time"), type = "integer"),
#  make_option(c("--quiet"), type = "logical"),
#  make_option(c("-i", "--input"), type = "character")
#);

#opt_parser = OptionParser(option_list=option_list);
#opt = parse_args(opt_parser);
#set.seed(opt$seed)

#datos <- read.csv(opt$input)

#resultado <- DE(datos ,opt$tamanoP, opt$fMutacion, opt$CR)

#cat(resultado)