source("EDAFinal.R")

datos <- read.csv("test1.csv")


ejeEDARandom <- DEAFinal(datos, 66)
#cluster1 <- ejeEDARandom[[1]]
#cluster2 <- ejeEDARandom[[2]]
#cluster3 <- ejeEDARandom[[3]]
#cluster4 <- ejeEDARandom[[4]]
#jaccard <- vectorJaccard(cluster1, cluster2, cluster3, cluster4)
#jaccard

ejeEDARandom <- DEAFinal(datos, 66)
#cluster1 <- ejeEDARandom[[1]]
#cluster2 <- ejeEDARandom[[2]]
#cluster3 <- ejeEDARandom[[3]]
#cluster4 <- ejeEDARandom[[4]]
#jaccard <- vectorJaccard(cluster1, cluster2, cluster3, cluster4)
#jaccard
