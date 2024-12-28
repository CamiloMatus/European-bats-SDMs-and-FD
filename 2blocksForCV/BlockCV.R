#Here spatial blocks are created by species for cross-validation.

require(blockCV)
require(sf)
library(doParallel)
require(tidyverse)
cl<-makeCluster(4);doParallel::registerDoParallel(cl)


db <- read_sf("../1dataBase/bioHistoric_SampRate.shp")
dbPres <- read.csv("../1dataBase/speciesPresences.csv")
db2 <- merge(db,dbPres,by="idg")%>%sf::st_centroid()
spps <- names(db2)[22:59]

foreach(i = spps,.packages = "blockCV",.combine = "return" ) %dopar% {
db.h <- db2[,i]
sb1 <- cv_spatial(x = db.h,
                  column = i,
                  size = 200000,
                  k = 5,plot = FALSE,
                  selection = "random",
                  iteration = 1000)

db.h2 <- data.frame(sb1$biomod_table)
db.h2$idg <- db2$idg
write.csv(db.h2,file = paste("blocksCV",i,"csv",sep = "."),row.names = F)
}
