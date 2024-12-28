require(embarcadero)
require(PresenceAbsence)
require(sf)
require(dplyr)
require(ecospat)
require(foreign)
source("funBART_HyperOpti.R")


#reading data
dbPres <- read.csv("../1dataBase/speciesPresences.csv")
biosSampRate <- read.dbf("../1dataBase/bioHistoric_SampRate.dbf")
explainVars <- biosSampRate[,c("idg","bio1","bio2","bio7","bio8","bio12","bio15","bio18","sampRate")]
db0 <- merge(dbPres,explainVars,by="idg")



species <- names(dbPres)[1:38]
predictors <- c("idg","bio1","bio2","bio7","bio8","bio12","bio15","bio18","sampRate")



for(i in 1:length(species)){

# library(doParallel)
# cl<-makeCluster(4);doParallel::registerDoParallel(cl)  
# foreach(i = 1:length(species),
#         .packages = c("embarcadero","ecospat","raster","PresenceAbsence"),
#         .combine = "return" ) %dopar% {
  spp.h <- species[i]
  blocks.h <-  read.csv(paste("../2blocksForCV/blocksCV.",spp.h,".csv",sep = ""))
  db <- merge(db0,blocks.h,by="idg")
  
  
  # Function to test hyperparameters of the BART (ntree, k, power and base)
  # Records are the records of the species in vector class
  # ExplainVars are the predictors to use in data.frame format
  # Folds for cross-validation must be in BIOMOD2 format
  bartHO_Eq(records = db[,spp.h],explainVars = db[,predictors],folds = db[,c("RUN1", "RUN2",  "RUN3", "RUN4", "RUN5")],name = spp.h,
           ntree = c(100,150,200,250,300),k = c(1,2,3),power = c(1,2,3),base = c(0.75,0.85,0.95))
           
           }

# stopCluster(cl)


#predictive performance final models by species
for(i in 1:length(species)){
name <- species[i]
d.h <- read.csv(paste("hyperParamsEval/hyperEval",name,".csv",sep = ""))
d.h <- d.h[order(d.h$BOYCE,decreasing = T),][1,]
d.h$Species <- name
if(i==1){d.out <- d.h}else{d.out <- rbind(d.out,d.h)}}


mean(d.out$BOYCE)%>%round(2)
mean(d.out$AUC)%>%round(2)


boxplot(d.out[,c(1,3)])

sink("summary.txt")
paste("boyce mean")
summary(d.out$BOYCE)
paste("boyce sd")
sd(d.out$BOYCE)%>%round(2)

paste("AUC mean")
summary(d.out$AUC)
paste("AUC sd")
sd(d.out$AUC)%>%round(2)
sink()



table(d.out$ntree)
table(d.out$k)
table(d.out$base)
table(d.out$power)

d.out[,c(1,2,3,4)] <- d.out[,c(1,2,3,4)]%>%round(2)

write.csv(d.out,file = "summaySDM_PredictivePerformance.csv",row.names = F)
