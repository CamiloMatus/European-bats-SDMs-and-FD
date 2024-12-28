require(embarcadero)
require(PresenceAbsence)
require(sf)
require(dplyr)
require(ecospat)
require(foreign)

#reading data
dbPres <- read.csv("../1dataBase/speciesPresences.csv")
predictors <- read.dbf("../1dataBase/bioHistoric_SampRate.dbf")
modsEval <- read.csv("../3hiperparametersOptimization/summaySDM_PredictivePerformance.csv")


climateProyections <- list.files("../1dataBase/")[str_detect(string = list.files("../1dataBase/"),pattern = "ssp")]
climateProyections <- climateProyections[str_detect(string = climateProyections,pattern = "dbf")]
climateProyections <- str_remove(climateProyections,".dbf")

clpr <- c()#climate proyection
for (s in climateProyections) {
  clpr.h <- read.dbf(paste("../1dataBase/",s,".dbf",sep = ""))
  clpr.h$clpr <- s
  clpr <- rbind(clpr,clpr.h)
}

predictors$clpr <-  "historic"
datAll <- rbind(predictors[,c(1:20,22)],clpr[,1:21])
namesEV <- c("bio1","bio2","bio7","bio8","bio12","bio15","bio18","sampRate")
explainVars <- predictors[,namesEV]
explainVars$idg <- predictors$idg

dat<- merge(dbPres,explainVars,by="idg")
species <- names(dbPres)[1:38]


#doing predictions by species
for(i in species){
  message(paste("working on ",i,"   ",which(species==i ),"of",length(species)))
  
  #hyperparameters that reached the best predictive performance
  hyperParams.h <- modsEval[modsEval$Species==i,]

  #condition for species with few records
  if(sum(dat[,i])<500){nrep <- 8}else{nrep <- 1}
  
  for (replication in 1:nrep) {
    if(nrep!=1){
      message(paste("Fitting for replication",replication,"of",nrep))}else{NULL}
    
    #random pseudo-absences used in model evaluation
    rp <- read.csv(paste("../3hiperparametersOptimization/randomPseudo-Absences/randomPseudoAbs_Replication_",replication,"_",i,".csv",sep = ""))
    
    #data with pseudo-absences used in model evaluation
    dat.1 <- dat[row.names(dat)%in%rp$x,]
    
    #data with only presences records
    dat.2 <- dat[dat[,i]==1,]
    
    #data for model fitting
    dat.h <- rbind(dat.1,dat.2)
    
    #fitting model
    message("fitting model")
    mod <- bart2(formula =  dat.h[,namesEV],
                 k = paste(hyperParams.h[1,"k"])%>%as.numeric(),
                 power = paste(hyperParams.h[1,"power"])%>%as.numeric(),
                 base = paste(hyperParams.h[1,"base"])%>%as.numeric(),
                 n.tree = paste(hyperParams.h[1,"ntree"])%>%as.numeric(),
                 n.chains = 1,
                 data =  dat.h[,i], 
                 keepTrees = TRUE, verbose = FALSE)

    #creating data for projections
    
    #------------------------------------------------------------------------------------
    # This would be the simple way to make the predictions, but it requires a lot of RAM, 
    # for this reason I have split the data for the projections into 12 parts.
    #------------------------------------------------------------------------------------
    
    # s <- list() 
    # for (e in 1:length(namesEV)) {
    #   
    #   if(namesEV[e]!="sampRate"){
    #     r.h<- raster(nrow=nrow(datAll),ncol=1)
    #     values(r.h) <- datAll[,namesEV[e]]
    #     s[[e]] <- r.h
    #   }else{
    #     r.h<- raster(nrow=nrow(datAll),ncol=1)
    #     values(r.h) <- rep(0,nrow(datAll))## sampRate equal to 0 for sampling bias mitigation
    #     s[[e]] <- r.h } } 
    # 
    # #data for proyections
    # s <- stack(s)
    # names(s) <-  namesEV
    # 
    # #doing projections
    # message("doing projections 1")
    # pred1 <- predict2.bart(object = mod,x.layers =s,quiet = F)%>% values() %>% as.numeric()
    # gc()
    # 
    # 
    # pred <- c(pred1)
    # if(replication==1){pred.h <- pred}else{pred.h <- cbind(pred.h,pred)}

    
    s <- list()
    for (e in 1:length(namesEV)) {

      if(namesEV[e]!="sampRate"){
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- datAll[1:165529,namesEV[e]]
        s[[e]] <- r.h
      }else{
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- rep(0,165529)## sampRate equal to 0 for sampling bias mitigation
        s[[e]] <- r.h } }

    s2 <- list()
    for (e in 1:length(namesEV)) {

      if(namesEV[e]!="sampRate"){
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- datAll[165530:331058,namesEV[e]]
        s2[[e]] <- r.h
      }else{
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- rep(0,165529)## sampRate equal to 0 for sampling bias mitigation
        s2[[e]] <- r.h } }

    s3 <- list()
    for (e in 1:length(namesEV)) {

      if(namesEV[e]!="sampRate"){
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- datAll[331059:496587,namesEV[e]]
        s3[[e]] <- r.h
      }else{
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- rep(0,165529)## sampRate equal to 0 for sampling bias mitigation
        s3[[e]] <- r.h } }

    s4 <- list()
    for (e in 1:length(namesEV)) {

      if(namesEV[e]!="sampRate"){
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- datAll[496588:662116,namesEV[e]]
        s4[[e]] <- r.h
      }else{
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- rep(0,165529)## sampRate equal to 0 for sampling bias mitigation
        s4[[e]] <- r.h } }


    s5 <- list()
    for (e in 1:length(namesEV)) {
      
      if(namesEV[e]!="sampRate"){
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- datAll[662117:827645,namesEV[e]]
        s5[[e]] <- r.h
      }else{
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- rep(0,165529)## sampRate equal to 0 for sampling bias mitigation
        s5[[e]] <- r.h } }
  
    
    s6 <- list()
    for (e in 1:length(namesEV)) {
      
      if(namesEV[e]!="sampRate"){
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- datAll[827646:993174,namesEV[e]]
        s6[[e]] <- r.h
      }else{
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- rep(0,165529)## sampRate equal to 0 for sampling bias mitigation
        s6[[e]] <- r.h } }
    
   
    
    s7 <- list()
    for (e in 1:length(namesEV)) {
      
      if(namesEV[e]!="sampRate"){
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- datAll[993175:1158703,namesEV[e]]
        s7[[e]] <- r.h
      }else{
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- rep(0,165529)## sampRate equal to 0 for sampling bias mitigation
        s7[[e]] <- r.h } }
    
    
    
    s8 <- list()
    for (e in 1:length(namesEV)) {
      
      if(namesEV[e]!="sampRate"){
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- datAll[1158704:1324232,namesEV[e]]
        s8[[e]] <- r.h
      }else{
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- rep(0,165529)## sampRate equal to 0 for sampling bias mitigation
        s8[[e]] <- r.h } }
  
    
    s9 <- list()
    for (e in 1:length(namesEV)) {
      
      if(namesEV[e]!="sampRate"){
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- datAll[1324233:1489761,namesEV[e]]
        s9[[e]] <- r.h
      }else{
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- rep(0,165529)## sampRate equal to 0 for sampling bias mitigation
        s9[[e]] <- r.h } }
    
    
    
    s10 <- list()
    for (e in 1:length(namesEV)) {
      
      if(namesEV[e]!="sampRate"){
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- datAll[1489762:1655290,namesEV[e]]
        s10[[e]] <- r.h
      }else{
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- rep(0,165529)## sampRate equal to 0 for sampling bias mitigation
        s10[[e]] <- r.h } }
    
    
    
    s11 <- list()
    for (e in 1:length(namesEV)) {
      
      if(namesEV[e]!="sampRate"){
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- datAll[1655291:1820819,namesEV[e]]
        s11[[e]] <- r.h
      }else{
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- rep(0,165529)## sampRate equal to 0 for sampling bias mitigation
        s11[[e]] <- r.h } }
    
    
    
    s12 <- list()
    for (e in 1:length(namesEV)) {
      
      if(namesEV[e]!="sampRate"){
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- datAll[1820820:1986348,namesEV[e]]
        s12[[e]] <- r.h
      }else{
        r.h<- raster(nrow=165529,ncol=1)
        values(r.h) <- rep(0,165529)## sampRate equal to 0 for sampling bias mitigation
        s12[[e]] <- r.h } }
    
    
    #data for proyections
    s <- stack(s)
    names(s) <-  namesEV

    s2 <- stack(s2)
    names(s2) <-  namesEV
    
    s3 <- stack(s3)
    names(s3) <-  namesEV
    
    s4 <- stack(s4)
    names(s4) <-  namesEV
    
    s5 <- stack(s5)
    names(s5) <-  namesEV
    
    s6 <- stack(s6)
    names(s6) <-  namesEV
    
    s7 <- stack(s7)
    names(s7) <-  namesEV
    
    s8 <- stack(s8)
    names(s8) <-  namesEV
    
    s9 <- stack(s9)
    names(s9) <-  namesEV
    
    s10 <- stack(s10)
    names(s10) <-  namesEV
    
    s11 <- stack(s11)
    names(s11) <-  namesEV
    
    s12 <- stack(s12)
    names(s12) <-  namesEV
    
    
    #doing projections
    message("doing projections 1")
    pred1 <- predict2.bart(object = mod,x.layers =s,quiet = F)%>% values() %>% as.numeric()
    gc()
    
    message("doing projections 2")
    pred2 <- predict2.bart(object = mod,x.layers =s2,quiet = F)%>% values() %>% as.numeric()
    gc()
    
    message("doing projections 3")
    pred3 <- predict2.bart(object = mod,x.layers =s3,quiet = F)%>% values() %>% as.numeric()
    gc()
    
    message("doing projections 4")
    pred4 <- predict2.bart(object = mod,x.layers =s4,quiet = F)%>% values() %>% as.numeric()
    gc()
    
    
    message("doing projections 5")
    pred5 <- predict2.bart(object = mod,x.layers =s5,quiet = F)%>% values() %>% as.numeric()
    gc()
    
    message("doing projections 6")
    pred6 <- predict2.bart(object = mod,x.layers =s6,quiet = F)%>% values() %>% as.numeric()
    gc()
    
    message("doing projections 7")
    pred7 <- predict2.bart(object = mod,x.layers =s7,quiet = F)%>% values() %>% as.numeric()
    gc()
    
    message("doing projections 8")
    pred8 <- predict2.bart(object = mod,x.layers =s8,quiet = F)%>% values() %>% as.numeric()
    gc()
    
    message("doing projections 9")
    pred9 <- predict2.bart(object = mod,x.layers =s9,quiet = F)%>% values() %>% as.numeric()
    gc()
    
    message("doing projections 10")
    pred10 <- predict2.bart(object = mod,x.layers =s10,quiet = F)%>% values() %>% as.numeric()
    gc()
    
    message("doing projections 11")
    pred11 <- predict2.bart(object = mod,x.layers =s11,quiet = F)%>% values() %>% as.numeric()
    gc()
    
    message("doing projections 12")
    pred12 <- predict2.bart(object = mod,x.layers =s12,quiet = F)%>% values() %>% as.numeric()
    gc()
    
    pred <- c(pred1,pred2,pred3,pred4,pred5,pred6,pred7,pred8,pred9,pred10,pred11,pred12)
    
    
    if(replication==1){pred.h <- pred}else{pred.h <- cbind(pred.h,pred)}
  
  }
    
    if(nrep==1){pred.mean <- pred.h}else{
    #mean predictions of replicated psudo-absences
    pred.mean <- rowMeans(pred.h,na.rm = T)}

      dat.h2 <- data.frame(datAll$idg[1:nrow(datAll[datAll$clpr=="historic",])],pred.mean[1:nrow(datAll[datAll$clpr=="historic",])])
      
      names(dat.h2) <- c("idg","pred")
      
      dat.h3 <- merge(dbPres,dat.h2,by="idg")
      
      d.h0 <- data.frame(1:nrow(dat.h3),dat.h3[,i],dat.h3$pred)
      
      message("computing optimal threshold")
      th <- optimal.thresholds(d.h0,opt.methods = 3,threshold = 1000)[2]%>%as.numeric()
      d.h <- data.frame(datAll$idg,pred.mean)
      d.h$Bin <- 0
      names(d.h)[1] <- "idg"
      names(d.h)[2] <- paste("Prob",i,sep = ".")
      names(d.h)[3] <- paste("Bin",i,sep = ".")
      d.h[d.h[,2]>=th,3] <- 1
      d.h$clpr <- datAll$clpr
      message("saving results")
      write.csv(d.h,file = paste("predictions",i,"csv",sep = "."),row.names = F)
}


# Creating final data for suitable area predictions
for (i in 1:length(species)) {
  message(paste("working on species",i,"of",length(species)))
  if(i==1){
    d.h <- read.csv(paste("predictions",species[i],"csv",sep = "."))[,c(1,4,3)]
    gc()
    names(d.h)[3] <-   species[i]
  }else{d.h <-  cbind(d.h,read.csv(paste("predictions",species[i],"csv",sep = "."))[,c(3)])
   gc()
   names(d.h)[ncol(d.h)] <-  species[i]}}

write.csv(d.h,file = "speciesPredictionsDat.csv",row.names = F)
