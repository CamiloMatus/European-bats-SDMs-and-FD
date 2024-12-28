#This function is used to optimize the hyperparameters of the BART algorithm. 
#The same number of pseudo absences (random) and presences are used 
#as recommended for machine learning techniques based on decision trees. 
#Also, mean predictions are made based on 10 repetitions 
#(each repetition with random pseudo absences), 
#in order to achieve the best predictive performances as recommended.
#Both recommendations come from the study of:
#https://doi.org/10.1111/j.2041-210X.2011.00172.x


bartHO_Eq<- function(records,explainVars,folds,ntree=200,k=2,power=2,base=0.95,name=NULL){

if(sum(records[records==1])<500){nrep <- 10}else{nrep <- 1}

if(file.exists("randomPseudo-Absences")==FALSE){
  dir.create(file.path(getwd(), "randomPseudo-Absences"))}

if(file.exists("hyperParamsEval")==FALSE){
  dir.create(file.path(getwd(), "hyperParamsEval"))}
  
if(is.null(name)==T){name <- "speciesX"}  

#hyperparameters combinations to test  
HyptoTest <- tidyr::expand_grid(ntree, k,power,base) %>% data.frame()

#Defining random pseudo-absences - background - random points 
#As equal number of pseudo-absences and presences in each fold
dat.h00 <-  data.frame(records,explainVars,folds)
for (r in 1:nrep) {
  rp.out <- c()
  for (f in 1:ncol(folds)) {
    
  #number of presence records at f fold
  n.p <- nrow(dat.h00[dat.h00$records==1&dat.h00[,((ncol(explainVars)+1)+f)]==FALSE,])
  
  #random points of pseudo-absences at f fold (same number as presences records in that fold) 
  rp.h <- sample(rownames(dat.h00[dat.h00$records==0&dat.h00[,((ncol(explainVars)+1)+f)]==FALSE,]),
                  size = n.p)
  rp.out <- c(rp.out,rp.h)}
  
  assign(paste("randomPoint_Replication_",r,sep = ""),rp.out)
  write.csv(get(paste("randomPoint_Replication_",r,sep = "")),
  file = paste("randomPseudo-Absences/randomPseudoAbs_Replication_",r,"_",name,".csv",sep = ""))  
}


#testing hyperparameters 
for (h in 1:nrow(HyptoTest)) {
auc <- c()
boyce <- c()  
for (run in 1:ncol(folds)) {  
      message(paste("Testing hyperparameter in Fold", run))
      message(paste("ntrees =", HyptoTest[h,"ntree"],
                    "; k =",HyptoTest[h,"k"],"; base =",
                    HyptoTest[h,"base"],"; power =",
                    HyptoTest[h,"power"],
                    "...","combination",h,"of",nrow(HyptoTest),sep=" "))

dat.h0 <-  data.frame(records,explainVars,folds[,run])    

#This data (dat.h.RP.meanPred) of pseudo-absences is used for 
#mean predictions for model validation when replications are employed.
#The data contains all pseudo-absences from all replications
if(nrep!=1){
for (replication in 1:nrep) {
        d.h <- 
        dat.h0[row.names(dat.h0)%in%get(paste("randomPoint_Replication_",replication,sep = "")),]
        if(replication==1){dat.h.RP.meanPred <- d.h
        }else{dat.h.RP.meanPred <- rbind(dat.h.RP.meanPred,d.h)}}}
#####      
    

    for (replication in 1:nrep) {
    if(nrep!=1){
    message(paste("Replication",replication,"of",nrep))}else{NULL}

      dat.h <- rbind(dat.h0[dat.h0$records==1,],
                  dat.h0[row.names(dat.h0)%in%get(paste("randomPoint_Replication_",replication,sep = "")),])
      
      ####condition for make validation on fold 
      ####if fold do not have observations, NA value for boyce and AUC is assigned
      ####this occur when blockCV function fail in putting presences records in all folds
      if(nrow(dat.h[dat.h[,"folds...run."]==FALSE,])!=0|nrow(dat.h[dat.h[,"folds...run."]==FALSE,])==1){
      
      #model fitting using same number of presence records and pseudo-absences
      mod <- bart(x.train =  dat.h[dat.h[,"folds...run."]==TRUE,names(explainVars)], 
                  y.train=dat.h[dat.h[,"folds...run."]==TRUE,"records"],
                  k = paste(HyptoTest[h,"k"])%>%as.numeric(),
                  power = paste(HyptoTest[h,"power"])%>%as.numeric(),
                  base = paste(HyptoTest[h,"base"])%>%as.numeric(),
                  ntree = paste(HyptoTest[h,"ntree"])%>%as.numeric(),
                  nchain = 1, 
                  keeptrees  = TRUE, verbose  = F)
     
  
      #data for validation when replications are not used
      if(nrep==1){
      s <- list() 
      for (e in 1:length(explainVars)) {
        r.h<- raster(nrow=nrow(dat.h[dat.h[,"folds...run."]==FALSE,]),ncol=1)
        values(r.h) <- dat.h[dat.h[,"folds...run."]==FALSE,names(explainVars)[e]]
        s[[e]] <- r.h} 
      s <- stack(s)
      names(s) <-  names(explainVars)}else{
      
      #data for validation when replications are used  
      s <- list() 
        for (e in 1:length(explainVars)) {
          dat.h <- rbind(dat.h0[dat.h0$records==1,],
                         dat.h.RP.meanPred)
          r.h<- raster(nrow=nrow(dat.h[dat.h[,"folds...run."]==FALSE,]),ncol=1)
          values(r.h) <- dat.h[dat.h[,"folds...run."]==FALSE,names(explainVars)[e]]
          s[[e]] <- r.h} 
        s <- stack(s)
        names(s) <-  names(explainVars)}
      
        
      #predictions for validation
      pred <- predict2.bart(object = mod,x.layers =s)%>% values() %>% as.numeric()
      if(replication==1){pred.h <- pred}else{pred.h <- cbind(pred.h,pred)}
      }else{NULL}
      
      }
      
      ####condition for make validation on fold 
      ####if fold do not have observations, NA value for Boyce and AUC is assign 
      if(nrow(dat.h[dat.h[,"folds...run."]==FALSE,])!=0|nrow(dat.h[dat.h[,"folds...run."]==FALSE,])==2){
      if(nrep==1){pred.out <- pred.h}else{
      #mean predictions of replicated psudo-absences
      pred.out <- rowMeans(pred.h,na.rm = T)}

      auc.h <- auc(DATA = data.frame(1:nrow(dat.h[dat.h[,"folds...run."]==FALSE,names(explainVars)]),
                                     dat.h[dat.h[,"folds...run."]==FALSE,"records"],pred.out))[1]%>%as.numeric()
      dat.h3 <- data.frame(dat.h[dat.h[,"folds...run."]==FALSE,"records"],pred.out)
      boyce.h <- ecospat.boyce(fit = pred.out,obs = dat.h3[dat.h3[,1]==1,2],PEplot = F)$cor %>% as.numeric()
      
      }else{auc.h <- NA;boyce.h <- NA}
      
      auc <- c(auc,auc.h)
      boyce <- c(boyce,boyce.h)
      if(nrep!=1){
      message(paste("Results of mean predictions"))}
      message(paste("auc = ",round(auc.h,2)))
      message(paste("boyce = ",round(boyce.h,2)))
              }   
    
    if(h==1){
      dat.out <- data.frame(auc,boyce)
      dat.out$ntree <- HyptoTest[h,"ntree"]
      dat.out$k <- HyptoTest[h,"k"]
      dat.out$power <- HyptoTest[h,"power"]
      dat.out$base <- HyptoTest[h,"base"]
      dat.out$species <- name
      dat.out$Fold <- 1:ncol(folds)
      dat.out.f <- dat.out  
      }else{
        dat.out <- data.frame(auc,boyce)
        dat.out$ntree <- HyptoTest[h,"ntree"]
        dat.out$k <- HyptoTest[h,"k"]
        dat.out$power <- HyptoTest[h,"power"]
        dat.out$base <- HyptoTest[h,"base"]
        dat.out$species <-name
        dat.out$Fold <- 1:ncol(folds)
        dat.out.f <- rbind(dat.out.f,dat.out)  
      }
}
  dat.out.f2 <- dat.out.f
  dat.out.f2$hp <- paste(dat.out.f2$ntree,dat.out.f2$k,dat.out.f2$power,dat.out.f2$base)
  
  BOYCE <- tapply(dat.out.f2$boyce, dat.out.f2$hp, mean,na.rm = T)
  BOYCE.SD <- tapply(dat.out.f2$boyce, dat.out.f2$hp, sd,na.rm = T)
  AUC <- tapply(dat.out.f2$auc, dat.out.f2$hp, mean,na.rm = T)
  AUC.SD <- tapply(dat.out.f2$auc, dat.out.f2$hp, sd,na.rm = T)
  
  out <- data.frame(BOYCE,BOYCE.SD,AUC,AUC.SD)
  out$hp <- rownames(out)
  
  out2 <- merge(x=out,y=dat.out.f2[,c(3:6,9)],by="hp",all.y=F)
  out3 <- distinct(out2, .keep_all = TRUE)
  
  write.csv(out3[,2:9],file = paste("hyperParamsEval/hyperEval",name,".csv",sep = ""),row.names = F)
  write.csv(out3[,2:9],file = paste("hyperParamsEval/hyperEval",name,".csv",sep = ""),row.names = F)
}