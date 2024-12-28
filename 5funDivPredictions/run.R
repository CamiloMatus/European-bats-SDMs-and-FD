require(ade4)
require(sf)
require(leaflet)
require(dplyr)
require(stringr)
library(dplyr)
require(foreach)
require(doParallel)
library("magrittr")

db <- read.csv("../4SDM_predictions/speciesPredictionsDat.csv")
dbPres <- read.csv("../1dataBase/speciesPresences.csv")
species <- names(db)[3:40]

dt2 <- read.csv("../1dataBase/funcTraitEuroBats.csv",stringsAsFactors=TRUE)


# Preparing functional trait data
str(dt2)
names(dt2)
dt2[,names(dt2)[14:23]] <- lapply(dt2[,names(dt2)[14:23]],as.numeric)#Must be numeric to be treated as dichotomous variables
dt2[,names(dt2)[35:40]] <- lapply(dt2[,names(dt2)[35:40]],as.character)#They are actually factors
dt2[,names(dt2)[35:40]] <- lapply(dt2[,names(dt2)[35:40]],as.factor)

##
levels(dt2$roos.OvergroundRoostDependence.Castle) <- c("1","2","3","4","5")
levels(dt2$roos.OvergroundRoostDependence.Church) <- c("1","2","3","4","5")
levels(dt2$roos.OvergroundRoostDependence.House) <- c("1","2","3","4","5")
levels(dt2$roos.OvergroundRoostDependence.Barn) <- c("1","2","3","4","5")
levels(dt2$roos.OvergroundRoostDependence.Bridge) <- c("1","2","3","4","5")
levels(dt2$roos.OvergroundRoostDependence.Tree) <- c("1","2","3","4","5")
str(dt2)


trait_type <- lapply(dt2,class)[2:ncol(dt2)]%>%data.frame()%>%t()
trait_name <- names(dt2)[2:ncol(dt2)]
fuzzy_name <- substr(trait_name,1,4)

# Weight of traits
dw <- data.frame((1/length(unique(fuzzy_name)))/table(fuzzy_name))

dt3 <- data.frame(trait_name,trait_type,fuzzy_name)
row.names(dt3)<- 1:nrow(dt3)
dt4 <- merge(dt3,dw,by="fuzzy_name")
dt5 <- dt4[,2:4]
names(dt5)[3] <- "trait_weight"

dt5[dt5$trait_type=="numeric","trait_type"] <- "Q"
dt5[dt5$trait_type=="factor","trait_type"] <- "N"
dt5[dt5$trait_type=="integer","trait_type"] <- "C"


# Compute functional distance
row.names(dt2) <- dt2$species
distTrait <- mFD::funct.dist(sp_tr = dt2[2:ncol(dt2)],
                             tr_cat = dt5,
                             metric = "gower",
                             # scale_euclid = "scale_center",
                             # ordinal_var = "classic",
                             # weight_type = "equal",
                             stop_if_NA = F)

# trait distance
MatDis <- as.matrix(distTrait)[1:38, 1:38]

# Retrieving principal coordinates
pco <- dudi.pco(distTrait, scann = FALSE)
summary(pco)
inertia.dudi(pco)
namesPCO <- row.names(pco$tab)
sp_faxes <- as.matrix(pco$tab)
row.names(sp_faxes) <- row.names(pco$tab)




#making DF estimates in each geographic grid for each time period and GEIs emissions scenario

registerDoParallel(cl <- makeCluster(5))
db$sp_richn <- NA
db$fdis     <- NA
db$fmpd      <- NA
db$fnnd     <- NA
db$feve <- NA
db$fric      <- NA
db$fdiv      <- NA
db$fori      <- NA
db$fspe<- NA



# Estimating functional diversity in each geographic grid of Europe (5 minutes) 
# according to several metrics
out<- foreach(i = 1:nrow(db),.packages = "dplyr",.combine = "rbind" ) %dopar% {

  a <- rep(0,38)
  a <- data.frame(a)
  row.names(a) <- row.names(pco$tab)
  names.h <- names(colSums(db[i,3:40])[colSums(db[i,3:40])>0])
  a[rownames(a)%in%names.h,"a"] <- 1
  
  
  if(length(names.h)>5){
    alpha_fd_indices <- mFD::alpha.fd.multidim(
      sp_faxes_coord = sp_faxes[,1:5],
      asb_sp_w = t(a),
      ind_vect = c("fdis", "fmpd", "fnnd", "feve", "fric", "fdiv",
                   "fori", "fspe"),
      scaling = TRUE,
      check_input = TRUE,
      details_returned = TRUE,
      verbose=F)

    
    db[i,"sp_richn"] <- length(names.h)
    db[i,"fdis"]<- alpha_fd_indices$functional_diversity_indices$fdis
    db[i,"fmpd"]<- alpha_fd_indices$functional_diversity_indices$fmpd
    db[i,"fnnd"]<- alpha_fd_indices$functional_diversity_indices$fnnd
    db[i,"feve"]<- alpha_fd_indices$functional_diversity_indices$feve
    db[i,"fric"]<- alpha_fd_indices$functional_diversity_indices$fric
    db[i,"fdiv"]<- alpha_fd_indices$functional_diversity_indices$fdiv
    db[i,"fori"]<- alpha_fd_indices$functional_diversity_indices$fori
    db[i,"fspe"]<- alpha_fd_indices$functional_diversity_indices$fspe
  }else{  db[i,"sp_richn"] <- length(names.h)}
  
  db[i,]
  

}

write.csv(out,file = "FuncDivPredictionsDat.csv",row.names = F)


