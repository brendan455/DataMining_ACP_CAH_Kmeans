source("http://www.bioconductor.org/biocLite.R")
setwd("C:/Users/brdeniau/Documents/Projet R")
RequiredfunctionsDir<-"requiredFiles/"
source(paste0(RequiredfunctionsDir,"/kmeans_cah_acp.r"))

###########################
####### Data IMPORT #######
###########################

legumineuses<-read.table(paste(getwd(),"/Data/Legumes.csv",sep=""),header=T,row.names=1,sep=";")

###########################
## Vérification données ###
###########################

boxplot(legumineuses$Eau)
boxplot(legumineuses$ENG)
boxplot(legumineuses$PROT)
boxplot(legumineuses$GLUC)
boxplot(legumineuses$LIPI)
boxplot(legumineuses$SUC)
boxplot(legumineuses$FIBR)

###########################
##### Calcul moyennes #####
###########################

summary(legumineuses)

###########################
### Clacul de Variances ###
###########################

VarLegumes<-(var(legumineuses)*(nrow(legumineuses)-1))/nrow(legumineuses)
VarLegumes

###########################
### Coef de Correlation ###
###########################

corleg<-cor(legumineuses)
corleg

###########################
# Données centr. réduit. ##
###########################

CRlegumineuses<-centreduire(legumineuses)
CRlegumineuses

###########################
## Calcul distance Eucl. ##
###########################

dlegumes<-dist(CRlegumineuses)
dlegumes

###########################
########## ACP ############
###########################

ACP<-ACPN(legumineuses)
plot(ACP)
ACPP<-VP(ACP) # part d'intertie + coeff de kaiser
coordonneespts<-ACP$score # coordonnées points
matrice<-AXEVAR(legumineuses,ACP) #Matrice de corrélation
matrice

###########################
########## Cos² ###########
###########################

cos2<-COS2TETA(ACP,3) #contrib des variablas à 3 axes sélectionnées
cos2

contrib<-CTR(ACP,3)#Contribution des individus aux axes
contrib

planacp<-PLAN(resacp = ACP, i = 1, j = 3)

###########################
######### Kmeans ##########
###########################

kmeans(legumineuses,7)
kleg<-kmeans(legumineuses,10)
kleg

kgrp<-kclass(kleg)
as.data.frame(kgrp)
getCount <- function(x) {
  u <- unique(x);
  data.frame(
    value=u,
    count=sapply(u, function(v) { length(which(x==v)) } )
  )
}

###########################
########### CAH ###########
###########################

cahward <- ward(legumineuses)
plot(cahward)
plot(cahward$height)

cah.min <- hclust(dlegumes, method = "single")
plot(cah.min, hang = -1)
partitions<-couper(cahward,5)
as.data.frame(partitions)
Max_classes<-getCount(partitions)
Max_classes[which(Max_classes$count %in% max(Max_classes$count)),]
rap_inertie(cahward,5)
cdgcl1(legumineuses,cahward,5) # Centre gravité
ctrcl(legumineuses,cahward,5) # Contributions relatives
ctrng(legumineuses,cahward,5) # Contribution CTR
rho2(legumineuses,cahward,5) # Carré des distances
iintra(legumineuses,cahward,5,3) # Inertie intraclasse
rect.hclust(cahward,k=5)

