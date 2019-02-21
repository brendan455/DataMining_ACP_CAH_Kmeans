source("http://www.bioconductor.org/biocLite.R")
setwd("Current Working Directory")
RequiredfunctionsDir<-"requiredFiles/"
source(paste0(RequiredfunctionsDir,"/kmeans_cah_acp.r"))

###########################
####### Data IMPORT #######
###########################

legumes<-read.table(paste(getwd(),"/Data/Legumes.csv",sep=""),header=T,row.names=1,sep=";")

###########################
## Vérification données ###
###########################

boxplot(legumes$Eau)
boxplot(legumes$ENG)
boxplot(legumes$PROT)
boxplot(legumes$GLUC)
boxplot(legumes$LIPI)
boxplot(legumes$SUC)
boxplot(legumes$FIBR)

###########################
##### Calcul moyennes #####
###########################

summary(legumes)

###########################
### Clacul de Variances ###
###########################

VarLegumes<-(var(legumes)*(nrow(legumes)-1))/nrow(legumes)
VarLegumes

###########################
### Coef de Correlation ###
###########################

corleg<-cor(legumes)
corleg

###########################
# Données centr. réduit. ##
###########################

CRlegumes<-centreduire(legumes)
CRlegumes

###########################
## Calcul distance Eucl. ##
###########################

dlegumes<-dist(CRlegumes)
dlegumes

###########################
########## ACP ############
###########################

ACP<-ACPN(legumes)
plot(ACP)
ACPP<-VP(ACP) #Part d'intertie + coeff de kaiser
coordonneespts<-ACP$score #Coordonnées points
matrice<-AXEVAR(legumes,ACP) #Matrice de corrélation
matrice

###########################
########## Cos² ###########
###########################

cos2<-COS2TETA(ACP,3) #Contribution des variablas à 3 axes sélectionnées
cos2

contrib<-CTR(ACP,3) #Contribution des individus aux axes
contrib

planacp<-PLAN(resacp = ACP, i = 1, j = 3)

###########################
######### Kmeans ##########
###########################

kmeans(legumes,7)
kleg<-kmeans(legumes,10)
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

cahward <- ward(legumes)
plot(cahward)
plot(cahward$height)

cah.min <- hclust(dlegumes, method = "single")
plot(cah.min, hang = -1)
partitions<-couper(cahward,5)
as.data.frame(partitions)
Max_classes<-getCount(partitions)
Max_classes[which(Max_classes$count %in% max(Max_classes$count)),]
rap_inertie(cahward,5)
cdgcl1(legumes,cahward,5) #Centre gravité
ctrcl(legumes,cahward,5) #Contributions relatives
ctrng(legumes,cahward,5) #Contribution CTR
rho2(legumes,cahward,5) #Carré des distances
iintra(legumes,cahward,5,3) #Inertie intraclasse
rect.hclust(cahward,k=5)
