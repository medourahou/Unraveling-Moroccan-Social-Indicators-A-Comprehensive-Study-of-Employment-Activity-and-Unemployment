#Premièrement on va importer les bibliothèques qu'on va utiliser
library("FactoMineR") # ici on a la fonction qui réalise l'ACP
library("xlsx") # elle permet de manipuler les fichier xlsx
library("writexl") #fcts permettant d'exporter des data frames en tant que fichier excel
library("factoextra")
library("corrplot")


#aller au dossier contenant les données
setwd("C:\\Users\\hp\\Desktop\\M1S1\\analyse_des_données\\projet")

#charger les données
dataset2 <- read.xlsx("dataset2.xlsx",sheetIndex = 1, header=TRUE, col_names=TRUE)

#Pour utiliser les années en tant que noms des individus au lieu de 
#1,2,...,7
rownames(dataset2) <- dataset2[,1]



resultACP <- PCA(dataset2[1:7, 2:22], ncp=21, scale.unit = TRUE, graph = TRUE)

resultACP$var
resultACP$ind

fviz_contrib(resultACP,choice = "ind", axes=1)
fviz_contrib(resultACP,choice = "ind", axes=2)
fviz_contrib(resultACP,choice = "ind", axes=3)

fviz_cos2(resultACP,choice = "ind", axe=1)
fviz_cos2(resultACP,choice = "ind", axe=2)
fviz_cos2(resultACP,choice = "ind", axe=3)

fviz_contrib(resultACP,choice = "var", axes=1)
fviz_contrib(resultACP,choice = "var", axes=2)
fviz_contrib(resultACP,choice = "var", axes=3)

fviz_cos2(resultACP,choice = "var", axe=1)
fviz_cos2(resultACP,choice = "var", axe=2)
fviz_cos2(resultACP,choice = "var", axe=3)

fviz_eig(resultACP, addlabels = TRUE, ylim = c(0, 60))

corrplot(resultACP$var$cos2, is.corr=FALSE)

corrplot(resultACP$var$contrib, is.corr=FALSE)    

fviz_pca_ind(resultACP, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
