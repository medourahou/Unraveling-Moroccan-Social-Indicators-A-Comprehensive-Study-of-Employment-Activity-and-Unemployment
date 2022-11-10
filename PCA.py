import math
import numpy as np
from numpy.linalg import eig
import time

def ACP(matrice):
    
    fichierResultat = open("./resultatACP.txt","w")
    
    # cacule des moyennes et des écart-types de chaque colones ( de chaque variable )
    nbVars = len(matrice[0]) -1 # car la première colonne contient les noms des individus
    nbInds = len(matrice) - 1 # car la première ligne contient les noms des variables
    
    nomVars = matrice[0][:]
    
    nomInds=[]
    for i in range(1,nbInds+1):
        nomInds.append(matrice[i][0])
    
    
    
    moyennes = []
    ecarts_type=[]
    
    for j in range(1,nbVars+1):
        moyennes.append(0)
        for i in range(1,nbInds+1):
            moyennes[j-1]+=matrice[i][j]
        moyennes[j-1]/=nbInds
        
    
    
    for j in range(1,nbVars+1):
        ecarts_type.append(0)
        for i in range(1,nbInds+1):
            ecarts_type[j-1]+=(matrice[i][j]-moyennes[j-1])**2
        ecarts_type[j-1]/=nbInds
        ecarts_type[j-1]=math.sqrt(ecarts_type[j-1])
    
    """
    CentreDeGraviteG = moyennes[:]
    InertieTotaleI = 0
    
    for ecart_type in ecarts_type:
        InertieTotaleI += ecart_type**2
    """    
    
    
    
    #Matrice de corrélation
    
    
    
    
    
    
    matriceDeCorrelation=[]

    
    for j in range(0,nbVars):
        matriceDeCorrelation.append(list())
        for k in range(0,nbVars):
            matriceDeCorrelation[j].append(0)
            for i in range(0,nbInds):            
                matriceDeCorrelation[j][k]+= (matrice[i+1][j+1]-moyennes[j]) * (matrice[i+1][k+1]-moyennes[k])
            
            
            matriceDeCorrelation[j][k] /= nbInds
            matriceDeCorrelation[j][k] /= ecarts_type[j]*ecarts_type[k]
            
            if k==j :
                matriceDeCorrelation[j][k] = 1.0
                
    #Affichage de la matrice de corrélation
    
    print("----------------------------------------------------------------------")
    fichierResultat.write("----------------------------------------------------------------------\n")
    print("Matrice de corrélation des variables")
    fichierResultat.write("Matrice de corrélation des variables\n")
    
    tmpNomVars = ""
    for i in range(nbVars+1):
        if i == 0:
            tmpNomVars += " -------- "
        else:
            tmpNomVars += " " + nomVars[i] + " -------- "
    
    print(tmpNomVars)
    fichierResultat.write(tmpNomVars+"\n")
    
    for i in range(nbVars):
        print(nomVars[i+1]," ",matriceDeCorrelation[i][:])
        tmpp = nomVars[i+1]+" "+str(matriceDeCorrelation[i][:])+"\n"
        fichierResultat.write(tmpp)
    
    print("----------------------------------------------------------------------")
    fichierResultat.write("----------------------------------------------------------------------\n")
    
    
    # centrer et réduire
    
    matriceCentreReduite = matrice[:]
    #moyennesCentreReduite = [0]*nbVars
    #ecarts_typeCentreReduite = [1]*nbVars
    
    InertieTotaleCentreReduite = nbVars
    
    for j in range(1,nbVars+1):
        for i in range(1,nbInds+1):
            tmp = matriceCentreReduite[i][j]-moyennes[j-1]
            tmp /= ecarts_type[j-1]
            matriceCentreReduite[i][j]= tmp
    
    
    
    matriceDeCorrNumpy = [None]*nbVars
    for i in range(0,nbVars):
        matriceDeCorrNumpy[i]=np.array(matriceDeCorrelation[i][:])
        #print(i,"---",matriceCentreReduiteNumpy[i-1])
        
    matriceDeCorrNumpy = np.array(matriceDeCorrNumpy)
    #print("#########################\nmatriceCentreReduiteNumpy",matriceCentreReduiteNumpy)
    valspropres, vectsPropres = eig(matriceDeCorrNumpy)
    """ Attention: pour une même valeur propre est associé plusieurs vecteurs propres, donc les axes peuvent différer car ils dépendent du vecteur propre choisi"""
    #transformer en liste python
    
    valspropres = list(valspropres)
    
    valspropres = [ float(v) for v in valspropres ]
    
    vectsPropres = [list(vect) for vect in vectsPropres]
    
    for vect in vectsPropres:
        vect = [float(v) for v in vect ]
    
    #print("valspropres, ",valspropres,"\n vectsPropres ", vectsPropres)
    
    
    VectsETValsPropres = list()
    
    nbValsPrps = len(valspropres)
    
    for i in range(nbValsPrps):
        VectsETValsPropres.append( ( float(valspropres[i]), [ float(v) for v in vectsPropres[i]] ) )
        
    #sort vectsEtValsPropres decending order
    
    VectsETValsPropres = [elm for elm in sorted(VectsETValsPropres, key=lambda elm: elm[0], reverse=True)]
    
    print("----------------------------------------------------------------------")
    fichierResultat.write("----------------------------------------------------------------------\n")
    print("Valeurs et vecteurs propres")
    fichierResultat.write("Valeurs et vecteurs propres\n")
    
    for i in range(nbValsPrps):
        print("Valeur prore",i," : ",  VectsETValsPropres[i][0] )
        tmpp= "Valeur prore"+ str(i) + " : " +  str(VectsETValsPropres[i][0])+"\n"
        fichierResultat.write(tmpp)
        print("Vecteur prore",i," : ", VectsETValsPropres[i][1] )
        tmpp= "Vecteur prore"+ str(i) + " : " +  str(VectsETValsPropres[i][1])+"\n"
        fichierResultat.write(tmpp)
    
    print("----------------------------------------------------------------------")
    fichierResultat.write("----------------------------------------------------------------------\n")
    
    # cherchons les axes
    
    
    MatriceDesAxes =[] # [[0]*nbValsPrps]*nbInds
    
    for i in range(nbInds):
        MatriceDesAxes.append(list())
        for j in range(nbValsPrps):
            tmp=0
            for k in range(nbValsPrps):
                #if( i == 0 and j==0 ):
                    #print("matriceCentreReduite[i+1][k+1] ",matriceCentreReduite[i+1][k+1]," * ",VectsETValsPropres[j][1][k])
                tmp += matriceCentreReduite[i+1][k+1] * VectsETValsPropres[j][1][k]
            MatriceDesAxes[i].append(tmp)
    
    
    #affichage des coordonnées des individus
    
    nomAxes = [ ("Axe "+str(i+1)+" ") for i in range(nbValsPrps)]
    print("----------------------------------------------------------------------")
    fichierResultat.write("----------------------------------------------------------------------\n")
    print("Coordonnées des individus")
    fichierResultat.write("Coordonnées des individus\n")
    
    tmpNomsAxes = "--------"
    for i in range(nbVars):
        tmpNomsAxes += " " + nomAxes[i] + " -------- "
    
    print(tmpNomsAxes)
    fichierResultat.write(str(tmpNomsAxes)+"\n")
    for i in range(nbInds):
        print(nomInds[i]," ", MatriceDesAxes[i][:])
        tmpp= nomInds[i]+" "+ str(MatriceDesAxes[i][:])+"\n"
        fichierResultat.write(tmpp)
    print("----------------------------------------------------------------------")
    fichierResultat.write("----------------------------------------------------------------------\n")
    
    #% d'inertie pour les valeurs propres
    
    TableuPourcentage=[None]*3
    
    TableuPourcentage[0]=[ val[0] for val in  VectsETValsPropres]
    
    TableuPourcentage[1]=[ val[0]/InertieTotaleCentreReduite for val in  VectsETValsPropres]
    


    

    
    TableuPourcentage[2]=[  ]
    
    for i in range(0,nbValsPrps):
        if i==0:
            TableuPourcentage[2].append(TableuPourcentage[1][0])
        else:
            TableuPourcentage[2].append(TableuPourcentage[2][i-1]+TableuPourcentage[1][i])
    
    

    # Affichage des % d'inerties cumulé
    print("----------------------------------------------------------------------")
    fichierResultat.write("----------------------------------------------------------------------\n")
    print("Pourcentage d'inertie des valeurs propres et pourcentage d'inertie cumulé")
    fichierResultat.write("Pourcentage d'inertie des valeurs propres et pourcentage d'inertie cumulé\n")
    
    print("Valeurs Propres: ", TableuPourcentage[0][:] )
    fichierResultat.write("Valeurs Propres: "+str(TableuPourcentage[0][:])+"\n" )
    
    print("% inertie: ",TableuPourcentage[1][:])
    fichierResultat.write("% inertie: "+str(TableuPourcentage[1][:])+"\n")
    
    print("% inertie cumulé: ",TableuPourcentage[2][:])
    fichierResultat.write("% inertie cumulé: "+str(TableuPourcentage[2][:])+"\n")
    
    print("----------------------------------------------------------------------")
    fichierResultat.write("----------------------------------------------------------------------\n")
    # Contrib Individus
    
    SommeCarreCoord = []
    
    for j in range(nbVars):
        
        tmp = 0
        for i in range(nbInds):
            tmp+= MatriceDesAxes[i][j]**2
        SommeCarreCoord.append(tmp)
        
    #print(SommeCarreCoord)
    
    ContribIndvs = []
    
    for i in range(nbInds):
        ContribIndvs.append(list())
        for j in range(nbVars):
            ContribIndvs[i].append(((MatriceDesAxes[i][j]**2)/SommeCarreCoord[j])*100)
    
    
    #affichage des contributions des individus
    print("----------------------------------------------------------------------")
    fichierResultat.write("----------------------------------------------------------------------\n")
    print("contributions des individus")
    fichierResultat.write("contributions des individus\n")
    print(tmpNomsAxes)
    fichierResultat.write(str(tmpNomsAxes)+"\n")
    for i in range(nbInds):
        print(nomInds[i]," ", ContribIndvs[i][:])
        fichierResultat.write(nomInds[i]+" "+str(ContribIndvs[i][:])+"\n")
        
    
    print("----------------------------------------------------------------------")
    fichierResultat.write("----------------------------------------------------------------------\n")
    

    
    # Cos2 des individus
    
    SommeCarreCos2 = []
    
    for i in range(nbInds):
        
        tmp = 0
        for j in range(nbVars):
            tmp+= MatriceDesAxes[i][j]**2
        SommeCarreCos2.append(tmp)
        
    
    Cos2Indvs = []
    
    for i in range(nbInds):
        Cos2Indvs.append(list())
        for j in range(nbVars):
            Cos2Indvs[i].append(((MatriceDesAxes[i][j]**2)/SommeCarreCos2[i])*100)
    
    #affichage des qualités de représentation
    print("----------------------------------------------------------------------")
    fichierResultat.write("----------------------------------------------------------------------\n")
    print("Qualités de représentation des individus")
    fichierResultat.write("Qualités de représentation des individus\n")
    print(tmpNomsAxes)
    for i in range(nbInds):
        print(nomInds[i]," ", Cos2Indvs[i][:])
        fichierResultat.write(nomInds[i]+" "+ str(Cos2Indvs[i][:]) + "\n")
    
    
    print("----------------------------------------------------------------------")
    fichierResultat.write("----------------------------------------------------------------------\n")
    
    
    # Pour les Vars
        # coord <==> matricedesaxes
    
    MatriceDesAxesVars =[] 
    
    for i in range(nbValsPrps):
        MatriceDesAxesVars.append(list())
        for j in range(nbValsPrps):
            tmp=0
            for k in range(nbValsPrps):
                #if( i == 0 and j==0 ):
                    #print("matriceCentreReduite[i+1][k+1] ",matriceCentreReduite[i+1][k+1]," * ",VectsETValsPropres[j][1][k])
                tmp += matriceDeCorrelation[i][k] * VectsETValsPropres[j][1][k]
            MatriceDesAxesVars[i].append(tmp)
    
    print("----------------------------------------------------------------------")
    fichierResultat.write("----------------------------------------------------------------------\n")
    print("Coordonnées des Variables")
    fichierResultat.write("Coordonnées des Variables\n")
    print(tmpNomsAxes)
    fichierResultat.write(str(tmpNomsAxes)+"\n")
    
    for i in range(nbVars):
        print(nomVars[i+1]," ", MatriceDesAxesVars[i][:])
        fichierResultat.write(nomVars[i+1]+" "+ str(MatriceDesAxesVars[i][:])+"\n")
            
    print("----------------------------------------------------------------------")
    fichierResultat.write("----------------------------------------------------------------------\n")
    
        # Contrib vars
    SommeCarreCoordVars = []
    
    for j in range(nbVars):
        
        tmp = 0
        for i in range(nbVars):
            tmp+= MatriceDesAxesVars[i][j]**2
        SommeCarreCoordVars.append(tmp)
        
    
    
    ContribVars = []
    
    for i in range(nbVars):
        ContribVars.append(list())
        for j in range(nbVars):
            ContribVars[i].append(((MatriceDesAxesVars[i][j]**2)/SommeCarreCoordVars[j])*100)
    
    print("----------------------------------------------------------------------")
    fichierResultat.write("----------------------------------------------------------------------\n")
    print("Contribution des Variables")
    fichierResultat.write("Contribution des Variables\n")
    print(tmpNomsAxes)
    
    for i in range(nbVars):
        print(nomVars[i+1]," ", ContribVars[i][:])
        fichierResultat.write(nomVars[i+1]+" "+ str(ContribVars[i][:])+"\n")
            
    print("----------------------------------------------------------------------")
    fichierResultat.write("----------------------------------------------------------------------\n")
        # cos2 Vars
    
    SommeCarreCos2Vars = []
    
    for i in range(nbVars):
        
        tmp = 0
        for j in range(nbVars):
            tmp+= MatriceDesAxesVars[i][j]**2
        SommeCarreCos2Vars.append(tmp)
        
    
    
    Cos2Vars = []
    
    for i in range(nbVars):
        Cos2Vars.append(list())
        for j in range(nbVars):
            Cos2Vars[i].append(((MatriceDesAxesVars[i][j]**2)/SommeCarreCos2Vars[i])*100)
    
    print("----------------------------------------------------------------------")
    fichierResultat.write("----------------------------------------------------------------------\n")
    print("Qualité de représentation des Variables")
    fichierResultat.write("Qualité de représentation des Variables\n")
    print(tmpNomsAxes)
    fichierResultat.write(str(tmpNomsAxes)+"\n")
    
    for i in range(nbVars):
        print(nomVars[i+1]," ", Cos2Vars[i][:])
        fichierResultat.write(nomVars[i+1]+" "+ str(Cos2Vars[i][:])+"\n")
            
    print("----------------------------------------------------------------------")
    fichierResultat.write("----------------------------------------------------------------------\n")
        
    

matriceTest = [
["None","PA15+","PA15+Urb","PA15+Rur","PAO","PAOUrb","PAORur","AgrForPêche","Indus","BatTrav","comm","TraEntCommu","Serv","MalDes","PAC","TFPAC","ChSansDipUrb","ChNivMoyUrb","ChNivSupUrb","ChDipUrb","ChSansDipRur","ChDipRur"],
["2002-08",10849,5532,5318,9710,4586,5124,5.3,21.9,10.0,20.9,6.2,35.6,0.1,1139,27.7,9.8,23.6,23.5,23.7,2.4,11.6],
["2008",11267,5874,5393,10189,5013,5176,5.5,20.9,11.2,19.9,6.8,35.6,0.2,1078,27.5,8.3,20.0,20.6,19.0,2.6,12.2],
["2009",11314,5916,5398,10284,5101,5184,5.0,20.2,11.8,19.9,6.6,36.0,0.2,1029,27.6,7.7,18.6,19.1,17.8,2.5,11.8],
["2010",11415,5966,5449,10405,5169,5235,4.8,20.2,12.4,20.2,6.9,35.3,0.2,1037,28.3,8.1,18.1,18.5,17.5,2.4,11.4],
["2011",11538,5553,5237,10509,5272,5237,4.9,20.2,12.4,20.6,7.1,35.1,0.2,1028,30.6,7.0,18.3,17.8,19.0,2.3,11.1],
["2012",11549,6655,5404,10511,5320,5190,5.1,18.7,11.9,20.5,6.8,36.8,0.2,1038,29.1,6.9,18.2,17.9,18.6,2.4,10.6],
["2013",11705,6217,5488,10625,5346,5278,4.9,18.4,10.9,21.0,6.6,38.0,0.2,1081,27.8,8.1,18.2,18.1,18.4,2.4,9.8]
]

started = time.time()

ACP(matriceTest)

print("temps d'execution : ",time.time()-started)
