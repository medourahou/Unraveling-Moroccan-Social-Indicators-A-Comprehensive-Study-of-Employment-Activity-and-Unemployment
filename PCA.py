import math
import numpy as np
from numpy.linalg import eig
import time

def PCA(matrix):
    resultFile = open("./PCA_result.txt", "w")

    # Calculate means and standard deviations for each column (variable)
    numVars = len(matrix[0]) - 1  # The first column contains the names of individuals
    numInds = len(matrix) - 1  # The first row contains the names of variables

    varNames = matrix[0][:]
    
    indNames = []
    for i in range(1, numInds + 1):
        indNames.append(matrix[i][0])

    means = []
    std_deviations = []

    for j in range(1, numVars + 1):
        means.append(0)
        for i in range(1, numInds + 1):
            means[j - 1] += matrix[i][j]
        means[j - 1] /= numInds

    for j in range(1, numVars + 1):
        std_deviations.append(0)
        for i in range(1, numInds + 1):
            std_deviations[j - 1] += (matrix[i][j] - means[j - 1]) ** 2
        std_deviations[j - 1] /= numInds
        std_deviations[j - 1] = math.sqrt(std_deviations[j - 1])

    # Correlation matrix

    correlationMatrix = []

    for j in range(0, numVars):
        correlationMatrix.append(list())
        for k in range(0, numVars):
            correlationMatrix[j].append(0)
            for i in range(0, numInds):
                correlationMatrix[j][k] += (matrix[i + 1][j + 1] - means[j]) * (matrix[i + 1][k + 1] - means[k])

            correlationMatrix[j][k] /= numInds
            correlationMatrix[j][k] /= std_deviations[j] * std_deviations[k]

            if k == j:
                correlationMatrix[j][k] = 1.0

    # Display the correlation matrix

    print("----------------------------------------------------------------------")
    resultFile.write("----------------------------------------------------------------------\n")
    print("Correlation matrix of variables")
    resultFile.write("Correlation matrix of variables\n")

    tmpVarNames = ""
    for i in range(numVars + 1):
        if i == 0:
            tmpVarNames += " -------- "
        else:
            tmpVarNames += " " + varNames[i] + " -------- "

    print(tmpVarNames)
    resultFile.write(tmpVarNames + "\n")

    for i in range(numVars):
        print(varNames[i + 1], " ", correlationMatrix[i][:])
        tmpp = varNames[i + 1] + " " + str(correlationMatrix[i][:]) + "\n"
        resultFile.write(tmpp)

    print("----------------------------------------------------------------------")
    resultFile.write("----------------------------------------------------------------------\n")

    # Center and reduce the data

    standardizedMatrix = matrix[:]

    totalInertia = numVars

    for j in range(1, numVars + 1):
        for i in range(1, numInds + 1):
            tmp = (standardizedMatrix[i][j] - means[j - 1]) / std_deviations[j - 1]
            standardizedMatrix[i][j] = tmp

    correlationMatrixNumpy = [None] * numVars
    for i in range(0, numVars):
        correlationMatrixNumpy[i] = np.array(correlationMatrix[i][:])

    correlationMatrixNumpy = np.array(correlationMatrixNumpy)

    valsEigen, eigenvectors = eig(correlationMatrixNumpy)
    """Note: For the same eigenvalue, several eigenvectors may be associated, 
    so the axes may differ as they depend on the chosen eigenvector"""

    # Convert to Python lists

    valsEigen = list(valsEigen)

    valsEigen = [float(v) for v in valsEigen]

    eigenvectors = [list(vect) for vect in eigenvectors]

    for vect in eigenvectors:
        vect = [float(v) for v in vect]

    VectsAndValsEigen = list()

    numValsEigen = len(valsEigen)

    for i in range(numValsEigen):
        VectsAndValsEigen.append((float(valsEigen[i]), [float(v) for v in eigenvectors[i]]))

    # Sort VectsAndValsEigen in descending order based on eigenvalues

    VectsAndValsEigen = [elm for elm in sorted(VectsAndValsEigen, key=lambda elm: elm[0], reverse=True)]

    print("----------------------------------------------------------------------")
    resultFile.write("----------------------------------------------------------------------\n")
    print("Eigenvalues and eigenvectors")
    resultFile.write("Eigenvalues and eigenvectors\n")

    for i in range(numValsEigen):
        print("Eigenvalue", i, " : ", VectsAndValsEigen[i][0])
        tmpp = "Eigenvalue" + str(i) + " : " + str(VectsAndValsEigen[i][0]) + "\n"
        resultFile.write(tmpp)
        print("Eigenvector", i, " : ", VectsAndValsEigen[i][1])
        tmpp = "Eigenvector" + str(i) + " : " + str(VectsAndValsEigen[i][1]) + "\n"
        resultFile.write(tmpp)

    print("----------------------------------------------------------------------")
    resultFile.write("----------------------------------------------------------------------\n")

    # Calculate the coordinates of individuals in the new space (axes)

    AxesMatrix = []

    for i in range(numInds):
        AxesMatrix.append(list())
        for j in range(numValsEigen):
            tmp = 0
            for k in range(numValsEigen):
                tmp += standardizedMatrix[i + 1][k + 1] * VectsAndValsEigen[j][1][k]
            AxesMatrix[i].append(tmp)

    # Display the coordinates of individuals

    axisNames = [("Axis " + str(i + 1) + " ") for i in range(numValsEigen)]
    print("----------------------------------------------------------------------")
    resultFile.write("----------------------------------------------------------------------\n")
    print("Coordinates of individuals")
    resultFile.write("Coordinates of individuals\n")

    tmpAxisNames = "--------"
    for i in range(numVars):
        tmpAxisNames += " " + axisNames[i] + " -------- "

    print(tmpAxisNames)
    resultFile.write(str(tmpAxisNames) + "\n")

    for i in range(numInds):
        print(indNames[i], " ", AxesMatrix[i][:])
        tmpp = indNames[i] + " " + str(AxesMatrix[i][:]) + "\n"
        resultFile.write(tmpp)

    print("----------------------------------------------------------------------")
    resultFile.write("----------------------------------------------------------------------\n")

    # Calculate the percentage of inertia for eigenvalues

    inertiaPercentages = [None] * 3

    inertiaPercentages[0] = [val[0] for val in VectsAndValsEigen]

    inertiaPercentages[1] = [val[0] / totalInertia for val in VectsAndValsEigen]

    inertiaPercentages[2] = []

    for i in range(0, numValsEigen):
        if i == 0:
            inertiaPercentages[2].append(inertiaPercentages[1][0])
        else:
            inertiaPercentages[2].append(inertiaPercentages[2][i - 1] + inertiaPercentages[1][i])

    # Display cumulative inertia percentages

    print("----------------------------------------------------------------------")
    resultFile.write("----------------------------------------------------------------------\n")
    print("Percentage of inertia for eigenvalues and cumulative inertia percentage")
    resultFile.write("Percentage of inertia for eigenvalues and cumulative inertia percentage\n")

    print("Eigenvalues: ", inertiaPercentages[0][:])
    resultFile.write("Eigenvalues: " + str(inertiaPercentages[0][:]) + "\n")

    print("% inertia: ", inertiaPercentages[1][:])
    resultFile.write("% inertia: " + str(inertiaPercentages[1][:]) + "\n")

    print("% cumulative inertia: ", inertiaPercentages[2][:])
    resultFile.write("% cumulative inertia: " + str(inertiaPercentages[2][:]) + "\n")

    print("----------------------------------------------------------------------")
    resultFile.write("----------------------------------------------------------------------\n")

    # Individual contributions

    sumSquareCoords = []

    for j in range(numVars):

        tmp = 0
        for i in range(numInds):
            tmp += AxesMatrix[i][j] ** 2
        sumSquareCoords.append(tmp)

    contributionsInd = []

    for i in range(numInds):
        contributionsInd.append(list())
        for j in range(numVars):
            contributionsInd[i].append(((AxesMatrix[i][j] ** 2) / sumSquareCoords[j]) * 100)

    # Display contributions of individuals

    print("----------------------------------------------------------------------")
    resultFile.write("----------------------------------------------------------------------\n")
    print("Contributions of individuals")
    resultFile.write("Contributions of individuals\n")
    print(tmpAxisNames)
    resultFile.write(str(tmpAxisNames) + "\n")

    for i in range(numInds):
        print(indNames[i], " ", contributionsInd[i][:])
        resultFile.write(indNames[i] + " " + str(contributionsInd[i][:]) + "\n")

    print("----------------------------------------------------------------------")
    resultFile.write("----------------------------------------------------------------------\n")

    # Cos2 of individuals

    sumSquareCos2 = []

    for i in range(numInds):

        tmp = 0
        for j in range(numVars):
            tmp += AxesMatrix[i][j] ** 2
        sumSquareCos2.append(tmp)

    cos2Ind = []

    for i in range(numInds):
        cos2Ind.append(list())
        for j in range(numVars):
            cos2Ind[i].append(((AxesMatrix[i][j] ** 2) / sumSquareCos2[i]) * 100)

    # Display qualities of representation

    print("----------------------------------------------------------------------")
    resultFile.write("----------------------------------------------------------------------\n")
    print("Qualities of representation of individuals")
    resultFile.write("Qualities of representation of individuals\n")
    print(tmpAxisNames)
    for i in range(numInds):
        print(indNames[i], " ", cos2Ind[i][:])
        resultFile.write(indNames[i] + " " + str(cos2Ind[i][:]) + "\n")

    print("----------------------------------------------------------------------")
    resultFile.write("----------------------------------------------------------------------\n")

    # For Variables
    # Coords <==> AxesMatrixVars

    AxesMatrixVars = []

    for i in range(numValsEigen):
        AxesMatrixVars.append(list())
        for j in range(numValsEigen):
            tmp = 0
            for k in range(numValsEigen):
                tmp += correlationMatrix[i][k] * VectsAndValsEigen[j][1][k]
            AxesMatrixVars[i].append(tmp)

    print("----------------------------------------------------------------------")
    resultFile.write("----------------------------------------------------------------------\n")
    print("Coordinates of Variables")
    resultFile.write("Coordinates of Variables\n")
    print(tmpAxisNames)
    resultFile.write(str(tmpAxisNames) + "\n")

    for i in range(numVars):
        print(varNames[i + 1], " ", AxesMatrixVars[i][:])
        resultFile.write(varNames[i + 1] + " " + str(AxesMatrixVars[i][:]) + "\n")

    print("----------------------------------------------------------------------")
    resultFile.write("----------------------------------------------------------------------\n")

    # Contributions of Variables

    sumSquareCoordsVars = []

    for j in range(numVars):

        tmp = 0
        for i in range(numVars):
            tmp += AxesMatrixVars[i][j] ** 2
        sumSquareCoordsVars.append(tmp)

    contributionsVars = []

    for i in range(numVars):
        contributionsVars.append(list())
        for j in range(numVars):
            contributionsVars[i].append(((AxesMatrixVars[i][j] ** 2) / sumSquareCoordsVars[j]) * 100)

    print("----------------------------------------------------------------------")
    resultFile.write("----------------------------------------------------------------------\n")
    print("Contributions of Variables")
    resultFile.write("Contributions of Variables\n")
    print(tmpAxisNames)

    for i in range(numVars):
        print(varNames[i + 1], " ", contributionsVars[i][:])
        resultFile.write(varNames[i + 1] + " " + str(contributionsVars[i][:]) + "\n")

    print("----------------------------------------------------------------------")
    resultFile.write("----------------------------------------------------------------------\n")

    # Cos2 of Variables

    sumSquareCos2Vars = []

    for i in range(numVars):

        tmp = 0
        for j in range(numVars):
            tmp += AxesMatrixVars[i][j] ** 2
        sumSquareCos2Vars.append(tmp)

    cos2Vars = []

    for i in range(numVars):
        cos2Vars.append(list())
        for j in range(numVars):
            cos2Vars[i].append(((AxesMatrixVars[i][j] ** 2) / sumSquareCos2Vars[i]) * 100)

    print("----------------------------------------------------------------------")
    resultFile.write("----------------------------------------------------------------------\n")
    print("Qualities of representation of Variables")
    resultFile.write("Qualities of representation of Variables\n")
    print(tmpAxisNames)
    resultFile.write(str(tmpAxisNames) + "\n")

    for i in range(numVars):
        print(varNames[i + 1], " ", cos2Vars[i][:])
        resultFile.write(varNames[i + 1] + " " + str(cos2Vars[i][:]) + "\n")

    print("----------------------------------------------------------------------")
    resultFile.write("----------------------------------------------------------------------\n")


matrixTest = [
    ["None", "PA15+", "PA15+Urb", "PA15+Rur", "PAO", "PAOUrb", "PAORur", "AgrForPêche", "Indus", "BatTrav", "comm",
     "TraEntCommu", "Serv", "MalDes", "PAC", "TFPAC", "ChSansDipUrb", "ChNivMoyUrb", "ChNivSupUrb", "ChDipUrb",
     "ChSansDipRur", "ChDipRur"],
    ["2002-08", 10849, 5532, 5318, 9710, 4586, 5124, 5.3, 21.9, 10.0, 20.9, 6.2, 35.6, 0.1, 1139, 27.7, 9.8, 23.6, 23.5,
     23.7, 2.4, 11.6],
    ["2008", 11267, 5874, 5393, 10189, 5013, 5176, 5.5, 20.9, 11.2, 19.9, 6.8, 35.6, 0.2, 1078, 27.5, 8.3, 20.0, 20.6, 19.0,
     2.6, 12.2],
    ["2009", 11314, 5916, 5398, 10284, 5101, 5184, 5.0, 20.2, 11.8, 19.9, 6.6, 36.0, 0.2, 1029, 27.6, 7.7, 18.6, 19.1, 17.8,
     2.5, 11.8],
    ["2010", 11415, 5966, 5449, 10405, 5169, 5235, 4.8, 20.2, 12.4, 20.2, 6.9, 35.3, 0.2, 1037, 28.3, 8.1, 18.1, 18.5, 17.5,
     2.4, 11.4],
    ["2011", 11538, 5553, 5237, 10509, 5272, 5237, 4.9, 20.2, 12.4, 20.6, 7.1, 35.1, 0.2, 1028, 30.6, 7.0, 18.3, 17.8, 19.0,
     2.3, 11.1],
    ["2012", 11549, 6655, 5404, 10511, 5320, 5190, 5.1, 18.7, 11.9, 20.5, 6.8, 36.8, 0.2, 1038, 29.1, 6.9, 18.2, 17.9, 18.6,
     2.4, 10.6],
    ["2013", 11705, 6217, 5488, 10625, 5346, 5278, 4.9, 18.4, 10.9, 21.0, 6.6, 38.0, 0.2, 1081, 27.8, 8.1, 18.2, 18.1, 18.4,
     2.4, 9.8]
]

started = time.time()

ACP(matrixTest)

print("Execution time: ", time.time() - started)
