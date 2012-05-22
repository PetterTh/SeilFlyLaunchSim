#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      skatun
#
# Created:     21/05/2012
# Copyright:   (c) skatun 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import numpy as np
from pylab import *
from SeilSim import *
from selFunc import *



def plotSensitivity(sensitivity,keys,exploded,showMy,save):
##    pie(sensitivity, explode=None, labels=keys)


    for i in range(0,2):
        figureNumber = getFigureNumber()
        sensitivityCol = sensitivity[:,i]

        if exploded:
            explodedMy = np.zeros(len(sensitivityCol))
            minindex = max((v,i) for i,v in enumerate(sensitivityCol))[1]
            explodedMy[minindex] = 0.05
        else:
            explodedMy = None

        figure(figureNumber,figsize=(8,8))
        ax = axes([0.1, 0.1, 0.8, 0.8])
        pie(sensitivityCol,explode=explodedMy, labels=keys)
        setFigureNumber()
    show()

