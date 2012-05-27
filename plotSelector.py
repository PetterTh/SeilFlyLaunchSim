#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      skatun
#
# Created:     25/05/2012
# Copyright:   (c) skatun 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from selFunc import *
import SeilSim as ss
import numpy as np
from pylab import *

figNum = 1
run = 2
def plotSelector():
    global figNum,logg,scalePlot,fileName,x,y,yLabel,titleLabel,paraMeterArray
    global saveOn,showOn,run
    paraMeterArray = ss.init()
    print ss.simulate(paraMeterArray)
    logg,setArr =ss.saveLogg()
    inp = setArr['plotVer']
    plotKey = setArr['plotKey']
    plotKeyArray = setArr['plotKeyArray']
    loopArray = logg
    colors = ["r","b","m","y","black","o"]
    showOn = setArr['showOn']
    saveOn = setArr['saveOn']
    scalePlot = setArr['scalePlot']
    numbers = setArr['numbersOfVariationStepForPlots']
    counterPhase = setArr['counterPhase']
    xLabel = "X-Position [m]"
    yLabel = "Y-Position [m]"
    titleLabel = 'Height plot'
    x = logg['x']
    y = logg['y']
    plotPhases = 1
    run = run-1

    if inp == 0:
        showOn = 0
        saveOn = 1
        inp = 6
        if run:
            ss.sensitivity(paraMeterArray)
            loopArray = paraMeterArray
    elif inp==1: # Single normal x,y plot
        fileName = 'Figures/Single/XY-Plot-'  + '.png'
    elif inp==2:
        generateYArray(plotKey)
    elif inp==5 or inp==6:
        loopArray = paraMeterArray
    elif inp==7:
        ss.sensitivity(plotKeyArray)
    else:
        fileName = 'Figures/Single/XY-Plot-'  + '.png'

    for myKey, value in loopArray.items():
        figure(figNum)
        if inp == 3 or inp==4:
            scalePlot = 1
            generateYArray(myKey)


        if (inp==4 and myKey in plotKeyArray) or inp <= 3 or (inp==5 and myKey==plotKey) or inp ==6 or inp==0:

            if inp == 5 or inp == 6:
                x,y = varyVariabel(myKey,paraMeterArray[myKey]*0.7,paraMeterArray[myKey]*1.3,numbers)
                plotPhases = 0
                fileName = 'Figures/Vary/'+ myKey +'Y-Plot-'  + '.png'
                yLabel = "Height [m]"
                legend(['Phase 2','Phase 3','Phase 4'])
            if plotPhases:
                for index in range(0,len(counterPhase)-1):
                    plot(x[counterPhase[index]:counterPhase[index+1]],y[counterPhase[index]:counterPhase[index+1]],colors[index])
            else:
                plot(x,y)

            xlabel(xLabel)
            ylabel(yLabel)
            title(titleLabel)
            grid()

            if saveOn:
                savefig(fileName)

            figNum = figNum +1

            if inp<= 2:
                break

    if showOn:
        show()

    if run and inp == 0:
        plotSelector()

    print run

def generateYArray(key):
    global y,fileName,yLabel,titleLabel,scalePlot
    y = array(logg[key]) * scalePlot
    fileName = 'Figures/Single/X' + key + '-Plot-'  + '.png'
    yLabel = key
    titleLabel = key + ' plot'
    # getDateString()

def varyVariabel(key,startValue,stopValue,numbers):
    global paraMeterArray
    ss.init()
    stepSize = ((float(stopValue-startValue)/numbers))
    x = []
    y = []
    i = startValue
    while i <stopValue:
        paraMeterArray[key] = i
        ss.setParametersArray(paraMeterArray)
        x.append(i)
        y.append(ss.simulate(paraMeterArray))
        i = i + stepSize
    return x,y

def sensitivity(keyArray):
    #global paraMeterArray
    #planeParameters0,flightParameters0, winchParameters0,lineParameters0,flighConditionsParameters0
    #paraMeterArray = init()
    #testArray = paraMeterArray.copy()
    tempArray = {}
    for i in keyArray:
        tempArray[i]= paraMeterArray[i]

    keys,minVal,maxVal,minRes,maxRes = sensitivityCheck(tempArray)
    sensitivityMy,state,keys2 = sensitivityDelta(keys,minVal,maxVal,minRes,maxRes)
    plotSensitivity(sensitivityMy,keys2,state,1)
    #printSensitivity(keys,minVal,maxVal,minRes,maxRes)

def sensitivityDelta(keys,minVal,maxVal,minRes,maxRes):

    maxValues = max(max(minRes),max(maxRes))

    tempArray = abs(subtract(minRes,maxRes))
    tempArray2 = []
    keys2 = []
    s = ''
    s2 = ''
    test =len(tempArray)
    for i in range(0,len(tempArray)):
        if tempArray[i,0] == 0:
            s = s + keys[i] + ', '
        elif tempArray[i,0] > 0 and tempArray[i,0] < .05*max(tempArray[i,:]):
            s2 = s2 + keys[i] + ', '
        else:
            tempArray2.append(tempArray[i])
            keys2.append(keys[i])

    tempArray2 = array(tempArray2)
    if len(s)>0:
        s = s + 'not implemented'
    if len(s2)>0 and len(s)>0:
        s = s + '\n'
    if len(s2)>0:
        s= s + s2 +'minor effect on result'

    #senstivityArrayTemp= (abs(subtract(minRes,maxRes))/maxValues)
    senstivityArrayTemp= tempArray2/maxValues
    summOfArray = senstivityArrayTemp.sum(axis=0)

    return senstivityArrayTemp/summOfArray,s,keys2

def plotSensitivity(sensitivity,keys,state,exploded):
    global figNum,fileName

    for i in range(0,3,2):

        sensitivityCol = sensitivity[:,i]

        if exploded:
            explodedMy = np.zeros(len(sensitivityCol))
            minindex = max((v,i) for i,v in enumerate(sensitivityCol))[1]
            explodedMy[minindex] = 0.05
        else:
            explodedMy = None

        figure(figNum,figsize=(8,8))
        title('Height in phase ' + str(i+2))
        ax = axes([0.1, 0.1, 0.8, 0.8])
        pie(sensitivityCol,explode=explodedMy, labels=keys)
        if len(state)>0:
            text(-1.2, -1.2, state, bbox=dict(facecolor='red', alpha=0.1))

        fileName = 'Figures/Sensitivity/Height' + str(i+2) + '-Plot-'  + '.png'
        figNum = figNum + 1

        if saveOn:
            savefig(fileName)
        if showOn:
            show()

def sensitivityCheck(changeArray0):

    paraMeterArray = ss.init()

    keys   = []
    minVal = []
    maxVal = []
    minRes = []
    maxRes = []

    for key, value in changeArray0.items():
        paraMeterArray[key] = value*.9
        minVal.append(paraMeterArray[key])
        minRes.append(ss.simulate(paraMeterArray))

        paraMeterArray[key] = value*1.1
        maxVal.append(paraMeterArray[key])
        maxRes.append(ss.simulate(paraMeterArray))
        keys.append(key)

    return keys,minVal,maxVal,minRes,maxRes


def printSensitivity(keys,minVal,maxVal,minRes,maxRes):

    sentivity,s = sensitivityDelta(keys,minVal,maxVal,minRes,maxRes)*100
    print "Parameter    Height2 Height3 Height5 Min Val Max Val "

    for i in range(0,len(sentivity)):
        print keys[i],sentivity[i],"%",minVal[i],maxVal[i]


def saveFigMy(fileName):
    _dictonary = getParametersArray()
    myString = getDateString()
    #savefig('Figures/' + myString + '.png')
    writeDict(_dictonary,fileName)