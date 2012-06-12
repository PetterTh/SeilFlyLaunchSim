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
legendLabel = []
numPlot = 0
newFig = 1
exclusiveSen = 0
senstivityName = ''
def plotSelector():
    """
    0: All available plot methods
    1: x,y plot

    2: x,plotkey
    3: x, all pararmeters in logg array
    4: x, plotKeyArray

    5: plotKey, height
    6: all pararmeters in paraMeterArray, height
    7: plotKeyArray,height

    8: all sensitivity plot
    9: planeParameters sensitivity plot
    10: flightParameters sensitivity plot
    11: winchParameters sensitivity plot
    12: lineParameters sensitivity plot
    13: flighConditionsParameters sensitivity plot
    14: plotKeyArray sensitivity plot

    """
    global figNum,logg,scalePlot,fileName,x,y,yLabel,titleLabel,paraMeterArray
    global saveOn,showOn,run,legendLabel,numPlot,newFig,exclusive,exclusiveSen
    global plotKeyArray,inp,senstivityName

    paraMeterArray = ss.init()
    setArr = ss.getPlotSettings()
    inp = setArr['plotVer']
    exclusive = setArr['exclusive']
    plotKey = setArr['plotKey']
    plotKeyArray = setArr['plotKeyArray']
    colors = ["r","b","m","y","black","o"]
    showOn = setArr['showOn']
    saveOn = setArr['saveOn']
    scalePlot = setArr['scalePlot']
    numbers = setArr['numbersOfVariationStepForPlots']
    counterPhase = setArr['counterPhase']
    xLabel = "X-Position [m]"
    yLabel = "Y-Position [m]"
    titleLabel = 'Height plot'
    colorIndex = 0

    plotPhases = 1
    run = run-1
    legendOn = 0
    loopInp = 1

    if inp == 0:
        showOn = 0
        saveOn = 1
        inp = 6
        if run:
            ss.sensitivity(paraMeterArray)
            loopInp = 2

    elif inp==1: # Single normal x,y plot
        fileName = 'Figures/Single/XY-Plot-'  + '.png'
    elif inp==2:
        loopInp = 3
    elif inp==5 or inp==6 or inp == 7:
        loopInp = 2
    elif inp==8:
        sensitivity(paraMeterArray)
    elif inp==9:
        sensitivity(ss.getPlaneParameters())
    elif inp==10:
        sensitivity(ss.getFlightParameters())
    elif inp==11:
        sensitivity(ss.getWinchParameters())
    elif inp==12:
        sensitivity(ss.getLineParameters())
    elif inp==13:
        sensitivity(ss.getFlighConditionsParameters())
    elif inp==14:
        exclusiveSen = exclusive
        sensitivity(plotKeyArray)
    elif inp==17:
        senstivityName = 'planePara'
        sensitivity(ss.getPlaneParameters())
        senstivityName = 'flightPara'
        sensitivity(ss.getFlightParameters())
        senstivityName = 'winchPara'
        sensitivity(ss.getWinchParameters())
        senstivityName = 'linePara'
        sensitivity(ss.getLineParameters())
        senstivityName = 'flightCondPara'
        sensitivity(ss.getFlighConditionsParameters())

    elif inp==15:
        colorIndex = numPlot
        value = paraMeterArray[plotKey]
        stepSize = ((float(value*1.5-value*.5)/numbers))
        if run and (value*1.5 - (value*0.5+(numPlot+1)*stepSize))>0.00001 :
            paraMeterArray[plotKey] = value*0.5+numPlot*stepSize
            run = run+1
            saveOn = 0
            showOn = 0
            if newFig:
                newFig=0
        else:
            paraMeterArray[plotKey] = value*0.7+numPlot*stepSize
            legendOn = 1
            numPlot= -1

        plotPhases = 0
        legendLabel.append(plotKey + ':' + str(round(paraMeterArray[plotKey]*scalePlot,2)))
        numPlot = numPlot + 1
        fileName = 'Figures/Alt/' + plotKey + '-Plot-'  + '.png'

    elif inp==16:
        plotPhases = 0
        colorIndex = numPlot

        if run:
            saveOn = 0
            showOn = 0
            if newFig:
                newFig=0
                paraMeterArray['preTensionOfLine'] = 0
                run = run +1
        else:
            legendOn = 1

        ss.setAlt(numPlot+1)
        if numPlot==0:
            legendLabel.append('Rolling on ground')
        if numPlot==1:
            legendLabel.append('Thrown in the air')
        if numPlot==2:
            legendLabel.append('No zooming')

        numPlot = numPlot + 1

        fileName = 'Figures/Alt/Alt-Plot-'  + '.png'

    else:
        fileName = 'Figures/Single/XY-Plot-'  + '.png'


    ss.simulate(paraMeterArray)
    logg =ss.saveLogg()
    loopArray = logg
    counterPhase = ss.getCounterPhases()

    x = logg['x']
    y = logg['y']

    if loopInp ==1:
        loopArray = logg
    elif loopInp == 2:
        loopArray = paraMeterArray
    elif loopInp == 3:
        generateYArray(plotKey)

    for myKey, value in loopArray.items():
        plotThisKey = 0
        if inp == 3 or inp==4:
            scalePlot = 1
            generateYArray(myKey)
        if myKey in plotKeyArray:
            plotThisKey = 1
        if exclusive:
            plotThisKey = abs(plotThisKey-1)

        if inp <= 3 or (inp==4 and plotThisKey ) or (inp==5 and myKey in plotKey) or inp==6 or (inp==7 and  myKey in plotKeyArray) or inp >=15 :
            if (run and inp>=15 and newFig) or inp <=14:
                figure(figNum)
                hold(True)  # hold is on

            if inp == 5 or inp == 6 or inp == 7:
                x,y = varyVariabel(myKey,paraMeterArray[myKey]*0.7,paraMeterArray[myKey]*1.3,numbers)
                if inp==5:
                    x=array(x)*scalePlot
                plotPhases = 0
                fileName = 'Figures/Vary/'+ myKey +'Y-Plot-'  + '.png'
                titleLabel = 'Height plot for differnet values of ' + myKey
                xLabel = myKey
                yLabel = "Height [m]"
                legendLabel = ['Phase 2','Phase 3','Phase 4','Phase 5']
                legendOn = 1
                colorIndex = [1,2,3,4]
            if plotPhases:
                for index in range(0,len(counterPhase)-1):
                    plot(x[counterPhase[index]:counterPhase[index+1]],y[counterPhase[index]:counterPhase[index+1]],colors[index])
            elif inp == 5 or inp == 6 or inp == 7 or inp==15:
                plot(x,y)
            else:
                plot(x,y,colors[colorIndex])

            xlabel(xLabel)
            ylabel(yLabel)
            title(titleLabel)
            grid()

            if legendOn:
                legend(legendLabel,loc='upper left')

            if saveOn:
                savefig(fileName)
            if newFig:
                figNum = figNum +1

            if inp<= 2 or inp>=15:
                break

    if showOn:
        show()

    if run and (inp == 0 or inp >= 16 or numPlot):
        plotSelector()



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
    tempArray = {}
    if exclusive and inp!=14:
        for key in keyArray:
            plotThisKey = 0
            if key in plotKeyArray:
                plotThisKey = 1
            if exclusive:
                plotThisKey = abs(plotThisKey-1)
            if plotThisKey:
                tempArray[key]= paraMeterArray[key]
    else:
        for i in keyArray:
            tempArray[i]= paraMeterArray[i]


    keys,minVal,maxVal,minRes,maxRes = sensitivityCheck(tempArray)
    sensitivityMy,state,keys2 = sensitivityDelta(keys,minVal,maxVal,minRes,maxRes)
    plotSensitivity(sensitivityMy,keys2,state,1)
    #printSensitivity(keys,minVal,maxVal,minRes,maxRes)

def sensitivityDelta(keys,minVal,maxVal,minRes,maxRes):
    temp = zeros([2,4], float)
    temp[0,:] = (minRes.max(axis=0))
    temp[1,:] = (maxRes.max(axis=0))
    maxValues = temp.max(axis=0)

    tempArray = abs(subtract(minRes,maxRes))
    tempArray2 = []
    keys2 = []
    s = ''
    s2 = ''
    test =len(tempArray)
    # Here we loop through the difference array to check if we do have any differences..
    # if there is no differences the parameter is most likely not be implemented
    # if there is minor differences its taged minor and removed from array to make pie plot nicer
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

        fileName = 'Figures/Sensitivity/Height' + str(i+2) + '-Plot-' + senstivityName + '.png'
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


    return keys,minVal,maxVal,array(minRes),array(maxRes)



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

def writeLatexTable(_dict):
    dictlist = []
    for key, value in _dict.iteritems():
        temp = [key,value]
        dictlist.append(temp)

    f = open('test.txt','w')

    print f,"""
    \documentclass{article}

    \usepackage{siunitx}

    \begin{document}
        \begin{table}
            \begin{center}
            \begin{tabular}{SSSS}
    """

    print f, " \\\\\n".join([" & ".join(map(str,line)) for line in dictlist])

    print """
            \end{tabular}
      \end{center}
    \end{table}
    \end{document}
    """

    f.close()

