#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      skatun
#
# Created:     20/05/2012
# Copyright:   (c) skatun 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from selFunc import *
from SeilSim import *
from plotSensitivity import *


def main():
    planeParameters,flightParameters,winchParameters,lineParameters,flighConditionsParameters = init()
    reset(planeParameters,flightParameters,winchParameters,lineParameters,flighConditionsParameters)
    loggingReset()
    print simulate([3])

def testSensitivity():
    planeParameters,flightParameters,winchParameters,lineParameters,flighConditionsParameters = init()
    reset(planeParameters,flightParameters,winchParameters,lineParameters,flighConditionsParameters)
    keys,minVal,maxVal,minRes,maxRes = sensitivityCheck(planeParameters)
    sensitivity = sensitivityDelta(keys,minVal,maxVal,minRes,maxRes)


    plotSensitivity(sensitivity,keys,1,1,1)
    printSensitivity(keys,minVal,maxVal,minRes,maxRes)

def test():
    y = 150
    p0 = p0ISA()
    T0 = T0ISA()
    print density(p0,T0,y)
    print densityWithHumidity(0,p0,T0,y)


def test3():
    print getDateString()

if __name__ == '__main__':
    #testSensitivity()
    test3()
