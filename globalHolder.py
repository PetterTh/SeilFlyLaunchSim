#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      skatun
#
# Created:     29/06/2012
# Copyright:   (c) skatun 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
plotVer = 1
plotKey = 'energy'
counterPhase = []


def getPlotVer():
    return plotVer

def setPlotVer(_plotVer):
    global plotVer
    plotVer = _plotVer

def getPlotKey():
    return plotKey

def setPlotKey(_plotKey):
    global plotKey
    plotKey = _plotKey

def getCounterPhases():
    return counterPhase

def setCounterPhase(_counterPhase):
    global counterPhase
    counterPhase = _counterPhase