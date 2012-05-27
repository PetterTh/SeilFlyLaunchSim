#-------------------------------------------------------------------------------
# Name: selFunc
# Purpose: Contains a selection of simple functions
#
# Author: Petter
#
# Created: 14.05.2012
# Copyright: (c) Petter 2012
# Licence: <your licence>
#-------------------------------------------------------------------------------

import numpy as np
import re
from datetime import datetime
#from pylab import *


figureNumber = 1
paraMeterArray = {}

def deg(a):
    """
    Returns the argument in degrees
    """
    return float(a)*180/np.pi

def rad(a):
    """
    Returns the argument in radians
    """
    return float(a)/180*np.pi

def rpm(a):
    """
    Returns the argument in rpm
    """
    return float(a)*30/np.pi

def radPerS(a):
    """
    Returns the argument in rad/s
    """
    return float(a)*np.pi/30

def temperature(y,T0):
    return T0-0.0065*y

def pressure(p0,T0,y):
    return p0*(1-0.0065*y/T0)**(g()/(287.058*0.0065))

def density(p0,T0,y):
    return pressure(p0,T0,y)/(287.058*temperature(y,T0))

def g():
    return 9.80665

def p0ISA():
    return 101325

def T0ISA():
    return 288.15

def pSat(T):
    return (6.1078*10**(7.5*(T-273.15)/((T-273.15)+237.3)))*100

def pVapour(T0,y,humidity):
    return humidity*pSat(temperature(y,T0))

def pDry(humidity,p0,T0,y):
    return pressure(p0,T0,y)-pVapour(T0,y,humidity)

def densityWithHumidity(humidity,p0,T0,y):
    return (pressure(p0,T0,y) / (287.058*temperature(y,T0)))*(1 + humidity)/(1 + humidity* 461.495 / 287.058)


def dynamicViscosity(T):
    return 	1.51204129*T**3/2/(T+120)

def kinematicViscocity(my,rho):
    return my/rho

def reynoldsNumber(nu,velocity,c):
    return c*velocity/nu

def setFigureNumber():
    global figureNumber
    figureNumber = figureNumber+1

def getFigureNumber():
    return figureNumber

def getDateString():
    collection= datetime.now()
    myString = collection.isoformat() #collection.year + collection.month +collection.day
    myString = re.sub("\:", "-", myString)
    myString = re.sub("\.", "MS", myString)
    return myString

def writeDict(_dictonary,_fileName):
    f = open(_fileName,'w')
    temp = str(_dictonary)
    temp = re.sub(",", "\n", temp)
    f.write(str(temp))
    f.close()

def readDict(_dictonary,_fileName):
    f = open(_fileName,'r')
    my_dict = eval1(f.read())
    f.close()

def eval1(x):
    try: return eval(x)
    except: return x
    A = { eval1(y[0]) : eval1(y[1]) for y in [x.split() for x in open(filename).readlines()] }





def lengthToPlaneFromPulley(x,y):
    """
    Returns the actual length of the line
    """
    return (x**2+y**2)**0.5

