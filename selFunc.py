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
figureNumber = 1
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