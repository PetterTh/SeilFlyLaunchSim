#-------------------------------------------------------------------------------
# Name:        winch
# Purpose:
#
# Author:      skatun
#
# Created:     20/05/2012
# Copyright:   (c) skatun 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from selFunc import *
import numpy as np

def Swinch(_lineForce,_drumDiameter,_winchStallTorque,_winchZeroSpeed,
        _momentOnWinchDrum,_dt):
    """
    Returns the length of wire which the winch collects during 1 dt
    If the torque is bigger than the stall torque, It is assumed that
    the winch stops and does not reverse
    """
    momentOnDrum = momentOnWinchDrum(_lineForce, _drumDiameter)# Torque acting on the cylinder [Nm]

    omega = radPerS(speedOfWinchDrum(_winchStallTorque,_winchZeroSpeed,
                    momentOnDrum))# Rotational speed of the winch [rad/s]

    return omega*(_drumDiameter/2)*_dt # The amount of line the winch collects [m]
    #lw.append(lw[-1]+S)


def momentOnWinchDrum(_lineForce,_drumDiameter):
    return _lineForce*_drumDiameter/2  # Torque acting on the cylinder [Nm]

def speedOfWinchDrum(_winchStallTorque,_winchZeroSpeed,_momentOnWinchDrum):
    return rpm(min(max(0,(1-_momentOnWinchDrum/_winchStallTorque)),1)
            *_winchZeroSpeed) # Rotational speed of the winch [rpm]

def kLine(E,lineDiameter):
    """
    Returns the spring constant of the line
    As the line is shortened will the springconstant increase
    Assumes the constant is reduced inverse to the relative length
    """
    # Springcoefficient of the line [N] E*pi*d^2/4
    return  E*np.pi*(lineDiameter)**2/4


def lineLength(x,y,L0):
    """
    Returns the actual length of the line
    """
    return (x**2+y**2)**0.5+L0

def deltaLineForce(_lineLength,_deltaLineLength,_kLine):
    return _deltaLineLength/_lineLength*_kLine


def lineForce(_phase,_lineForceOld,_deltaLineForce):
    """
    Returns the force in the line
    dF=(dL+winchspeed)/k
    """
    if _phase >= 4:
        return 0
    else:
        return max(0,(_lineForceOld+_deltaLineForce))


def deltaLineLength(lineLengthOld,x,y,L0):
    return lineLength(x,y,L0)-lineLengthOld

def lineDiameter(phase,D0,l0,deltaL):
    if phase >= 4:
        return D0
    else:
        return np.sqrt(l0/(l0+deltaL))*D0

def diameter():
    global drumDiameter,layersOnDrum
    if lw[-1]>50*n:
        layersOnDrum
        drumDiameter = drumDiameter+lineDiameter

    return drumDiameter