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

def Swinch(lineForce,drumDiameter,winchStallTorque,winchZeroSpeed,
        momentOnWinchDrum):
    """
    Returns the length of wire which the winch collects during 1 dt
    If the torque is bigger than the stall torque, It is assumed that
    the winch stops and does not reverse
    """
    momentOnDrum = momentOnWinchDrum(lineForce,
                drumDiameter)# Torque acting on the cylinder [Nm]

    omega = radPerS(speedOfWinchDrum(winchStallTorque,winchZeroSpeed,
                    momentOnWinchDrum))# Rotational speed of the winch [rad/s]

    return omega*(drumDiameter/2)*dt # The amount of line the winch collects [m]
    #lw.append(lw[-1]+S)
    return S

def momentOnWinchDrum(lineForce,drumDiameter):
    return lineForce*drumDiameter/2  # Torque acting on the cylinder [Nm]

def speedOfWinchDrum(winchStallTorque,winchZeroSpeed,momentOnWinchDrum):
    return rpm(min(max(0,(1-momentOnWinchDrum/winchStallTorque)),1)
            *winchZeroSpeed) # Rotational speed of the winch [rpm]

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

def lineForce(phase,lineForce,drumDiameter,winchStallTorque,
            winchZeroSpeed,momentOnWinchDrum,lineLengthOld,x,y,L0,
            E,lineDiameter):
    """
    Returns the force in the line
    dF=(dL+winchspeed)/k
    """
    if phase == 4:
        return 0
    else:
        return max(0,lineForce+(deltaLineLength(lineLengthOld,x,y,L0)
                +Swinch(lineForce,drumDiameter,winchStallTorque,winchZeroSpeed,
                momentOnWinchDrum))/lineLength(x,y,L0)*kLine(E,lineDiameter))

def deltaLineLength(lineLengthOld,x,y,L0):
    return lineLength(x,y,L0)-lineLengthOld


def diameter():
    global drumDiameter,layersOnDrum
    if lw[-1]>50*n:
        layersOnDrum
        drumDiameter = drumDiameter+lineDiameter

    return drumDiameter