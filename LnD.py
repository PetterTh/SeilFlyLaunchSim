#-------------------------------------------------------------------------------
# Name:        Lnd.py
# Purpose:     Calculates all Lift and drag for the SeilSim.py
#
# Author:      Petter
#
# Created:     14.05.2012
# Copyright:   (c) Petter 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import numpy as np
from selFunc import *

def calcCd(phase,cd0_2,cd0_3,cl,AR):
    """
    Returns the drag coefficient
    """
    if phase <3:
        cd0 =cd0_2
    else:
        cd0 = cd0_3

    cli = cdInduced(cl,AR)
    return cli + cd0

def clO(flapPos):
    #-2.5 0.111 stall occurs at 10deg
    #10 0.9 stall occurs at 5deg
    y1 = 0.111
    x1 = -2.5
    y2 = 0.9
    x2 = 10
    return (y2-y1)/(x2-x1)*(flapPos-x1)+y1


def clAlpha():
    #0deg 0.278
    #3.5deg 0.615

    y1 = 0.278
    x1 = 0
    y2 = 0.615
    x2 = 3.3
    return (y2-y1)/(x2-x1)*180/np.pi


def clCD():
    cl=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0 ]
    cdp_n5_Re150000 = [.0127,.0113,.0108,.0108,.0107,0.0106,.0113,.0123,.0153,.0215,.0320]
    #cdp_n0_Re200000 = [.0113 0.0098 .0088 .0085 .0084  .0088 .0098 .0115 .0145  .0185 .00265]
    #cdp_n_25_Re300000 = [.0098 0.0078 .0073 .0078 0.0090 .0113 0.0170 0.0220]

    # ved 5deg flap cwTot~0.049 ved 0 deg 0.031
    return 3


def calcCl(phase,cl0_2,cl0_3,alpha):
    """
    Returns the lift coefficient
    """
    if phase <3:
        cl0 =cl0_2
    else:
        cl0 = cl0_3
    return 2*np.pi*rad(alpha)+cl0



def cdInterference(Re):
    return 0.01*(150000/Re)^.5



def cdInduced(cl,AR):
    return cl**2/np.pi/AR


