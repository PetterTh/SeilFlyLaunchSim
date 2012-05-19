#-------------------------------------------------------------------------------
# Name: Lnd.py
# Purpose: Calculates all Lift and drag for the SeilSim.py
#
# Author: Petter
#
# Created: 14.05.2012
# Copyright: (c) Petter 2012
# Licence: <your licence>
#-------------------------------------------------------------------------------
import numpy as np
from math import *
from pylab import *

from selFunc import *

# Drag formulas
def calcCd(AR,cdInducedFactor,cl,cdParasiticSpeedFlap,cdParasiticStartFlap,speedFlapPos,startFlapPos,flapPos,Re,refRe,ReCoeff,cdInference):
    """
    Returns the total drag coefficient
    """ 
    cli = cdInduced(AR,cdInducedFactor,cl)
    cdp = cdParasitic(cdParasiticSpeedFlap,cdParasiticStartFlap,speedFlapPos,startFlapPos,flapPos,Re,refRe,ReCoeff)
    cdf = cdInterference(cdInference,Re,refRe,ReCoeff)

    return cli + cdp +cdf

def cdParasitic(cdParasiticSpeedFlap,cdParasiticStartFlap,speedFlapPos,startFlapPos,flapPos,Re,refRe,ReCoeff):
    """
    Returns the parasitic drag coefficient
    """
    #-2.5 0.025 drag coeff speed flap
    #  10 0.04  drag coeff start flap
    y1 = cdParasiticSpeedFlap
    x1 = speedFlapPos
    y2 = cdParasiticStartFlap
    x2 = startFlapPos
    return (y2-y1)/(x2-x1)*(flapPos-x1)+y1

def cdInterference(cdInference,Re,refRe,ReCoeff):
    """
    Returns the interference drag coefficient
    """
    return cdInference*(refRe/Re)**ReCoeff

def cdInduced(AR,cdInducedFactor,cl):
    """
    Returns the induced drag coefficient
    """
    # 0.95 < cdInducedFactor < 1 
    return cl**2/(np.pi*AR*cdInducedFactor)


# Lift formulas
def calcCl(attAng,clAlphaCoeff,speedFlapPos,startFlapPos,speedFlapCl0,startFlapCl0,flapPos):
    """
    Returns the lift coefficient
    """
    return clAlpha(clAlphaCoeff)*rad(attAng)+calcClO(speedFlapPos,startFlapPos,speedFlapCl0,startFlapCl0,flapPos)

def clAlpha(clAlphaCoeff):
    return 2*np.pi*clAlphaCoeff

def calcClO(speedFlapPos,startFlapPos,speedFlapCl0,startFlapCl0,flapPos):
    #-2.5 0.111 stall occurs at 10deg
    #10 0.9 stall occurs at 5deg
    y1 = speedFlapCl0
    x1 = speedFlapPos
    y2 = startFlapCl0
    x2 = startFlapPos
    return (y2-y1)/(x2-x1)*(flapPos-x1)+y1

# These are just for control or future uses
def clAlpha2():
    #0deg 0.278
    #3.5deg 0.615

    y1 = 0.278
    x1 = 0
    y2 = 0.615
    x2 = 3.3
    return (y2-y1)/(x2-x1)*180/np.pi


def clCD():
##    cl[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 ]
##    cdp_n5_Re150000 = [.0127 .0113 .0108 .0108 .0107 0.0106 .0113 .0123 .0153 .0215 .0320]
##    cdp_n0_Re200000 = [.0113 0.0098 .0088 .0085 .0084  .0088 .0098 .0115 .0145  .0185 .00265]
##    cdp_n_25_Re300000 = [.0098 0.0078 .0073 .0078 0.0090 .0113 0.0170 0.0220]

    # ved 5deg flap cwTot~0.049 ved 0 deg 0.031 
    return 3




