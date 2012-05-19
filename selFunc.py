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
