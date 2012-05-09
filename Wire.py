#-------------------------------------------------------------------------------
# Name:        Wire
# Purpose:     Simulate a Wire under gravity
#
# Author:      Petter TKØ
#
# Created:     07.05.2012
# Copyright:   (c) Petter TKØ 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------

from math import *
from pylab import *
import time

np = 20 # Points on line
m = 1.0/np # Mass of each element on the string
l = 1.0/(np-1) # Original distance between each element
g = 1 # Gravity
x = linspace(0.0,1.0,np) # Position of each point on the wire.
y = zeros_like(x)
v = zeros_like(x) # Velocity y of each point on the wire.
u = zeros_like(x) # Velocity in x direction
Fx   = zeros_like(x) # Forces for each element in x direction.
Fy  = zeros_like(x) # Forces for each element in y direction
T   = [0] # Time
dt  = 1.0/(np*2) # Timestep
Tmax = 20    # Maximal simulation time
k   = 1     # Spring constant
d   = 0.05    # Dampning constant

def dist(P1,P2):
    """
    Returns the distance between P1 and P2
    """
    return ((P1[0]-P2[0])**2+(P1[1]-P2[1])**2)**0.5

def sign(a):
    """
    Returns +1/-1 depending on the sign of a
    """
    if a==0:
        return 1
    else:
        return a/abs(a)

def calcF():
    """
    Calculates the forces acting on each wire element
    """
    for i in range(1,len(x)-1):
        l1t=dist([x[i-1],y[i-1]],[x[i],y[i]])
        l1x=x[i-1]-x[i]
        l1y=y[i-1]-y[i]
        l2t=dist([x[i+1],y[i+1]],[x[i],y[i]])
        l2x=x[i+1]-x[i]
        l2y=y[i+1]-y[i]
        Fy[i] = 0
        Fy[i] += -m*g
        Fy[i] += (1-l/l1t)*k*l1y/l
        Fy[i] += (1-l/l2t)*k*l2y/l
        Fy[i] += -v[i]*d
        Fx[i] = 0
        Fx[i] += (1-l/l1t)*k*l1x/l
        Fx[i] += (1-l/l2t)*k*l2x/l
        Fx[i] += -u[i]*d



def euler():

    for i in range(1,len(x)-1):

        accY        =   Fy[i]/m
        accX        =   Fx[i]/m
        v[i]        =   v[i]+accY*dt
        u[i]        =   u[i]+accX*dt
        y[i]        =   y[i]+v[i]*dt
        x[i]        =   x[i]+u[i]*dt

def simDyn(line):

    n=0
    while T[-1]<Tmax:
        n+=1
        calcF()
        euler()
        T.append(T[-1]+dt)
        if(n==10):
            line.set_ydata(y)
            line.set_xdata(x)
            title("T: %s"%T[-1])
            draw()
            n=0




def main():
    t1=time.time()
    ion()
    line, = plot(x,y)
    ylim(-1,1)
    xlim(0,1)
    draw()
    simDyn(line)
    t2=time.time()
    print t1-t2





if __name__ == '__main__':
    main()