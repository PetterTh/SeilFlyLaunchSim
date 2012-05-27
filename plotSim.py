#-------------------------------------------------------------------------------
# Name:plotSim
# Purpose:
#
# Author:      skatun
#
# Created:     20/05/2012
# Copyright:   (c) skatun 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import numpy as np
from pylab import *
from selFunc import *


def initPlot(_logg,_counterPhase):
    global logg,counterPhase,colors

    colors = ["r","b","m","y","black","o"]
    logg = _logg
    counterPhase = _counterPhase

def plotXY(showOn,saveOn):

    grid()
    x= logg['x']
    teller = 1
    for key, value in logg.items():
        #if key!='x'
        figure(teller)
        print key
        for index in range(0,len(counterPhase)-1):
            plot(x[counterPhase[index]:counterPhase[index+1]],value[counterPhase[index]:counterPhase[index+1]],colors[index])
            xlabel("X-Position [m]")
            ylabel(key)
            title(key)
            grid()
        myString = getDateString()
        saveFigMy()
        savefig('Figures/' + myString + '.png')
        teller = teller +1

    #plot(x,y)



    if showOn:
        show()
    if saveOn:
        myString = getDateString()
        saveFigMy()
        savefig('Figures/' + myString + '.png')

def plotXYSingle(key,showOn,saveOn):


    x = logg['x']
    y = logg[key]

    figure(1)
    for index in range(0,len(counterPhase)-1):
        plot(x[counterPhase[index]:counterPhase[index+1]],y[counterPhase[index]:counterPhase[index+1]],colors[index])
        xlabel("X-Position [m]")
        ylabel("Y-Position [m]")
        title('Height plot')
        grid()
    myString = getDateString()
    saveFigMy()
    savefig('Figures/' + myString + '.png')

def plotSingle(x,y,showOn,saveOn,figNum):

    figure(figNum)
    for index in range(0,len(counterPhase)-1):
        plot(x[counterPhase[index]:counterPhase[index+1]],y[counterPhase[index]:counterPhase[index+1]],colors[index])
        xlabel("X-Position [m]")
        ylabel("Y-Position [m]")
        title('Height plot')
        grid()
    myString = getDateString()
    saveFigMy()
    savefig('Figures/' + myString + '.png')

    #plot(x,y)



    if showOn:
        show()
    if saveOn:
        myString = getDateString()
        saveFigMy()
        savefig('Figures/' + myString + '.png')

def plotSingleNoColors(x,y,showOn,saveOn,figNum,key):

    figure(figNum)
    plot(x,y)
    xlabel(key)
    ylabel("Y-Position [m]")
    title('Height different phases')
    legend(['phase 2','phase 3','phase 4'])
    grid()

    if showOn:
        show()
    if saveOn:
        myString = getDateString()
        saveFigMy()
        savefig('Figures/' + myString + '.png')


def plotSim(save,plotOn):
    """
    Plots some graphs providing some information
    """
    colors = ["r","b","m","y","black","o"]

    #Plot for position and force,energy and velocity
    figure(1)
    subplot(2,2,1)
    grid()
    for index in range(0,len(counterPhase)-1):
       plot(x[counterPhase[index]:counterPhase[index+1]],y[counterPhase[index]:counterPhase[index+1]],colors[index])
    xlabel("X-Position [m]")
    ylabel("Y-Position [m]")


    subplot(2,2,2)
    grid()
    for index in range(0,len(counterPhase)-1):
        plot(x[counterPhase[index]:counterPhase[index+1]],lf[counterPhase[index]:counterPhase[index+1]],colors[index])
    xlabel("X-Position [m]")
    ylabel("Force in wire [N]")

    subplot(2,2,3)
    grid()
    for index in range(0,len(counterPhase)-1):
        plot(x[counterPhase[index]:counterPhase[index+1]],velocity[counterPhase[index]:counterPhase[index+1]],colors[index])
    xlabel("X-Position [m]")
    ylabel("Velocity")

    subplot(2,2,4)
    grid()
    for index in range(0,len(counterPhase)-1):
        plot(x[counterPhase[index]:counterPhase[index+1]],E[counterPhase[index]:counterPhase[index+1]],colors[index])
    xlabel("X-Position [m]")
    ylabel("Energy [J]")

    if save:
        savefig('Figures/fig1.png')

# Angles plot
    figure(2)
    subplot(2,2,1)
    grid()
    for index in range(0,len(counterPhase)-1):
        plot(x[counterPhase[index]:counterPhase[index+1]],psiAng[counterPhase[index]:counterPhase[index+1]],colors[index])
    xlabel("X-Position [m]")
    ylabel("Psi angle [deg]")

    subplot(2,2,2)
    grid()
    for index in range(0,len(counterPhase)-1):
        plot(x[counterPhase[index]:counterPhase[index+1]],velAng[counterPhase[index]:counterPhase[index+1]],colors[index])
    xlabel("X-Position [m]")
    ylabel("Velocity angle [deg]")

    subplot(2,2,3)
    grid()
    for index in range(0,len(counterPhase)-1):
        plot(x[counterPhase[index]:counterPhase[index+1]],gamma[counterPhase[index]:counterPhase[index+1]],colors[index])
    xlabel("X-Position [m]")
    ylabel("Gamma angle [deg]")

    subplot(2,2,4)
    grid()
    for index in range(0,len(counterPhase)-1):
        plot(x[counterPhase[index]:counterPhase[index+1]],attAng[counterPhase[index]:counterPhase[index+1]],colors[index])
    xlabel("X-Position [m]")
    ylabel("Angle of attack [deg]")

    if save:
        savefig('Figures/fig2.png')


### Drag coeffcient plots
##    figure(3)
##    subplot(2,2,1)
##    grid()
##    for index in range(0,len(counterPhase)-1):
##        plot(x[counterPhase[index]:counterPhase[index+1]],clVal[counterPhase[index]:counterPhase[index+1]],colors[index])
##    xlabel("X-Position [m]")
##    ylabel("CL coeffcient [-]")
##
##    subplot(2,2,2)
##    grid()
##    for index in range(0,len(counterPhase)-1):
##        plot(x[counterPhase[index]:counterPhase[index+1]],cdVal[counterPhase[index]:counterPhase[index+1]],colors[index])
##    xlabel("X-Position [m]")
##    ylabel("CD coeffcient [-]")
##
##    subplot(2,2,3)
##    grid()
##    for index in range(0,len(counterPhase)-1):
##        plot(x[counterPhase[index]:counterPhase[index+1]],loadFactor[counterPhase[index]:counterPhase[index+1]],colors[index])
##    xlabel("X-Position [m]")
##    ylabel("Load Factor [-]")
##
##    subplot(2,2,4)
##    grid()
##    for index in range(0,len(counterPhase)-1):
##        plot(x[counterPhase[index]:counterPhase[index+1]],cdiVal[counterPhase[index]:counterPhase[index+1]],colors[index])
##    xlabel("X-Position [m]")
##    ylabel("CD induced coeffcient [-]")
##    if save:
##        savefig('Figures/fig3.png')


# Winch configurations  plot
    figure(4)

    subplot(2,2,1)
    grid()
##    print "len x", len(x)
##    print "len M", len(M)
##    print counterPhase[3]
    plot(x[0:len(M)],M,"r")
    xlabel("X-Position [m]")
    ylabel("Moment on winch [Nm]")

    subplot(2,2,2)
    grid()
    plot(x[0:len(M)],rpm)
    xlabel("X-Position [m]")
    ylabel("Speed of winch [rpm]")

    subplot(2,2,3)
    grid()
    plot(x,lf)
    xlabel("X-Position [m]")
    ylabel("Line on winch [m]")

    subplot(2,2,4)
    grid()
    plot(x,y)
    xlabel("X-Position [m]")
    ylabel("CD induced coeffcient [-]")
    if save:
        savefig('Figures/fig4.png')

    if plotOn:
        show()

