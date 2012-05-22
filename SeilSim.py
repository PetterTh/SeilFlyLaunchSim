import numpy as np
from scipy import optimize
from math import *
from pylab import *
from random import *

from Lnd import *
from selFunc import *
from winch import *
from plotSim import *
from plotSensitivity import *


# Some global parameters

def init():
    global figureNumber,alt
    alt = 2
    figureNumber = 1
    """
        1: all phase flown
        2: thrown in the air
        3: no zooming
        4: taking out all the energy of line
    """

    planeParameters0 =          {'wingSpan':3,
                                'wingLoading':3.5,
                                'aspectRatio':15,
                                'refRe':150000,
                                'cdInference':0.01,
                                'cdInducedFactor':0.9,
                                'cdParasiticSpeedFlap':0.025,
                                'cdParasiticStartFlap':0.035,
                                'clAlphaCoeff':0.9,
                                'maxLoadFactor':50,
                                'speedFlapCl0':0.111,
                                'startFlapCl0':0.9,
                                'ReCoeff':0.5}

    flightParameters0 =         {'preTensionOfLine':150,
                                'launchAngle':70,
                                'takeOffSpeed':10,
                                'setpointAOA':8,
                                'startFlapPos':10,
                                'gammaR0':200,
                                'gammaR1':200,
                                'gammaR2':600,
                                'gammaR3':400,
                                'diveStartAngle':75,
                                'speedFlapPos':-2.5,
                                'thermicFlapPos':5,
                                'climbAngle':80};

    lineParameters0 =           {'lineDiameter':1.4e-3,
                                'totalLineLength':400,
                                'EModule':2e9,
                                'lineDragCoeffsient':1.3,
                                'parachuteDragCoeffcient':3,
                                'parachuteArea':0.2};

    winchParameters0 =          {'drumDiameter':50e-3,
                                'drumLength':300e-3,
                                'winchStallTorque':9.8,
                                'winchZeroSpeed':3800,
                                'distanceToPulley':200};

    flighConditionsParameters0 = {'windSpeed':10,
                                'thermic':5,
                                'thermicCeil':600,
                                'temperatureAtGround':25,
                                'pressureAtGround':101325,
                                'humidity':0};

    paraMeterArray0 = planeParameters0.copy()
    paraMeterArray0.update(flightParameters0)
    paraMeterArray0.update(winchParameters0)
    paraMeterArray0.update(lineParameters0)
    paraMeterArray0.update(flighConditionsParameters0)


    #paraMeterArray0 = [planeParameters0,flightParameters0, winchParameters0,lineParameters0,flighConditionsParameters0]


    initSim()

    return paraMeterArray0,planeParameters0,flightParameters0, winchParameters0,lineParameters0,flighConditionsParameters0


def reset():
    """
    Resets all neccecary variables
    """

    #Plane parameters
    global AR,wingArea,pm,refRe,cdInference,cdInducedFactor,cdParasiticSpeedFlap
    global cdParasiticStartFlap,clAlphaCoeff,maxLoadFactor,vMinMy

    #Flight parameters
    global clAlphaCoeff,speedFlapPos,startFlapPos,speedFlapCl0,startFlapCl0,flapPos
    global gamma0,gammaR1,gammaR2,preTensionOfLine,takeOffSpeed,setpointAOA

    #Flight conditions parameters
    global humidity,pressureAtGround,tempAtGround,windSpeed,thermic,thermicCeil

    #Line parameters
    global totalLineLength,EModule

    #Winch parameters
    global distanceToPulley,drumDiameter,wzs,wst,_lineOnWinch


    # need to delete these
    global l,lw,lf,k,lf,gamma,psi,x,y,u,v,T,attAng,velAng,omega,E,phase, velocity,gammaDesiredAngle
    global heightPhase,integral,psiAng,cdiVal,loadFactor,cdVal,clVal,flapPos,M,rpm
    global wingArea,AR,wingSpan,pm,counterPhase
    global cdTotal,clTotal

    global AR,cdInducedFactor,clTotal,cdParasiticSpeedFlap,cdParasiticStartFlap,flapPos,Re,refRe,ReCoeff,cdInference
    global flapPosPhase2 , diveStartAngle,flapPosPhase3 ,climbAngle



    """
    ********* PLANE PARAMETERS *******************************
    """

    AR = float(paraMeterArray['aspectRatio'])                             # Aspect ratio [-]
    wingSpan = float(paraMeterArray['wingSpan'])                          # Wingspan [m]
    wingLoading = float(paraMeterArray['wingLoading'])                    # Wing loading [kg/m^2]
    wingArea  = float(wingSpan**2/AR)                                      # Planform area of the plane [m^2]
    pm  = float(wingLoading*wingArea)                                      # Mass of plane [kg]

    refRe = float(paraMeterArray['refRe'])
    cdInference = float(paraMeterArray['cdInference'])                    # Interfernece cd factor
    cdInducedFactor = float(paraMeterArray['cdInducedFactor'])            # Addition of cd induced factor
    cdParasiticSpeedFlap = float(paraMeterArray['cdParasiticSpeedFlap'])  # value for parasitic drag coeffcient when in speed flap mode
    cdParasiticStartFlap = float(paraMeterArray['cdParasiticStartFlap'])  # value for parasitic drag coeffcient when in speed flap mode
    clAlphaCoeff = float(paraMeterArray['clAlphaCoeff'])                  # Reduction of clAlpha from 2*pi
    maxLoadFactor = float(paraMeterArray['maxLoadFactor'])
    speedFlapCl0  = float(paraMeterArray['speedFlapCl0'])
    startFlapCl0 = float(paraMeterArray['startFlapCl0'])

    # Not implemented...
    ReCoeff = float(paraMeterArray['ReCoeff'])
    Re = float(150000)





    """
    ******************** Flight parameters ********************************
    """
    speedFlapPos = float(paraMeterArray['speedFlapPos'])
    startFlapPos = float(paraMeterArray['startFlapPos'])

    """
    ******** LAUNCHCONFIGURATION IN PHASE 0*****************************
    """
    # Preforce applied to the wire [N]
    preTensionOfLine = float(paraMeterArray['preTensionOfLine'])

    """
    ******** LAUNCHCONFIGURATION IN PHASE 1, in alt2 mode used as init...
    """
    # Takeoff speed [m/s]
    takeOffSpeed = float(paraMeterArray['takeOffSpeed'])

    # Launch angle [deg]
    gamma0 = float(paraMeterArray['launchAngle'])

    # Rate of gamma change [deg/s]. A maximal value of which the gamma
    # can change per second. Used to limit the turn rate
    gammaR0 = float(paraMeterArray['gammaR0'])

    """
    ******** LAUNCHCONFIGURATION IN PHASE 2*****************************
    """
    setpointAOA = float(paraMeterArray['setpointAOA'])
    flapPosPhase2 = float(paraMeterArray['startFlapPos'])


    """
    ******** LAUNCHCONFIGURATION IN PHASE 3*****************************
    """
    # Rate of gamma change [deg/s]. A maximal value of which the gamma
    # can change per second. Used to limit the turn rate
    gammaR1 = float(paraMeterArray['gammaR1'])
    diveStartAngle = float(paraMeterArray['diveStartAngle'])

    flapPosPhase3 = float(paraMeterArray['speedFlapPos'])


    """
    ******** LAUNCHCONFIGURATION IN PHASE 4*****************************
    """
    # Rate of gamma change [deg/s]. A maximal value of which the gamma
    # can change per second. Used to limit the turn rate
    gammaR2 = float(paraMeterArray['gammaR2'])
    flapPosPhase4 = float(paraMeterArray['speedFlapPos'])
    climbAngle = float(paraMeterArray['climbAngle'])


    """
    ******** LAUNCHCONFIGURATION IN PHASE 5 *****************************
    """
    flapPosPhase5 = float(paraMeterArray['thermicFlapPos'])
    vMinMy = float(velocityMin(flapPosPhase5))
    # Rate of gamma change [deg/s]. A maximal value of which the gamma
    # can change per second. Used to limit the turn rate
    gammaR3 = float(paraMeterArray['gammaR3'])

    """
    ********* WINCH PARAMETERS ********************************
    """
    # Winch Stall torque - Torque at zero speed [Nm]
    wst = float(paraMeterArray['winchStallTorque'])

    # Winch zero torque speed  - Speed where the winch has no torque[rad/s]
    wzs = radPerS(float(paraMeterArray['winchZeroSpeed']))

    # Diameter of the cylinger collecting the wire [m]
    drumDiameter=float(paraMeterArray['drumDiameter'])

    # Distance between the winch and the pulley. The plane is assumed to start
    # at the same location as the winch [m]
    distanceToPulley = float(paraMeterArray['distanceToPulley'])

    # Length of winch drum
    drumLength = float(paraMeterArray['drumLength'])

    """
    ********* LINE PARAMETERS ********************************
    """

     # Diameter of the line [m]
    lineDiameter = float(paraMeterArray['lineDiameter'])
    totalLineLength  = float(paraMeterArray['totalLineLength'])

    """
    Some Different E values:
        Steel 210e9
        Rubber 0.01e9-0.1e9
        Nylon 2e9-4e9
    """
    EModule    =  float(paraMeterArray['EModule'])
    lineDragCoeffsient = float(paraMeterArray['lineDragCoeffsient'])
    parachuteDragCoeffcient = float(paraMeterArray['parachuteDragCoeffcient'])
    parachuteArea = float(paraMeterArray['parachuteArea'])


    """
    ****************** WIND PARAMETERS ***********************************
    Not fully implemented yet
    """
    # The headwind meassured at 2 meters height, m/s
    windSpeed = float(paraMeterArray['windSpeed'])

    # The thermic upwind speed component, m/s
    thermic = float(paraMeterArray['thermic'])

    # where the thermic wind ceils, or where it is known,
    # its intrepated linear up to this height
    thermicCeil = float(paraMeterArray['thermicCeil'])

    tempAtGround = float(paraMeterArray['temperatureAtGround'])+273.15
    pressureAtGround = float(paraMeterArray['pressureAtGround'])
    humidity = float(paraMeterArray['humidity'])

    initSim()
    loggingReset()


def initSim():
    """
    Here goes all the other variables that we need
    """
    global Tmax,dt,phase,colors
    global integral,gammaDesiredAngle,previous_error,gammaDesiredAngle0,Kp,Kd,Ki
    Tmax = 50               # Maximal simulation time [s]
    dt  = 0.005              # Time step for the calculation [s]
    phase = 0               # initial phase

    integral  = 0           # used for the I controller
    previous_error = 0


    Kp = 1.1
    Ki = 0
    Kd = 0


    # rewrite all this ........
    gammaDesiredAngle0 = 100  # init for the climbangle
    if alt ==1:
        gammaDesiredAngle = gammaDesiredAngle0 # Launch angle init

    if alt == 2:
        gammaDesiredAngle = 0 #gamma0 # Launch angle init



    colors = ["r","b","m","y","black","o"]

def loggingReset():
    global _u,_v,_x,_y,_ax,_ay
    global _attAng,_velAng,_psi,_gamma
    global _flapPos
    global _lineLengthToPlane,_velocity,_rho,_clTotal,_cdTotal,_fDrag,_fLift
    global _lineForce,_lineDiameter,_kLine,_fx,_fy
    global E,T
    global _drumDiameter,_momentOnWinchDrum,_lineOnWinch
    """
    Each phase of the launch determines how the plane should behave:
        0: Preload the wire. The plane is stationary and the winch starts to tention the wire
        1: Takeoff. The plane is released and starts to accelerate along the ground. This phase starts when the line reaches the wanted preforce
        2: Liftoff and climb. The plane increases the angle of attack and leaves the ground. This phase starts when the takeoffspeed (v0) is reached
        3: Dive. At the peak height the plane starts to dive against the pulley to increase its energy
        4: Climb. The winch is released and the plane starts to climb.
        5: Climb is ended and plane is flying straight ahead with thermic setting.

    """

    # Logging arrays.......
    _flapPos = [0]
    _attAng = [0]       # Angle of attack [deg]
    _velAng = [gamma0]  # The planes velocity angle [deg]
    layersOnDrum = 1        # Layers of line on drum


    _clTotal = [calcCl(_attAng[-1],clAlphaCoeff,speedFlapPos,startFlapPos,speedFlapCl0,startFlapCl0,_flapPos[-1])]
    _cdTotal = [calcCd(AR,cdInducedFactor,_clTotal[-1],cdParasiticSpeedFlap,cdParasiticStartFlap,speedFlapPos,startFlapPos,_flapPos[-1],Re,refRe,ReCoeff,cdInference)]

    if alt==1:
        _u   = [0.0]      # Plane velocity in x direction [m/s]
        _v   = [0.0]      # Plane velocity in y direction [m/s]
        _gamma = [0]   # Angle between plane and ground [deg]

    elif alt == 2:
        _u   = [np.cos(rad(gamma0))*takeOffSpeed ]      # Plane velocity in x direction [m/s]
        _v   = [np.sin(rad(gamma0))*takeOffSpeed ]      # Plane velocity in y direction [m/s]
        _gamma = [gamma0]   # Angle between plane and ground [deg]

    #Not implemented
    else:
        _u   = [0.0]      # Plane velocity in x direction [m/s]
        _v   = [0.0]      # Plane velocity in y direction [m/s]
        _gamma = [0]   # Angle between plane and ground [deg]

    _velocity = [0] # Apperent velocity

    _x   = [-distanceToPulley]    # Positon of the plane in x direction [m]
    _y   = [0.0]        # Position of the plane in y direction [m]

    _ax = [0.0]
    _ay = [0.0]

    # Some global variables
    _lineLengthToPlane    = [totalLineLength]      # Length of the line between the winch and the plane [m]
    _lineOnWinch    = [0]        # Meters of line on the winch [m]
    _lineForce  = [0]          # Lineforce [N]


    _psi = [0.0]          # Angle between line and ground [deg]


    T  = [0.0]         # Accumulated time [s]

    omega = [0]        # Speed of the winch [rad/s]
    E = [0]            # Total energy of the plane [J]

    rpm = [0]          # speed of winch[rpm]
    _drumDiameter = [drumDiameter]

    # moment on winch cylinder [Nm]
    _momentOnWinchDrum = [momentOnWinchDrum(_lineForce[-1],_drumDiameter[-1])]

    phase = 0
    integral = 0

    flapPos = flapPosPhase2

    cdiVal = [0] #induced drag coefficient
    clVal  = [0] #lift coeffcient
    cdVal = [0] #drag coefficient
    loadFactor = [0] # loadfactor

    _lineDiameter = [lineDiameter()]
    _kLine = [kLine(EModule,_lineDiameter[-1])]


    """
    *********** STANDARD CONFIGURATION *********************************
    """
    # Airdensity [kg/m3]
    _rho = [densityWithHumidity(humidity,pressureAtGround,tempAtGround,_y[-1])]

    """
    Here we add some logging features ...
    """

    _fDrag = [Fdrag(_cdTotal[-1],_velocity[-1],_rho[-1],wingArea)]
    _fLift = [Flift(_clTotal[-1],_velocity[-1],_rho[-1],wingArea)]

    _fx = [0]
    _fy = [0]


def calcVelAng():
    """
    Returns the angle of the velocity vector of the plane
    """
    return deg(math.atan2(_v[-1],_u[-1]))

def  calcLoadFactor():
    r = (x[-1]**2+y[-1]**2)**0.5
    zentAcc = velocity[-1]**2/r
    return (Flift(velocity[-1])+loadVector *zentAcc)/g()*pm

def calcAttAng():
    """
    Returns the angle of attack for the plane
    """
    if _u[-1]==0 and _v[-1]==0: # If the plane is stationary the attack angle should be zero
        return 0


    return _gamma[-1]-_velAng[-1]

def calcGamma():
    global _flapPos, loadVector
    """
    Returns the plane angle.
    Assumes the plane flies with gamma0 degrees towards the line all the time
    """
    if phase <=3:
        gammaMyR = gammaR1 # radius of top of zoom
    else:
        gammaMyR = gammaR2 # radius of bottom of zoom


    goal = 0
    if phase == 0:
        if alt==2:
            goal=gamma0+setpointAOA # Included to simplify things for the governor
        else:
            goal=gamma0
    if phase == 1:
        goal=gamma0
    if phase == 2:
        goal=gammaDesired()-_psi[-1]
    if phase == 3:
       goal= -_psi[-1]
    if phase == 4:
        goal=climbAngle
        print "climbing"
    if goal < _gamma[-1]:
        loadVector = 1
        #flapPos = 0
        return max(goal,_gamma[-1]-gammaMyR*dt)
    elif goal > _gamma[-1]:
        #flapPos = 0
        loadVector = -1
        return min(goal,_gamma[-1]+gammaMyR*dt)

    else:
        return _gamma[-1]



def calcPsi():
    """
    Returns the line angle

    """

    return deg(math.atan2(_y[-1],(-_x[-1])))

def gammaDesired():
    global integral,gammaDesiredAngle,previous_error

    error = setpointAOA-calcAttAng()
    integral = integral + (error*dt)
    derivative = (error - previous_error)/dt

    gammaDesiredAngle =  gammaDesiredAngle +( Kp*error+ Ki*integral + Kd*derivative)

    previous_error = error
    if gammaDesiredAngle>95:
        gammaDesiredAngle = 95
    if gammaDesiredAngle<50:
        gammaDesiredAngle = 50


    #print "gammaDesiredAngle: ",gammaDesiredAngle ,"Error:", error
    return gammaDesiredAngle



def sumForces():
    """
    Calculates the resulting forces acting on the plane
    """
    global _lineLengthToPlane,_velocity,_rho,_clTotal,_cdTotal,_fDrag,_fLift
    global _lineForce,_fx,_fy,_kLine,_lineDiameter,_lineOnWinch

    _velocity.append( sqrt((_u[-1]+ windSpeed)**2+ (_v[-1])**2))
    _rho.append(densityWithHumidity(humidity,pressureAtGround,tempAtGround,_y[-1]))

    _psi.append(calcPsi())
    _velAng.append(calcVelAng())
    _gamma.append(calcGamma())
    _attAng.append(calcAttAng())
    _clTotal.append(calcCl(_attAng[-1],clAlphaCoeff,speedFlapPos,startFlapPos,speedFlapCl0,startFlapCl0,_flapPos[-1]))
    _cdTotal.append(calcCd(AR,cdInducedFactor,_clTotal[-1],cdParasiticSpeedFlap,cdParasiticStartFlap,speedFlapPos,startFlapPos,_flapPos[-1],Re,refRe,ReCoeff,cdInference))


    if phase>0:
        gravityForce = g()*pm
        _fDrag.append(Fdrag(_cdTotal[-1],_velocity[-1],_rho[-1],wingArea))
        _fLift.append(Flift(_clTotal[-1],_velocity[-1],_rho[-1],wingArea))
    else:
        _fDrag.append(0)
        _fLift.append(0)
        gravityForce = 0


    _lineLengthToPlane.append(lineLength(_x[-1],_y[-1],distanceToPulley))
    _lineDiameter.append(lineDiameter())
    _kLine.append(kLine(EModule,_lineDiameter[-1]))
    _momentOnWinchDrum.append(momentOnWinchDrum(_lineForce[-1],_drumDiameter[-1]))
    s = Swinch(_lineForce[-1],_drumDiameter[-1],wst,wzs,_momentOnWinchDrum[-1],dt)
    _lineOnWinch.append(_lineOnWinch[-1]+s)

    deltaLineLengthMy = _lineLengthToPlane[-1]-_lineLengthToPlane[-2] + s
    deltaLineForceMy = deltaLineForce(_lineLengthToPlane[-1],deltaLineLengthMy,_kLine[-1])

    _lineForce.append(lineForce(phase,_lineForce[-1],deltaLineForceMy))



    _fx.append(-_fDrag[-1]*np.cos(rad(_velAng[-1]))+_lineForce[-1]*np.cos(rad(_psi[-1]))-_fLift[-1]*np.sin(rad(_velAng[-1])))
    _fy.append(-_fDrag[-1]*np.sin(rad(_velAng[-1]))-_lineForce[-1]*np.sin(rad(_psi[-1]))+_fLift[-1]*np.cos(rad(_velAng[-1]))- gravityForce)





def euler():
    """
    Calculates the new position of the plane using forward Euler iteration
    Does not allow the plane to go below y=0 as this is ground level
    """
    global _u,_v,_ax,_ay
    if phase == 0:
        sumForces()
        _u.append(0)
        _v.append(0)
        _ax.append(0)
        _ay.append(0)
    elif phase ==1 and alt ==1:
        sumForces()
        _ax.append(_fx[-1]/pm)
        _ay.append(_fy[-1]/pm)
        _u.append(_u[-1]+_ax[-1]*dt)
        _v.append(0)

    else:
        sumForces()
        _ax.append(_fx[-1]/pm)
        _ay.append(_fy[-1]/pm)
        _u.append(_u[-1]+_ax[-1]*dt)
        _v.append(_v[-1]+_ay[-1]*dt)

    _x.append(_x[-1]+_u[-1]*dt)
    _y.append(max(_y[-1]+_v[-1]*dt,0)) #To make sure y always is positive




def runInit(plotOn):
    global planeParameters
    planeParameters = planeParameters0.copy()
    flightParameters = flightParameters0.copy()
    height = simulate([0])
    plotSim(1,plotOn)
    return height

def simulate(inp):
    global _flapPos,dt
    global T,E
    global AR,cdInducedFactor,clTotal,cdParasiticSpeedFlap,cdParasiticStartFlap,flapPos,Re,refRe,ReCoeff,cdInference
    global phase,counterPhase
    """
    Runs the simulation
    """
    reset()
    teller = 0
    heightPhase = [0.0]
    counterPhase = [0]
    while T[-1]<=Tmax and _y[-1] >= -10.0 and phase<5:



##        clVal.append(calcCl())
##        cdVal.append(calcCd())
##        cdiVal.append(cdInduced())
        #loadFactor.append(calcLoadFactor())

        euler()
        T.append(T[-1]+dt)
        E.append(_y[-1]*g()*pm+0.5*pm*_velocity[-1]**2)

        #print T[-1],attAng,gamma

        # Change phases
        if _lineForce[-1]>preTensionOfLine and phase==0:
            if alt==2: # Alt 2 does not contain any takeoff along the ground
               phase=2
            else:
                phase = 1

            heightPhase.append(_y[-1])
            counterPhase.append(teller)
            print "Phase 1: T:",T[-1],"X:",_x[-1]
        if _velocity[-1]>takeOffSpeed and phase==1:
            phase = 2
            heightPhase.append(_y[-1])
            counterPhase.append(teller)
            print "Phase 2: T:",T[-1],"X:",_x[-1]
        if _psi[-1] > diveStartAngle and phase ==2:
            phase = 3
            #dt = dt/5
            _flapPos[-1] = flapPosPhase3
            heightPhase.append(_y[-1])
            counterPhase.append(teller)
            print "Phase 3: T:",T[-1],"X:",_x[-1]

        if (_y[-1]<50 or _lineForce[-1]==0) and phase ==3:
            phase = 4
            heightPhase.append(_y[-1])
            counterPhase.append(teller)
            if _y[-1]<50:
                print "diveHeight reached"
            if _lineForce[-1]==0:
                print "all lineForce used"
            print "Phase 4: T:",T[-1],"X:",_x[-1]

        if _velocity[-1]<vMinMy and phase>3:
            phase = 5
            heightPhase.append(_y[-1])
            counterPhase.append(teller)
            _flapPos[-1]= flapPosPhase5
            print "Phase 5: T:",T[-1],"X:",_x[-1]

        teller=teller + 1

        if phase==4:
            print "x:",_x[-1],"y:",_y[-1]

    heightPhase.append(_y[-1])
    return heightPhase[2:-1]

##
##    if max(E)==0:
##        return 1000
##    else:
##        return 1/max(E)




def sensitivityCheck(changeArray0):
    global paraMeterArray
##    paraMeterArray = changeArray0.copy()


    keys   = []
    minVal = []
    maxVal = []
    minRes = []
    maxRes = []

    for key, value in changeArray0.items():
        paraMeterArray[key] = value*.9
        minVal.append(paraMeterArray[key])
        minRes.append(simulate(paraMeterArray))

        paraMeterArray[key] = value*1.1
        maxVal.append(paraMeterArray[key])
        maxRes.append(simulate(paraMeterArray))
        keys.append(key)

    return keys,minVal,maxVal,minRes,maxRes

def sensitivityDelta(keys,minVal,maxVal,minRes,maxRes):

    maxValues = max(max(minRes),max(maxRes))

    senstivityArrayTemp= (abs(subtract(minRes,maxRes))/maxValues)
    summOfArray = senstivityArrayTemp.sum(axis=0)

    return senstivityArrayTemp/summOfArray

def printSensitivity(keys,minVal,maxVal,minRes,maxRes):

    sentivity = sensitivityDelta(keys,minVal,maxVal,minRes,maxRes)*100
    print "Parameter    Height2 Height3 Height5 Min Val Max Val "

    for i in range(0,len(sentivity)):
        print keys[i],sentivity[i],"%",minVal[i],maxVal[i]




def simulateTemp(inp):

    return [float(randint(130,160)),float(randint(100,130)),float(randint(200,260))]

def sensitivity():
    global paraMeterArray
    paraMeterArray,planeParameters,flightParameters,winchParameters,lineParameters,flighConditionsParameters = init()
    reset()
    loggingReset()
    testArray = paraMeterArray.copy()
    tempArray = {'wingSpan':3,
                'wingLoading':3.5}
    keys,minVal,maxVal,minRes,maxRes = sensitivityCheck(testArray)
    sensitivityMy = sensitivityDelta(keys,minVal,maxVal,minRes,maxRes)


    plotSensitivity(sensitivityMy,keys,1,1,1)
    printSensitivity(keys,minVal,maxVal,minRes,maxRes)

if __name__=="__main__":

    #lim = ([0,300],[0,50])
    #res=optimize.brute(simulate,lim,Ns=4)
##    planeParameters = planeParameters0
##    refHeight = simulate([0])
##    plotSim(1)
    print "Start!!!!"
##    planeParameters,flightParameters,winchParameters,lineParameters,flighConditionsParameters = init()
##    reset(planeParameters,flightParameters,winchParameters,lineParameters,flighConditionsParameters)
##    loggingReset()
##    print simulate([3])
##    initPlot(_x,_y,counterPhase)
##    plotXY(1,1)



    sensitivity()
    #sensitivity(runInit(0))
    print "Done!!!!"



