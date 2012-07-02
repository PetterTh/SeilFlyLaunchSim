import numpy as np
from scipy import optimize
from math import *
from pylab import *
from random import *

from globalHolder import *
from plotSelector import *

from Lnd import *
from selFunc import *
from winch import *



paraMeterArray = {}
master = -1
exclusive = 0
alt = 2
saveOn= 1
showOn = 1
setPlotVer(9)
"""
    0: All available plot methods
    1: x,y plot

    2: x,plotkey
    3: x, all pararmeters in logg array
    4: x, plotKeyArray

    5: plotKey, height
    6: all pararmeters in paraMeterArray, height
    7: plotKeyArray,height

    8: all sensitivity plot
    9: planeParameters sensitivity plot
    10: flightParameters sensitivity plot
    11: winchParameters sensitivity plot
    12: lineParameters sensitivity plot
    13: flighConditionsParameters sensitivity plot
    14: plotKeyArray sensitivity plot
    15: x,y for different values of plotKey
    16: x,y plot for different alternatives(alt)
    17: all sensitivity

"""

setPlotKey('energy')
#plotKeyArray = ['lineDiameter','totalLineLength']
#plotKeyArray = ['lineDiameter','deltaLineLength','totalLineLength']
#plotKeyArray = ['distanceToPulley','totalLineLength']
plotKeyArray = ['distanceToPulley']
scalePlot = 1
numbersOfVariationStepForPlots = 1
counterPhase = []
debug = 0
# Some global parameters

def getPlotSettings():
    test = getPlotVer()
    bbfh = getPlotKey()
    setArr = {  'master':master,
                'counterPhase':getCounterPhases(),
                'saveOn':saveOn,
                'showOn':showOn,
                'scalePlot':scalePlot,
                'numbersOfVariationStepForPlots':numbersOfVariationStepForPlots,
                'plotKey':getPlotKey(),
                'plotKeyArray':plotKeyArray,
                'plotVer':getPlotVer(),
                'exclusive':exclusive}
    return setArr



def setAlt(type):
    global alt
    alt = type

def init():
    global figureNumber,alt,debug
    global planeParameters0,flightParameters0, winchParameters0,lineParameters0,flighConditionsParameters0
    global paraMeterArray


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
                                'gammaR2':550,
                                'gammaR3':550,
                                'diveStartAngle':75,
                                'speedFlapPos':-2.5,
                                'thermicFlapPos':5,
                                'climbAngle':80};

    lineParameters0 =           {'lineDiameter':1.4e-3,
                                'totalLineLength':400,
                                'EModule':2e9,
                                'lineDragCoeffsient':0.69,
                                'parachuteDragCoeffcient':3,
                                'parachuteArea':0.2};

    winchParameters0 =          {'drumDiameter':50e-3,
                                'drumLength':300e-3,
                                'winchStallTorque':9.8,
                                'winchZeroSpeed':3800,
                                'distanceToPulley':200};

    flighConditionsParameters0 = {'windSpeed':0,
                                'thermic':0,
                                'thermicCeil':600,
                                'temperatureAtGround':25,
                                'pressureAtGround':101325,
                                'humidity':0};

    paraMeterArray0 = planeParameters0.copy()
    paraMeterArray0.update(flightParameters0)
    paraMeterArray0.update(winchParameters0)
    paraMeterArray0.update(lineParameters0)
    paraMeterArray0.update(flighConditionsParameters0)


    paraMeterArray = paraMeterArray0.copy()
    reset()
    initSim()

    return paraMeterArray0


def reset():
    """
    Resets all neccecary variables
    """

    #Plane parameters
    global AR,wingArea,pm,refRe,cdInference,cdInducedFactor,cdParasiticSpeedFlap
    global cdParasiticStartFlap,clAlphaCoeff,maxLoadFactor,vMinMy

    #Flight parameters
    global clAlphaCoeff,speedFlapPos,startFlapPos,speedFlapCl0,startFlapCl0
    global gammaR0,gammaR1,gammaR2,gammaR3,preTensionOfLine,takeOffSpeed,setpointAOA
    global flapPosPhase,gamma0,diveStartAngle,climbAngle

    #Flight conditions parameters
    global humidity,pressureAtGround,tempAtGround,windSpeed,thermic,thermicCeil

    #Line parameters
    global totalLineLength,EModule,lineDiameterMy

    #Winch parameters
    global distanceToPulley,drumDiameter,wzs,wst,_lineOnWinch,drumLength


    # need to delete these
    global Re,refRe,ReCoeff


    setParametersArray(paraMeterArray)

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
    flapPosPhase0 = 0.0

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

    flapPosPhase1 = float(paraMeterArray['thermicFlapPos'])

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
    lineDiameterMy = float(paraMeterArray['lineDiameter'])
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

    flapPosPhase = [flapPosPhase0,flapPosPhase1,flapPosPhase2,flapPosPhase3,flapPosPhase4,flapPosPhase5]

    initSim()
    loggingReset()


def initSim():
    """
    Here goes all the other variables that we need
    """
    global Tmax,dt,phase,colors
    global integral,gammaDesiredAngle,previous_error,gammaDesiredAngle0,Kp,Kd,Ki
    Tmax = 100               # Maximal simulation time [s]
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
        gammaDesiredAngle = 0#gammaDesiredAngle0 # Launch angle init

    if alt >= 2:
        gammaDesiredAngle = gamma0 # Launch angle init


    layersOnDrum = 1        # Layers of line on drum

    colors = ["r","b","m","y","black","o"]

def loggingReset():
    global _u,_v,_x,_y,_ax,_ay
    global _attAng,_velAng,_psi,_gamma
    global _flapPos
    global _totalLineLength,_velocity,_rho,_clTotal,_cdTotal,_fDrag,_fLift
    global _lineForce,_lineDiameter,_kLine,_fx,_fy
    global E,T
    global _drumDiameter,_momentOnWinchDrum,_lineOnWinch,_lengthToPlaneFromPulley
    global _lineArea,_lineDrag,_fTotalDrag,_deltaLineLength

    # All arrays can be zero since we start with preloading the wire in time step 0,
    # except of x,k,lineDiameter,drumDiameter,total Line Length which is set to init valiues

    # Logging arrays.......
    _flapPos = [0.0]
    _clTotal = [0.0]
    _cdTotal = [0.0]

    _velocity = [0.0] # Apperent velocity
    _x   = [-distanceToPulley]# Positon of the plane in x direction [m]
    _y   = [0.0]# Position of the plane in y direction [m]
    _ax = [0.0]# Acceleration of the plane in x direction [m]
    _ay = [0.0]# Acceleration of the plane in y direction [m]
    _u   = [0.0]      # Plane velocity in x direction [m/s]
    _v   = [0.0]      # Plane velocity in y direction [m/s]
    _fx = [0.0]
    _fy = [0.0]
    _totalLineLength    = [totalLineLength] # Length of the line between the winch and the plane [m]
    _deltaLineLength =[0.0]
    _lineOnWinch    = [0.0]# Meters of line on the winch [m]
    _lineForce  = [0.0] # Lineforce [N]
    _drumDiameter = [drumDiameter]

    # angles
    _psi = [0.0]          # Angle between line and ground [deg]
    _gamma = [0.0]   # Angle between plane and ground [deg]
    _attAng = [0.0]       # Angle of attack [deg]
    _velAng = [0.0]       # The planes velocity angle [deg]
    T  = [0.0]         # Accumulated time [s]
    omega = [0.0]        # Speed of the winch [rad/s]
    E = [0.0]            # Total energy of the plane [J]
    rpm = [0.0]          # speed of winch[rpm]
    _momentOnWinchDrum = [0] # moment on winch cylinder [Nm]
    cdiVal = [0.0] #induced drag coefficient
    clVal  = [0.0] #lift coeffcient
    cdVal = [0.0] #drag coefficient
    loadFactor = [0.0] # loadfactor
    _lineDiameter = [lineDiameterMy]
    _kLine = [kLine(EModule,_lineDiameter[-1])]
    _rho = [densityWithHumidity(humidity,pressureAtGround,tempAtGround,_y[-1])]# Airdensity [kg/m3]
    _fDrag = [0.0]
    _fLift = [0.0]
    _lengthToPlaneFromPulley = [distanceToPulley]
    _lineArea = [lineArea(_lineDiameter[-1],_lengthToPlaneFromPulley[-1],_gamma[-1])]
    _lineDrag = [0.0]
    _fTotalDrag = [0.0]


def getPlaneParameters():
    return planeParameters0
def getFlightParameters():
    return flightParameters0

def getWinchParameters():
    return winchParameters0
def getLineParameters():
    return lineParameters0
def getFlighConditionsParameters():
    return flighConditionsParameters0

#########################################################################################################

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
    if phase == 2:
        gammaMyR = gammaR0 # radius of bottom of launch
    elif phase==3:
        gammaMyR = gammaR1 # radius of top of zoom
    elif phase==4:
        gammaMyR = gammaR2 # radius of bottom of zoom
    else:
        gammaMyR = gammaR3 # radius of flatten


    goal = 0
    if phase == 0:
        if alt>=2:
            goal=gamma0+setpointAOA # Included to simplify things for the governor
        else:
            goal=0#gamma0
    if phase == 1:
        goal=0#gamma0
    if phase == 2:
        goal=gammaDesired()-_psi[-1]
    if phase == 3:
        if alt==3:
            goal = climbAngle
        else:
            goal= -_psi[-1]
    if phase == 4:
        goal=climbAngle
    if phase == 5:
        goal=0
    if goal < _gamma[-1]:
        loadVector = 1
        #flapPos = 0
        return  max(goal,_gamma[-1]-gammaMyR*dt)
    elif goal > _gamma[-1]:
        #flapPos = 0
        loadVector = -1
##        temp = min(goal,_gamma[-1]-gammaMyR*dt)
##        if phase==5:
##            if temp<0:
##                return 0
##            else:
##                return temp

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
    global _totalLineLength,_velocity,_rho,_clTotal,_cdTotal,_fDrag,_fLift
    global _lineForce,_fx,_fy,_kLine,_lineDiameter,_lineOnWinch
    global _lengthToPlaneFromPulley,_lineArea,_lineDrag,_fTotalDrag,_deltaLineLength

    velX = abs(_u[-1]+ windSpeed)*(_u[-1]+ windSpeed)
    velY = abs(_v[-1])*_v[-1]

    vel2 = ( sqrt(velX+ velY))

    _velocity.append( sqrt((_u[-1]+ windSpeed)**2 + (_v[-1])**2))

    _rho.append(densityWithHumidity(humidity,pressureAtGround,tempAtGround,_y[-1]))

    _psi.append(calcPsi())
    _velAng.append(calcVelAng())
    _gamma.append(calcGamma())
    _attAng.append(calcAttAng())


    _clTotal.append(calcCl(_attAng[-1],clAlphaCoeff,speedFlapPos,startFlapPos,speedFlapCl0,startFlapCl0,_flapPos[-1]))
    _cdTotal.append(calcCd(AR,cdInducedFactor,_clTotal[-1],cdParasiticSpeedFlap,cdParasiticStartFlap,speedFlapPos,startFlapPos,_flapPos[-1],Re,refRe,ReCoeff,cdInference))
    cdLineMy = cdLine()

    if onLine:
        _totalLineLength.append(lineLength(_x[-1],_y[-1],distanceToPulley))
        _drumDiameter.append(drumDiameter2(_lineOnWinch[-1],_drumDiameter[-1],drumLength))
        _momentOnWinchDrum.append(momentOnWinchDrum(_lineForce[-1],_drumDiameter[-1]))
        s = Swinch(_lineForce[-1],_drumDiameter[-1],wst,wzs,_momentOnWinchDrum[-1],dt)
        _lineOnWinch.append(_lineOnWinch[-1]+s)

        _deltaLineLength.append(_totalLineLength[-1]-_totalLineLength[-2] + s)

        _lineDiameter.append(lineDiameter(phase,lineDiameterMy,_totalLineLength[-1],_deltaLineLength[-1]))
        _kLine.append(kLine(EModule,_lineDiameter[-1]))
        deltaLineForceMy = deltaLineForce(_totalLineLength[-1],_deltaLineLength[-1],_kLine[-1])

        _lineForce.append(lineForce(phase,_lineForce[-1],deltaLineForceMy))
        _drumDiameter.append(_drumDiameter[-1])


        _lengthToPlaneFromPulley.append(lengthToPlaneFromPulley(_x[-1],_y[-1]))
        _lineArea.append(lineArea(_lineDiameter[-1],_lengthToPlaneFromPulley[-1],_gamma[-1]))
    else:
        _totalLineLength.append(_totalLineLength[-1])
        _momentOnWinchDrum.append(0.0)
        _lineOnWinch.append(_lineOnWinch[-1])
        _deltaLineLength.append(0.0)
        _lineDiameter.append(lineDiameterMy)
        _kLine.append(_kLine[-1])
        _lineForce.append(_lineForce[-1])
        _drumDiameter.append(_drumDiameter[-1])
        _lengthToPlaneFromPulley.append(lengthToPlaneFromPulley(_x[-1],_y[-1]))
        _lineArea.append(_lineArea[-1])

    if phase>0:
        gravityForce = g()*pm

##        _fDrag.append(Fdrag(_cdTotal[-1],_velocity[-1],_rho[-1],wingArea))
##        _fLift.append(Flift(_clTotal[-1],_velocity[-1],_rho[-1],wingArea))
##        _lineDrag.append(FlineDrag(cdLineMy,_velocity[-1],_rho[-1],_lineArea[-1],phase))
        _fDrag.append(Fdrag(_cdTotal[-1],vel2,_rho[-1],wingArea))
        _fLift.append(Flift(_clTotal[-1],vel2,_rho[-1],wingArea))
        _lineDrag.append(FlineDrag(cdLineMy,vel2,_rho[-1],_lineArea[-1],phase))
    else:
        _fDrag.append(0)
        _fLift.append(0)
        _lineDrag.append(0)
        gravityForce = 0

    _fTotalDrag.append(_fDrag[-1] + _lineDrag[-1])
    _fx.append(-_fTotalDrag[-1]*np.cos(rad(_velAng[-1]))+_lineForce[-1]*np.cos(rad(_psi[-1]))-_fLift[-1]*np.sin(rad(_velAng[-1])))
    _fy.append(-_fTotalDrag[-1]*np.sin(rad(_velAng[-1]))-_lineForce[-1]*np.sin(rad(_psi[-1]))+_fLift[-1]*np.cos(rad(_velAng[-1]))- gravityForce)




##        clVal.append(calcCl())
##        cdVal.append(calcCd())
##        cdiVal.append(cdInduced())
        #loadFactor.append(calcLoadFactor())



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
    global phase,heightPhase,onLine,paraMeterArray
    global _gamma,_u,_v,change

    """
    Runs the simulation
    """
    paraMeterArray = inp
    reset()
    teller = 0
    heightPhase = []
    counterPhase = []
    change = 0
    onLine = 1
    i = 1
    while T[-1]<=Tmax and _y[-1] >= -10.0 and phase<6:

        _flapPos.append(flapPosPhase[phase])
        euler()
        T.append(T[-1]+dt)
        E.append(_y[-1]*g()*pm+0.5*pm*_velocity[-1]**2)


        """
        Each phase of the launch determines how the plane should behave:
            0: Preload the wire. The plane is stationary and the winch starts to tention the wire
            1: Takeoff. The plane is released and starts to accelerate along the ground. This phase starts when the line reaches the wanted preforce
            2: Liftoff and climb. The plane increases the angle of attack and leaves the ground. This phase starts when the takeoffspeed (v0) is reached
            3: Dive. At the peak height the plane starts to dive against the pulley to increase its energy
            4: Climb. The winch is released and the plane starts to climb.
            5: Climb is ended and plane is flying straight ahead with thermic setting.

        """
        if phase>2 and i:
            dt = dt * 5
            i = 0
        # Change phases
        if _lineForce[-1]>preTensionOfLine and phase==0:
            if alt>=2: # Alt 2 does not contain any takeoff along the ground
               phase=2
               _u   = [np.cos(rad(gamma0))*takeOffSpeed ]      # Plane velocity in x direction [m/s]
               _v   = [np.sin(rad(gamma0))*takeOffSpeed ]      # Plane velocity in y direction [m/s]
               _gamma = [gamma0]   # Angle between plane and ground [deg]
               change = 1
            else:
                phase = 1
                change = 1

        if _velocity[-1]>takeOffSpeed and phase==1:
            phase = 2
            change = 1

        if _psi[-1] > diveStartAngle and phase ==2:
            phase = 3
            change = 1
            #dt = dt/5


        if (_y[-1]<50 or _lineForce[-1]==0) and phase ==3:
            phase = 4
            change = 1
##            if _y[-1]<50:
##                print "diveHeight reached"
##            if _lineForce[-1]==0:
##                print "all lineForce used" ,heightPhase[-1]-_y[-1]
            onLine = 0

        if _velocity[-1]<vMinMy*1.1 and phase>=3 and phase< 5:
            phase = 5
            change = 1
##            print _velocity[-1]-vMinMy


        if (_gamma[-1]<=0 and phase>=5):
            phase = 6
            change = 1


        teller=teller + 1



        if change:
            heightPhase.append(_y[-1])
            counterPhase.append(teller)
            change = 0

            if debug:
                print "Phase ", phase, ": T:",T[-1],"X:",_x[-1],"Y:",_y[-1]

    setCounterPhase(counterPhase)

    # This code snippets makes sure we set all height of all phases except of phase0,
    # in case of failed launch, i.e line tear etc the returned array will have eqaul heights for the last elemnts.
    returnArray = []
    for i in range(0,4):
        if i>=len(heightPhase)-1:
            returnArray.append(returnArray[-1])
        else:
            returnArray.append(heightPhase[i+1])
    return returnArray


def simulateTemp(inp):

    return [float(randint(130,160)),float(randint(100,130)),float(randint(200,260))]




def saveLogg():
    loggArray = {'u':_u,
                'v':_v,
                'x':_x,
                'y':_y,
                'ax':_ax,
                'ay':_ay,
                'attAng':_attAng,
                'velAng':_velAng,
                'psi':_psi,
                'gamma':_gamma,
                'flapPos':_flapPos,
                'totalLineLength':_totalLineLength,
                'velocity':_velocity,
                'rho':_rho,
                'clTotal':_clTotal,
                'cdTotal':_cdTotal,
                'fDrag':_fDrag,
                'fLift':_fLift,
                'lineForce':_lineForce,
                'lineDiameter':_lineDiameter,
                'kLine':_kLine,
                'fx':_fx,
                'fy':_fy,
                'time':T,
                'energy':E,
                'drumDiameter':_drumDiameter,
                'momentOnWinchDrum':_momentOnWinchDrum,
                'lineOnWinch':_lineOnWinch,
                'deltaLineLength':_deltaLineLength}



    return loggArray


def getParametersArray():
    return paraMeterArray

def setParametersArray(paraMeterArray0):
    global paraMeterArray
    paraMeterArray = paraMeterArray0

def masterTables():
    tableSelector(0)
    tableSelector(1)
    tableSelector(2)
    tableSelector(3)
    tableSelector(4)
    tableSelector(5)
    tableSelector(6)


def masterPlot():

    setPlotVer(2)
    setPlotKey('energy')
    plotSelector()

    setPlotVer(16)
    plotSelector()

    setPlotVer(1)
    plotSelector()

    setPlotVer(8)
    plotSelector()

    setPlotVer(9)
    plotSelector()

    setPlotVer(10)
    plotSelector()

    setPlotVer(11)
    plotSelector()

    setPlotVer(12)
    plotSelector()

    setPlotVer(13)
    plotSelector()


if __name__=="__main__":

    #lim = ([0,300],[0,50])
    #res=optimize.brute(simulate,lim,Ns=4)

    print "Start!!!!"
    if master==1:
        masterPlot()
    elif master == -1:
        masterTables()
    else:
        #wplotSelector()
        tableSelector()


    print "Done!!!!"



