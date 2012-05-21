import numpy as np
from scipy import optimize
from math import *
from pylab import *
from random import *

from Lnd import *
from selFunc import *




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




    return planeParameters0,flightParameters0, winchParameters0,lineParameters0,flighConditionsParameters0


def reset(planeParameters,flightParameters,winchParameters,lineParameters,flighConditionsParameters):
    """
    Resets all neccecary variables
    """
    #Plane parameters
    global AR,wingArea,pm,refRe,cdInference,cdInducedFactor,cdParasiticSpeedFlap
    global cdParasiticStartFlap,clAlphaCoeff,maxLoadFactor

    #Flight parameters
    global clAlphaCoeff,speedFlapPos,startFlapPos,speedFlapCl0,startFlapCl0,flapPos
    global gamma0

    #Flight conditions parameters
    global humidity,pressureAtGround,tempAtGround,windSpeed,thermic,thermicCeil

    #Line parameters
    global totalLineLength

    global l,lw,lf,k,lf,gamma,psi,x,y,u,v,T,attAng,velAng,omega,E,phase, velocity,gammaDesiredAngle
    global heightPhase,integral,psiAng,cdiVal,loadFactor,cdVal,clVal,flapPos,M,rpm
    global wingArea,AR,wingSpan,pm,counterPhase
    global cdTotal,clTotal

    global AR,cdInducedFactor,clTotal,cdParasiticSpeedFlap,cdParasiticStartFlap,flapPos,Re,refRe,ReCoeff,cdInference
    global flapPosPhase2 , diveStartAngle,flapPosPhase3 ,climbAngle



    """
    ********* PLANE PARAMETERS *******************************
    """

    AR = float(planeParameters['aspectRatio'])                             # Aspect ratio [-]
    wingSpan = float(planeParameters['wingSpan'])                          # Wingspan [m]
    wingLoading = float(planeParameters['wingLoading'])                    # Wing loading [kg/m^2]
    wingArea  = float(wingSpan**2/AR)                                      # Planform area of the plane [m^2]
    pm  = float(wingLoading*wingArea)                                      # Mass of plane [kg]

    refRe = float(planeParameters['refRe'])
    cdInference = float(planeParameters['cdInference'])                    # Interfernece cd factor
    cdInducedFactor = float(planeParameters['cdInducedFactor'])            # Addition of cd induced factor
    cdParasiticSpeedFlap = float(planeParameters['cdParasiticSpeedFlap'])  # value for parasitic drag coeffcient when in speed flap mode
    cdParasiticStartFlap = float(planeParameters['cdParasiticStartFlap'])  # value for parasitic drag coeffcient when in speed flap mode
    clAlphaCoeff = float(planeParameters['clAlphaCoeff'])                  # Reduction of clAlpha from 2*pi
    maxLoadFactor = float(planeParameters['maxLoadFactor'])

    # Not implemented...
    ReCoeff = float(planeParameters['ReCoeff'])
    Re = float(150000)


    speedFlapCl0  = float(planeParameters['speedFlapCl0'])
    startFlapCl0 = float(planeParameters['startFlapCl0'])


    """
    ******************** Flight parameters ********************************
    """
    speedFlapPos = float(flightParameters['speedFlapPos'])
    startFlapPos = float(flightParameters['startFlapPos'])

    """
    ******** LAUNCHCONFIGURATION IN PHASE 0*****************************
    """
    # Preforce applied to the wire [N]
    pf = float(flightParameters['preTensionOfLine'])

    """
    ******** LAUNCHCONFIGURATION IN PHASE 1, in alt2 mode used as init...
    """
    # Takeoff speed [m/s]
    v0 = float(flightParameters['takeOffSpeed'])

    # Launch angle [deg]
    gamma0 = float(flightParameters['launchAngle'])

    # Rate of gamma change [deg/s]. A maximal value of which the gamma
    # can change per second. Used to limit the turn rate
    gammaR0 = float(flightParameters['gammaR0'])

    """
    ******** LAUNCHCONFIGURATION IN PHASE 2*****************************
    """
    setpointAOA = float(flightParameters['setpointAOA'])
    flapPosPhase2 = float(flightParameters['startFlapPos'])


    """
    ******** LAUNCHCONFIGURATION IN PHASE 3*****************************
    """
    # Rate of gamma change [deg/s]. A maximal value of which the gamma
    # can change per second. Used to limit the turn rate
    gammaR1 = float(flightParameters['gammaR1'])
    diveStartAngle = float(flightParameters['diveStartAngle'])

    flapPosPhase3 = float(flightParameters['speedFlapPos'])


    """
    ******** LAUNCHCONFIGURATION IN PHASE 4*****************************
    """
    # Rate of gamma change [deg/s]. A maximal value of which the gamma
    # can change per second. Used to limit the turn rate
    gammaR2 = float(flightParameters['gammaR2'])
    flapPosPhase4 = float(flightParameters['speedFlapPos'])
    climbAngle = float(flightParameters['climbAngle'])


    """
    ******** LAUNCHCONFIGURATION IN PHASE 5 *****************************
    """
    flapPosPhase5 = float(flightParameters['thermicFlapPos'])
    vMin = velocityMin(flapPosPhase5)
    # Rate of gamma change [deg/s]. A maximal value of which the gamma
    # can change per second. Used to limit the turn rate
    gammaR3 = float(flightParameters['gammaR3'])

    """
    ********* WINCH PARAMETERS ********************************
    """
    # Winch Stall torque - Torque at zero speed [Nm]
    wst = float(winchParameters['winchStallTorque'])

    # Winch zero torque speed  - Speed where the winch has no torque[rad/s]
    wzs = radPerS(float(winchParameters['winchZeroSpeed']))

    # Diameter of the cylinger collecting the wire [m]
    drumDiameter=float(winchParameters['drumDiameter'])

    # Distance between the winch and the pulley. The plane is assumed to start
    # at the same location as the winch [m]
    l0 = float(winchParameters['distanceToPulley'])

    # Length of winch drum
    drumLength = float(winchParameters['drumLength'])

    """
    ********* LINE PARAMETERS ********************************
    """

     # Diameter of the line [m]
    lineDiameter = float(lineParameters['lineDiameter'])
    totalLineLength  = float(lineParameters['totalLineLength'])

    """
    Some Different E values:
        Steel 210e9
        Rubber 0.01e9-0.1e9
        Nylon 2e9-4e9
    """
    EModule    =  float(lineParameters['EModule'])
    lineDragCoeffsient = float(lineParameters['lineDragCoeffsient'])
    parachuteDragCoeffcient = float(lineParameters['parachuteDragCoeffcient'])
    parachuteArea = float(lineParameters['parachuteArea'])


    """
    ****************** WIND PARAMETERS ***********************************
    Not fully implemented yet
    """
    # The headwind meassured at 2 meters height, m/s
    windSpeed = float(flighConditionsParameters['windSpeed'])

    # The thermic upwind speed component, m/s
    thermic = float(flighConditionsParameters['thermic'])

    # where the thermic wind ceils, or where it is known,
    # its intrepated linear up to this height
    thermicCeil = float(flighConditionsParameters['thermicCeil'])

    tempAtGround = float(flighConditionsParameters['temperatureAtGround'])+273.15
    pressureAtGround = float(flighConditionsParameters['pressureAtGround'])
    humidity = float(flighConditionsParameters['humidity'])


def loggingReset():
    """
    Each phase of the launch determines how the plane should behave:
        0: Preload the wire. The plane is stationary and the winch starts to tention the wire
        1: Takeoff. The plane is released and starts to accelerate along the ground. This phase starts when the line reaches the wanted preforce
        2: Liftoff and climb. The plane increases the angle of attack and leaves the ground. This phase starts when the takeoffspeed (v0) is reached
        3: Dive. At the peak height the plane starts to dive against the pulley to increase its energy
        4: Climb. The winch is released and the plane starts to climb.

    """






    # Some global variables

    heightPhase = [0]
    counterPhase =[0]
    counterPhase = [0]
    integral  = 0           # used for the I controller
    previous_error = 0
    gammaDesiredAngle0 = 100  # init for the climbangle
    Kp = 1.1
    Ki = 0
    Kd = 0

    flapPos = float(0)
    attAng = float(0)
    layersOnDrum = 1        # Layers of line on drum


    clTotal = [calcCl(attAng,clAlphaCoeff,speedFlapPos,startFlapPos,speedFlapCl0,startFlapCl0,flapPos)]
    cdTotal = [calcCd(AR,cdInducedFactor,clTotal[-1],cdParasiticSpeedFlap,cdParasiticStartFlap,speedFlapPos,startFlapPos,flapPos,Re,refRe,ReCoeff,cdInference)]




    # Some global variables
    l    = [totalLineLength]      # Length of the line between the winch and the plane [m]
    lw    = [0]        # Meters of line on the winch [m]
    lf  = [0]          # Lineforce [N]

    gamma = [gamma0]   # Angle between plane and ground [deg]
    psi = 0.0          # Angle between line and ground [deg]
    psiAng = [0]       # Angle between line and ground [deg]



    x   = [-l[0]/2]    # Positon of the plane in x direction [m]
    y   = [0.0]        # Position of the plane in y direction [m]

    if alt ==1:
        gammaDesiredAngle = gammaDesiredAngle0 # Launch angle init

    if alt == 2:
        gammaDesiredAngle = gamma0 # Launch angle init

    u   = [0]      # Plane velocity in x direction [m/s]
    v   = [0]      # Plane velocity in y direction [m/s]

    velocity = [0]     # Plane total velocity
    T  = [0.0]         # Accumulated time [s]
    attAng = [0]       # Angle of attack [deg]
    velAng = [gamma0]  # The planes velocity angle [deg]
    omega = [0]        # Speed of the winch [rad/s]
    E = [0]            # Total energy of the plane [J]
    M = [0]            # moment on winch cylinder [Nm]
    rpm = [0]          # speed of winch[rpm]
    phase = 0
    integral = 0

    flapPos = flapPosPhase2

    cdiVal = [0] #induced drag coefficient
    clVal  = [0] #lift coeffcient
    cdVal = [0] #drag coefficient
    loadFactor = [0] # loadfactor


    """
    *********** STANDARD CONFIGURATION *********************************
    """
    # Airdensity [kg/m3]
    rho = densityWithHumidity(humidity,pressureAtGround,tempAtGround,y[-1])
    Tmax = 50               # Maximal simulation time [s]
    dt  = 0.005              # Time step for the calculation [s]
    phase = 0               # initial phase

    """
    Here we add some logging features ...
    """
    heightPhase[phase] = y[-1]
    counterPhase[phase]= 0

def calcVelAng():
    """
    Returns the angle of the velocity vector of the plane
    """
    return deg(math.atan2(v[-1],u[-1]))

def  calcLoadFactor():
    r = (x[-1]**2+y[-1]**2)**0.5
    zentAcc = velocity[-1]**2/r
    return (Flift(velocity[-1])+loadVector *zentAcc)/g()*pm

def calcAttAng():
    """
    Returns the angle of attack for the plane
    """
    if u[-1]==0 and v[-1]==0: # If the plane is stationary the attack angle should be zero
        return 0


    return gamma[-1]-velAng[-1]

def calcGamma():
    global flapPos, loadVector
    """
    Returns the plane angle.
    Assumes the plane flies with gamma0 degrees towards the line all the time
    """
    if phase <=3:
        gammaMyR = gammaR3 # radius of top of zoom
    else:
        gammaMyR = gammaR4 # radius of bottom of zoom


    goal = 0
    if phase == 0:
        if alt==2:
            goal=gamma0+setpointAOA # Included to simplify things for the governor
        else:
            goal=gamma0
    if phase == 1:
        goal=gamma0
    if phase == 2:
        goal=gammaDesired()-psi
    if phase == 3:
       goal= -psi
    if phase == 4:
        goal=climbAngle
    if goal < gamma[-1]:
        loadVector = 1
        #flapPos = 0
        return max(goal,gamma[-1]-gammaMyR*dt)
    elif goal > gamma[-1]:
        #flapPos = 0
        loadVector = -1
        return min(goal,gamma[-1]+gammaMyR*dt)

    else:
        return gamma[-1]



def calcPsi():
    """
    Returns the line angle

    """

    return deg(math.atan2(y[-1],(-x[-1])))

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
    returns the fx and fy
    """
    velocity = sqrt((u[-1]+ windSpeed)**2+ (v[-1])**2)

    clTotal = [calcCl(attAng,clAlphaCoeff,speedFlapPos,startFlapPos,speedFlapCl0,startFlapCl0,flapPos)]
    cdTotal = [calcCd(AR,cdInducedFactor,clTotal[-1],cdParasiticSpeedFlap,cdParasiticStartFlap,speedFlapPos,startFlapPos,flapPos,Re,refRe,ReCoeff,cdInference)]

    fx=-Fdrag(cdTotal,velocity,rho,wingArea)*np.cos(rad(velAng[-1]))+lf[-1]*np.cos(rad(psi))-Flift(clTotal,velocity,rho,wingArea)*np.sin(rad(velAng[-1]))
    fy=-Fdrag(cdTotal,velocity,rho,wingArea)*np.sin(rad(velAng[-1]))-lf[-1]*np.sin(rad(psi))+Flift(clTotal,velocity,rho,wingArea)*np.cos(rad(velAng[-1]))-g()*pm

    return fx,fy


def Euler():
    """
    Calculates the new position of the plane using forward Euler iteration
    Does not allow the plane to go below y=0 as this is ground level
    """

    if phase == 0:
        u.append(0)
        v.append(0)
    elif phase ==1 and alt ==1:
        fx,fy = sumForces()
        ax=fx/pm
        ay=fy/pm
        u.append(u[-1]+ax*dt)
        v.append(0)

    else:
        fx,fy = sumForces()
        ax=fx/pm
        ay=fy/pm
        u.append(u[-1]+ax*dt)
        v.append(v[-1]+ay*dt)

    x.append(x[-1]+u[-1]*dt)
    y.append(max(y[-1]+v[-1]*dt,0))


def runInit(plotOn):
    global planeParameters
    planeParameters = planeParameters0.copy()
    flightParameters = flightParameters0.copy()
    height = simulate([0])
    plotSim(1,plotOn)
    return height

def simulate(inp):
    global gamma,psi,wf,velAng,attAng,cd,cl,phase,pf,v0,vMin,heightPhase,counterPhase,psiAng,flapPos,dt
    global clAlphaCoeff,speedFlapPos,startFlapPos,speedFlapCl0,startFlapCl0,flapPos
    global AR,cdInducedFactor,clTotal,cdParasiticSpeedFlap,cdParasiticStartFlap,flapPos,Re,refRe,ReCoeff,cdInference
    """
    Runs the simulation
    """
    reset()
    teller = 0
    while T[-1]<=Tmax and y[-1] >= -10.0 and phase<5:

        psi = calcPsi()

        psiAng.append(calcPsi())
        velAng.append(calcVelAng())
        gamma.append(calcGamma())
        attAng.append(calcAttAng())
        clTotal.append(calcCl(attAng[-1],clAlphaCoeff,speedFlapPos,startFlapPos,speedFlapCl0,startFlapCl0,flapPos))
        cdTotal.append(calcCd(AR,cdInducedFactor,clTotal[-1],cdParasiticSpeedFlap,cdParasiticStartFlap,speedFlapPos,startFlapPos,flapPos,Re,refRe,ReCoeff,cdInference))

##        clVal.append(calcCl())
##        cdVal.append(calcCd())
##        cdiVal.append(cdInduced())
        loadFactor.append(calcLoadFactor())


        l.append(lLine())
        lf.append(fLine())
        Euler()
        T.append(T[-1]+dt)
        E.append(y[-1]*g()*pm+0.5*pm*(u[-1]**2+v[-1]**2))
        velocity.append((u[-1]**2+v[-1]**2)**0.5)
        #print T[-1],attAng,gamma

        # Change phases
        if lf[-1]>pf and phase==0:
            if alt==2: # Alt 2 does not contain any takeoff along the ground
               phase=2
               u[-1]   = np.cos(rad(gammaDesiredAngle))*v0      # Plane velocity in x direction [m/s]
               v[-1]   = np.sin(rad(gammaDesiredAngle))*v0      # Plane velocity in y direction [m/s]
            else:
                phase = 1

            heightPhase.append(y[-1])
            counterPhase.append(teller)
            #print "Phase 1: T:",T[-1],"X:",x[-1]
        if (u[-1]**2+v[-1]**2)**0.5>v0 and phase==1:
            phase = 2
            heightPhase.append(y[-1])
            counterPhase.append(teller)
            #print "Phase 2: T:",T[-1],"X:",x[-1]
        if psi > diveStartAngle and phase ==2:
            phase = 3
            #dt = dt/5
            flapPos= flapPosPhase3
            heightPhase.append(y[-1])
            counterPhase.append(teller)
            #print "Phase 3: T:",T[-1],"X:",x[-1]

        if (y[-1]<50 or lf[-1]==0) and phase ==3:
            phase = 4
            heightPhase.append(y[-1])
            counterPhase.append(teller)
            #print "Phase 4: T:",T[-1],"X:",x[-1]

        if (u[-1]**2+v[-1]**2)**0.5<vMin and phase>3:
            phase = 5
            heightPhase.append(y[-1])
            counterPhase.append(teller)
            flapPos= flapPosPhase5

        teller=teller + 1

##
##    for index, item in enumerate(heightPhase):
##        print "Phase ",index, "Height :",item
##
##
##    print "Pre tension:",pf,"Max energy:",max(E)
##
##
##    print "clAlpha:", clAlpha()
    return heightPhase[-1]

##
##    if max(E)==0:
##        return 1000
##    else:
##        return 1/max(E)





def sensitivityPlane(refHeight):
    global planeParameters
    planeParameters = planeParameters0.copy()

    #refHeight = simulate([0])
    print refHeight
    svarMin = [refHeight]
    svarMax = [refHeight]
    senstivityArrayPlane = [0]
    keysPlane = ["Ref"]



    for key, value in planeParameters0.items():
        value = planeParameters0[key]
        refValue = value
        planeParameters[key] = value*.9
        minVal = planeParameters[key]
        minRes = simulate([0])

        value = planeParameters0[key]
        planeParameters[key] = value*1.1
        maxVal = planeParameters[key]
        maxRes = simulate([0])

        svarMin.append(minRes)
        svarMax.append(maxRes)
        keysPlane.append(key)

        senstivityArrayPlane.append(abs(svarMax[-1]-svarMin[-1])/refHeight)

        if abs(svarMin[-1]-svarMax[-1])<0.01:
            print key, " mostly liked to be not impmented"
        if abs(abs(svarMin[-1]-refHeight) - abs(svarMax[-1]-refHeight))>0.5:
            print key, " is not linear, asymmetric"

        print key, "Ref:",refValue," -- ", refHeight, " Min: ",minVal," -- " ,minRes,"Max: ",maxVal," -- ", maxRes

    totalSum = sum(senstivityArrayPlane)
    for i in range(1,len(senstivityArrayPlane)):
        prosent = senstivityArrayPlane[i]/totalSum*100
        print keysPlane[i], prosent

def sensitivityCheck(changeArray0):

    changeArray = changeArray0.copy()


    keys   = []
    minVal = []
    maxVal = []
    minRes = []
    maxRes = []

    for key, value in changeArray0.items():
        changeArray[key] = value*.9
        minVal.append(changeArray[key])
        minRes.append(simulateTemp(changeArray))

        changeArray[key] = value*1.1
        maxVal.append(changeArray[key])
        maxRes.append(simulateTemp(changeArray))
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

def sensitivity(refHeight):
    global planeParameters0,planeParameters
    sensitivity2(refHeight,planeParameters,planeParameters0)
    sensitivityPlane(refHeight)
    sensitivityFlight(refHeight)

if __name__=="__main__":

    #lim = ([0,300],[0,50])
    #res=optimize.brute(simulate,lim,Ns=4)
##    planeParameters = planeParameters0
##    refHeight = simulate([0])
##    plotSim(1)
    print "Start!!!!"
    sensitivity(runInit(0))
    print "Done!!!!"



