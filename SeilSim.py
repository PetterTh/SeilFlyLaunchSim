import numpy as np
from scipy import optimize
from math import *
from pylab import *




# Some global parameters
alt = 2
"""
    1: all phase flown
    2: thrown in the air
    3: no zooming
"""
"""
********* PLANE PARAMETERS *******************************
"""
A  = 0.6                # Planform area of the plane [m2]
pm  = 2                 # Mass of plane [kg]
AR = 13                 # Aspect ratio
wingSpan = 3            # Wingspan

"""
********* WINCH PARAMETERS ********************************
"""
wst = 9.8               # Winch Stall torque - Torque at zero speed [Nm]
wzs = 3800/60*2*np.pi   # Winch zero torque speed  - Speed where the winch has no torque[rad/s]
drumDiameter=0.055                 # Diameter of the cylinger collecting the wire [m]
l0 = 200                # Distance between the winch and the pulley. The plane is assumed to start at the same location as the winch [m]
drumLength = 0.3        # Length of winch drum
lineDiameter = 1.3*10**-3 # Diameter of the line
k0    =  2*10**9*np.pi*(lineDiameter)**2/4       # Springcoefficient of the line [N] E*pi*d^2/4
layersOnDrum = 1        # Layers of line on drum

"""
Some Different E values:
    Steel 210e9
    Rubber 0.01e9-0.1e9
    Nylon 2e9-4e9
"""


"""
****************** WIND PARAMETERS ***********************************
Not implemented yet
"""

windSpeed = 3       # The headwind meassured at 2 meters height, m/s
thermic = 5         # The thermic upwind speed component, m/s
thermicCeil = 600   # where the thermic wind ceils, or where it is known, its intrepated linear up to this height

"""
******** LAUNCHCONFIGURATION IN PHASE 0*****************************
"""
pf =100                 # Preforce applied to the wire [N]
"""
******** LAUNCHCONFIGURATION IN PHASE 1, in alt2 mode used as init...
"""
v0 = 10                 # Takeoff speed [m/s]
gamma0 = 75             # Launch angle
"""
******** LAUNCHCONFIGURATION IN PHASE 2*****************************
"""
cl0_2 = 0.1
cd0_2 = 0.003
setpointAOA = 8         # AOA during climb phase
integral  = 0           # used for the I controller
previous_error = 0
gammaDesiredAngle0 = 80  # init for the climbangle
Kp = 1
Ki = 0.1
Kd = 0.0
"""
******** LAUNCHCONFIGURATION IN PHASE 3*****************************
"""
gammaR3 = 200 # Rate of gamma change [deg/s]. A maximal value of which the gamma can change per second. Used to limit the turn rate
diveStartAngle = 75
cl0_3 = 0.05
cd0_3 = 0.002
"""
******** LAUNCHCONFIGURATION IN PHASE 4*****************************
"""
gammaR4 = 600 # Rate of gamma change [deg/s]. A maximal value of which the gamma can change per second. Used to limit the turn rate
climbAngle = 75
cl0_4 = 0.05
cd0_4 = 0.002

"""
******** LAUNCHCONFIGURATION IN PHASE 5*****************************
"""
vMin = 9
gammaR5 = 200 # Rate of gamma change [deg/s]. A maximal value of which the gamma can change per second. Used to limit the turn rate
finalHeight = 0

"""
*********** STANDARD CONFIGURATION *********************************
"""
g  = 9.81               # Gavitational acceleration [m/s2]
rho = 1.4               # Airdensity [kg/m3]
Tmax = 50               # Maximal simulation time [s]
dt  = 0.01              # Time step for the calculation [s]
phase = 0               # initial phase
"""
Each phase of the launch determines how the plane should behave:
    0: Preload the wire. The plane is stationary and the winch starts to tention the wire
    1: Takeoff. The plane is released and starts to accelerate along the ground. This phase starts when the line reaches the wanted preforce
    2: Liftoff and climb. The plane increases the angle of attack and leaves the ground. This phase starts when the takeoffspeed (v0) is reached
    3: Dive. At the peak height the plane starts to dive against the pulley to increase its energy
    4: Climb. The winch is released and the plane starts to climb.

"""
# Some global variables
cl = 0.0     # Lift Coefficient [-]
cd = 0.0     # Drag coefficient [-]


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

def reset():
    """
    Resets all neccecary variables
    """
    global l,lw,lf,k,lf,gamma,psi,x,y,u,v,T,attAng,velAng,omega,E,phase, velocity,gammaDesiredAngle,integral

    # Some global variables
    l    = [2*l0]      # Length of the line between the winch and the plane [m]
    lw    = [0]        # Meters of line on the winch [m]
    lf  = [0]          # Lineforce [N]
    k    = k0          # Actual spring force of the line [N/m]
    gamma = [gamma0]   # Angle between plane and ground [deg]
    psi = 0.0          # Angle between line and ground [deg]

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
    phase = 0
    integral = 0



def calcCd():
    """
    Returns the drag coefficient
    """
    if phase <3:
        cd0 =cd0_2
    else:
        cd0 = cd0_3

    return cl**2/AR + cd0

def calcCl():
    """
    Returns the lift coefficient
    """
    if phase <3:
        cl0 =cl0_2
    else:
        cl0 = cl0_3
    return 2*np.pi*rad(attAng[-1])+cl0

def calcVelAng():
    """
    Returns the angle of the velocity vector of the plane
    """
    return deg(math.atan2(v[-1],u[-1]))


def calcAttAng():
    """
    Returns the angle of attack for the plane
    """
    if u[-1]==0 and v[-1]==0: # If the plane is stationary the attack angle should be zero
        return 0
    return gamma[-1]-velAng[-1]

def calcGamma():
    """
    Returns the plane angle.
    Assumes the plane flies with gamma0 degrees towards the line all the time
    """
    if phase <=3:
        gammaMyR = gammaR3 # Does not work!!!!!!!
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

        return max(goal,gamma[-1]-gammaMyR*dt)
    elif goal > gamma[-1]:
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
    gammaDesiredAngle =  gammaDesiredAngle+( Kp*error+ Ki*integral + Kd*derivative)

    previous_error = error
    #if gammaDesiredAngle>90:
    #    gammaDesiredAngle = 90
    #if gammaDesiredAngle<75:
    #    gammaDesiredAngle = 75

    #print "gammaDesiredAngle: ",gammaDesiredAngle
    return gammaDesiredAngle

def Flift(vel):
    """
    Returns the lift force of the plane based on the input velocity
    Flift = cl*v^2*rho/2*A
    """

    return cl*vel**2*rho/2*A


def Fdrag(vel):
    """
    Returns the grad force of the plane based on the input velocity
    Flift = cd*v^2*rho/2*A
    """

    return cd*vel**2*rho/2*A

def diameter():
    global drumDiameter,layersOnDrum
    if lw[-1]>50*n:
        layersOnDrum
        drumDiameter = drumDiameter+lineDiameter

    return drumDiameter



def Swinch():
    """
    Returns the length of wire which the winch collects during 1 dt
    If the torque is bigger than the stall torque, It is assumed that the winch stops and does not reverse
    """
    M=lf[-1]*(drumDiameter/2) # Torque acting on the cylinder [Nm]
    omega.append(min(max(0,(1-M/wst)),1)*wzs) # Rotational speed of the winch [rad/s]
    S = omega[-1]*(drumDiameter/2)*dt # The amount of line the winch collects [m]
    lw.append(lw[-1]+S)
    return S

def kLine():
    """
    Returns the spring constant of the line
    As the line is shortened will the springconstant increase
    Assumes the constant is reduced inverse to the relative length
    """

    return  k0#*l[0]/(l[0]-lw[-1])

def lLine():
    """
    Returns the actual length of the line
    """
    return ((x[-1])**2+y[-1]**2)**0.5+l0

def fLine():
    """
    Returns the force in the line
    dF=(dL+winchspeed)/k
    """
    if phase == 4:
        return 0
    else:
        return max(0,lf[-1]+((l[-1]-l[-2])+Swinch())/l[-1]*kLine())


def sumForces():
    """
    Calculates the resulting forces acting on the plane
    returns the fx and fy
    """
    vel = sqrt(np.power(u[-1],2)+np.power(v[-1],2))


    fx=-Fdrag(vel)*np.cos(rad(velAng[-1]))+lf[-1]*np.cos(rad(psi))-Flift(vel)*np.sin(rad(velAng[-1]))
    fy=-Fdrag(vel)*np.sin(rad(velAng[-1]))-lf[-1]*np.sin(rad(psi))+Flift(vel)*np.cos(rad(velAng[-1]))-g*pm

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




def simulate(inp):
    global gamma,psi,wf,velAng,attAng,cd,cl,phase,pf,v0,vMIn
    """
    Runs the simulation
    """
    reset()
    while T[-1]<=Tmax and y[-1] >= -10.0 and phase<5:

        psi = calcPsi()
        velAng.append(calcVelAng())
        gamma.append(calcGamma())
        attAng.append(calcAttAng())
        cl=calcCl()
        cd=calcCd()
        l.append(lLine())
        lf.append(fLine())
        Euler()
        T.append(T[-1]+dt)
        E.append(y[-1]*g*pm+0.5*pm*(u[-1]**2+v[-1]**2))
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
            #print "Phase 1: T:",T[-1],"X:",x[-1]
        if (u[-1]**2+v[-1]**2)**0.5>v0 and phase==1:
            phase = 2
            #print "Phase 2: T:",T[-1],"X:",x[-1]
        if psi > diveStartAngle and phase ==2:
            phase = 3
            #print "Phase 3: T:",T[-1],"X:",x[-1]
        if (y[-1]<50 or lf[-1]==0) and phase ==3:
            phase = 4
            #print "Phase 4: T:",T[-1],"X:",x[-1]

        if (u[-1]**2+v[-1]**2)**0.5<vMin and phase>3:
            phase = 5

    print pf,max(E)
    if max(E)==0:
        return 1000
    else:
        return 1/max(E)




def plotSim():
    """
    Plots some graphs providing some information
    """
    subplot(2,2,1)
    xlabel("X-Position [m]")
    ylabel("Y-Position [m]")
    plot(x,y)
    #axes().set_aspect('equal', 'datalim')

    subplot(2,2,2)
    plot(T,lf)
    xlabel("Time [s]")
    ylabel("Force in wire [N]")

    subplot(2,2,4)
    plot(T,attAng)
    xlabel("Time [s]")
    ylabel("Angle of attack [deg]")

    subplot(2,2,3)
    plot(x,velocity)
    xlabel("X-Position [m]")
    ylabel("Velocity")
    show()
    """
    subplot(2,2,3)
    plot(x,E)
    xlabel("X-Position [m]")
    ylabel("Energy [J]")
    """
    show()


if __name__=="__main__":
    #lim = ([0,300],[0,50])
    #res=optimize.brute(simulate,lim,Ns=4)
    simulate([0])
    print "Hmax: ",max(y),"EnergyMax: ",max(E)
    plotSim()

