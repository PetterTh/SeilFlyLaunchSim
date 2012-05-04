import numpy as np
from math import *
from pylab import *


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

# Some global parameters

A  = 0.5     # Planform area of the plane [m2]
pm  = 1      # Mass of plane [kg]
k0    =  2*10**9*np.pi*0.005**2/4       # Springcoefficient of the line [N] E*pi*d^2/4
"""
Some Different E values:
    Steel 210e9
    Rubber 0.01e9-0.1e9
    Nylon 2e9-4e9
"""
g  = 9.81    # Gavitational acceleration [m/s2]
rho = 1.4    # Airdensity [kg/m3]
v0 = 10        # Launch speed [m/s]
gamma0 = 80 # Launch angle
Tmax = 150    # Maximal simulation time [s]
dt  = 0.005    # Time step for the calculation [s]
wst = 9.8    # Winch Stall torque - Torque at zero speed [Nm]
wzs = 3800/60*2*np.pi #Winch zero torque speed  - Speed where the winch has no torque[rad/s]
D=0.055 # Diameter of the cylinger collecting the wire [m]
l0 = 200 # Wirelength, witouth tention [m]
lf0 =600 # Preforce applied to the wire [N]

# Some global variables
cl = 0.0     # Lift Coefficient [-]
cd = 0.0     # Drag coefficient [-]
l    = [l0+lf0/k0*l0]       # Length of the line between the winch and the plane [m]
lw    = [0]        # Meters of line on the winch [m]
lf  = [lf0]    # Lineforce [N]
k    = k0        # Actual spring force of the line [N/m]
gamma = gamma0   # Angle between plane and ground [deg]
psi = 0.0      # Angle between line and ground [deg]
x   = [-l[0]]      # Positon of the plane in x direction [m]
y   = [0.0]      # Position of the plane in y direction [m]
u   = [v0*math.cos(rad(gamma))]      # Plane velocity in x direction [m/s]
v   = [v0*math.sin(rad(gamma))]      # Plane velocity in y direction [m/s]
T  = [0.0]     # Accumulated time [s]
attAng = [0] # Angle of attack [deg]
velAng = [gamma] # The planes velocity angle [deg]
omega = [0]    # Speed of the winch [rad/s]


def calcCd():
    """
    Returns the drag coefficient
    """
    return cl**2/13 + 0.002

def calcCl():
    """
    Returns the lift coefficient
    """
    return 2*np.pi*rad(attAng[-1])

def calcVelAng():
    """
    Returns the angle of the velocity vector of the plane
    """
    return deg(math.atan2(v[-1],u[-1]))


def calcAttAng():
    """
    Returns the angle of attack for the plane
    """
    return gamma-velAng[-1]

def calcGamma():
    """
    Returns the plane angle.
    Assumes the plane flies with gamma0 degrees towards the line all the time
    """


    return gamma0-psi#deg(math.atan2(v[-1],u[-1]))



def calcPsi():
    """
    Returns the line angle

    """

    return deg(math.atan2(y[-1],(-x[-1])))

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


def Swinch():
    """
    Returns the length of wire which the winch collects during 1 dt
    If the torque is bigger than the stall torque, It is assumed that the winch stops and does not reverse
    """
    M=lf[-1]*(D/2) # Torque acting on the cylinder [Nm]
    omega.append(min(max(0,(1-M/wst)),1)*wzs) # Rotational speed of the winch [rad/s]
    S = omega[-1]*(D/2)*dt # The amount of line the winch collects [m]
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
    return ((x[-1])**2+y[-1]**2)**0.5

def fLine():
    """
    Returns the force in the line
    dF=(dL+winchspeed)/k
    """
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
    fx,fy = sumForces()
    ax=fx/pm
    ay=fy/pm

    u.append(u[-1]+ax*dt)
    v.append(v[-1]+ay*dt)

    x.append(x[-1]+u[-1]*dt)
    y.append(max(y[-1]+v[-1]*dt,0))


def simulate():
    global gamma,psi,wf,velAng,attAng,cd,cl
    """
    Runs the simulation
    """
    while T[-1]<=Tmax and y[-1] >= -10.0:

        psi = calcPsi()
        gamma = calcGamma()
        velAng.append(calcVelAng())
        attAng.append(calcAttAng())
        cl=calcCl()
        cd=calcCd()
        l.append(lLine())
        lf.append(fLine())
        Euler()
        T.append(T[-1]+dt)
        #print T[-1],attAng,gamma
        if x[-1]>0:
            break


def plotSim():
    """
    Plots some graphs providing some information
    """
    subplot(3,1,1)
    xlabel("X-Position [m]")
    ylabel("Y-Position [m]")
    plot(x,y)
    #axes().set_aspect('equal', 'datalim')
    subplot(3,1,2)
    plot(T,lf)
    xlabel("Time [s]")
    ylabel("Force in wire [N]")
    subplot(3,1,3)
    plot(T,attAng)
    xlabel("Time [s]")
    ylabel("Angle of attack [deg]")
    show()


if __name__=="__main__":
    simulate()
    print "Hmax: ",max(y),"Energy: ",y[-1]*g*pm+0.5*pm*(u[-1]**2+v[-1]**2)**0.5
    plotSim()
