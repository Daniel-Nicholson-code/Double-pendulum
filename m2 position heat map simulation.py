import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from matplotlib.patches import Circle, Polygon
from scipy.integrate import odeint
from IPython import display

#Creating plot
fig, ax1 = plt.subplots(1, 1, figsize=(5, 5))








#   PRE-PROCESSING

#Defining constants:
#g := gravitational constant
g = 9.81
pi = 3.14
#t := array containg the elapsed time for each frame
t = np.arange(0, 1000, 0.03)
#L1,L2 := lenth of each rod
L1, L2 = 1, 1
#m1,m2 := mass of each joint
m1, m2 = 2, 2








# SIMULATION

#Return the first derivative of y = theta1, omega1, theta2, omega2
def deriv(z, t):

    #Unpacking z
    theta1, omega1, theta2, omega2 = z

    #Precomputing these to save computation time
    c, s = np.cos(theta1-theta2), np.sin(theta1-theta2)

    d_theta1 = omega1
    d_omega1 = (m2*g*np.sin(theta2)*c - m2*s*(L1*omega1**2*c + L2*omega2**2) -
             (m1+m2)*g*np.sin(theta1)) / L1 / (m1 + m2*s**2)
    d_theta2 = omega2
    d_omega2 = ((m1+m2)*(L1*omega1**2*s - g*np.sin(theta2) + g*np.sin(theta1)*c) + 
             m2*L2*omega2**2*s*c) / L2 / (m1 + m2*s**2)
    
    return d_theta1, d_omega1, d_theta2, d_omega2

# Initial conditions: theta1, dtheta1/dt, theta2, dtheta2/dt.
z0 = np.array([np.pi/2, 0, np.pi/2, 0])

# Do the numerical integration of the equations of motion
z = odeint(deriv, z0, t)

# Unpack z and theta as a function of time
theta1, theta2 = z[:,0], z[:,2]

#Turning the angles into positions so we can draw them
x1 = L1 * np.sin(theta1)
y1 = 0 - L1 * np.cos(theta1)
x2 = x1 + L2 * np.sin(theta2)
y2 = y1 - L2 * np.cos(theta2)








#   POST-PROCCESSING

#First, cleaning the data:
#Function that inputs an angle n and returns it in the form: -pi <= n <= pi
def cleanAngleData(n):
    
    #If n is outside the range, puts it back into range
    while n > pi: n -= 2*pi
    while n < -pi: n += 2*pi
    
    return n

#Using this function to put the angle data in the right format
theta1 = list(map(cleanAngleData, theta1))
theta2 = list(map(cleanAngleData, theta2))

#Detects if there is a jump between values and removes them from the plot
#This jumps are caused by the data cleaning we just did
for i in range(len(theta1) - 1):
    
    #If there is a sign change and a large jump in magnitude then there is a jump
    a, b = theta1[i], theta1[i+1]
    if (np.sign(a) != np.sign(b)) and (abs(a) > pi/2):
        #Removing the jump from the plot
        theta1[i] = np.nan

    #Doing the same for the other list of angles
    a, b = theta2[i], theta2[i+1]
    if (np.sign(a) != np.sign(b)) and (abs(a) > pi/2):
        theta2[i] = np.nan

#PLotting the data:
#Creating a graph showing the path of the second mass
step = 4
for i in range(0, len(t), step-1):
    ax1.plot(x2[i:i+step], y2[i:i+step], c="black", alpha=0.05)


#Presenting the plot nicely
fig.suptitle("Double pendulum m2 position heat map")
ax1.set(xlim=[-2.5, 2.5], ylim=[-2.5, 2.5], xlabel="x", ylabel="y")

plt.subplots_adjust(right=0.9, top=0.9, left=0.1, bottom=0.1)



plt.show()