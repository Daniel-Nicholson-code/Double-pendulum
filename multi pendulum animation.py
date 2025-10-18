import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from scipy.integrate import odeint








#   PRE-PROCESSING

#Defining constants:
#g := gravitational constant
g = 9.81
pi = 3.14
#t := array containg the elapsed time for each frame
t = np.arange(0, 20, 0.03)
#L1,L2 := lenth of each rod
L1, L2 = 1, 1
#m1,m2 := mass of each joint
m1, m2 = 2, 2
#pendulumCount is the amount of pendulums in the simulation
pendulumCount = 25

#Creating plot
fig, ax1 = plt.subplots(1, 1, figsize=(5, 5))








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

#Creating simulation variables
theta1 = np.zeros(pendulumCount).tolist()
theta2 = np.zeros(pendulumCount).tolist()
x1 = np.zeros(pendulumCount).tolist()
x2 = np.zeros(pendulumCount).tolist()
y1 = np.zeros(pendulumCount).tolist()
y2 = np.zeros(pendulumCount).tolist()

#Simulating each pendulum
for i in range(pendulumCount):

    # Initial conditions: theta1, dtheta1/dt, theta2, dtheta2/dt
    z0 = np.array([np.pi/2 + i*0.002, 0, np.pi/2 + i*0.002, 0])

    # Do the numerical integration of the equations of motion
    z = odeint(deriv, z0, t)

    # Unpack z and theta as a function of time
    theta1[i], theta2[i] = z[:,0], z[:,2]

    #Turning the angles into positions so we can draw them
    x1[i] = L1 * np.sin(theta1[i])
    y1[i] = 0 - L1 * np.cos(theta1[i])
    x2[i] = x1[i] + L2 * np.sin(theta2[i])
    y2[i] = y1[i] - L2 * np.cos(theta2[i])








#   POST-PROCCESSING

#PLotting the data:
#Creating graphs showing the positions of each pendulum
rod1 = np.zeros(pendulumCount).tolist()
rod2 = np.zeros(pendulumCount).tolist()

for i in range(pendulumCount):

    rod1[i] = ax1.plot([0, x1[i][0]], [0, y1[i][0]], c=(i/pendulumCount, 0, 0))[0]
    rod2[i] = ax1.plot([x1[i][0], x2[i][0]], [y1[i][0], y2[i][0]], c=(i/pendulumCount, 0, 0))[0]

#Defining behaviour of the animation
def update(frame):

    #Updating what is being displayed for each element
    for i in range(pendulumCount):

        rod1[i].set_xdata([0, x1[i][frame]])
        rod1[i].set_ydata([0, y1[i][frame]])

        rod2[i].set_xdata([x1[i][frame], x2[i][frame]])
        rod2[i].set_ydata([y1[i][frame], y2[i][frame]])

    return (rod1, rod2)

#Presenting the plot nicely
ax1.set(xlim=[-2.5, 2.5], ylim=[-2.5, 2.5], xlabel="x", ylabel="y")
ani = animation.FuncAnimation(fig=fig, func=update, frames=np.size(t), interval=30)
plt.subplots_adjust(right=0.9, top=0.9, left=0.1, bottom=0.1)

plt.show()