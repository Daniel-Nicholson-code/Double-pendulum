import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import time








#   PRE-PROCESSING

#Defining constants: 
g = 9.81
pi = 3.14
t = np.arange(0, 5, 0.03) #array containg the elapsed time for each frame
L1, L2 = 1, 1 #lenth of each rod
m1, m2 = 2, 2 #mass of each joint
N = 50 #size of the side lengths of the square grid of pixels

# Important: Complexity is O(N^2)

#keeping track of how long the program took to run
startTime = time.perf_counter()








# SIMULATION

#Return the first derivative of y = theta1, omega1, theta2, omega2
def deriv(z, t):

    #Unpacking z
    theta1, omega1, theta2, omega2 = z

    #Precomputing these to save computation time
    c, s = np.cos(theta1-theta2), np.sin(theta1-theta2)

    #Performing the computation
    d_theta1 = omega1
    d_omega1 = (m2*g*np.sin(theta2)*c - m2*s*(L1*omega1**2*c + L2*omega2**2) -
             (m1+m2)*g*np.sin(theta1)) / L1 / (m1 + m2*s**2)
    d_theta2 = omega2
    d_omega2 = ((m1+m2)*(L1*omega1**2*s - g*np.sin(theta2) + g*np.sin(theta1)*c) + 
             m2*L2*omega2**2*s*c) / L2 / (m1 + m2*s**2)
    
    return d_theta1, d_omega1, d_theta2, d_omega2

#Creating simulation variables
theta1 = np.zeros([N, N, len(t)])
theta2 = np.zeros([N, N, len(t)])

#Simulating each pendulum
progress = 0
print("starting...", end="\r")

for i in range(N):

    #Indicating how much progress has been made in computing the simulation
    temp = round(i*100/N)
    if progress != temp:
        progress = temp
        print(f"{progress}% complete...", end="\r")

    for j in range(N):

        #Initial conditions: theta1, dtheta1/dt, theta2, dtheta2/dt
        z0 = np.array([(i/N)*2*pi - pi, 0, (j/N)*2*pi - pi, 0])

        #Do the numerical integration of the equations of motion
        z = odeint(deriv, z0, t)

        #Unpack z and theta as a function of time
        theta1[i][j], theta2[i][j] = z[:,0], z[:,2]

#Displaying how long the program took to run
print(f"program finished in {(time.perf_counter()-startTime):.3f}s")






#   POST-PROCCESSING

#Setting up the plot
fig, ax = plt.subplots()
image = ax.pcolormesh(theta2[:, :, 0], cmap=plt.cm.CMRmap, vmin=0, vmax=4*pi)
ax.set(xlabel="theta2", ylabel="theta1")

#Displaying the animation
for frame in range(len(t)):

    # Updating the image with the next frame
    image.set_array(theta2[:, :, frame])
    
    plt.pause(0.05)

plt.show()