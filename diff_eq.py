"""
diff_eq.py
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from statistics import mean

"""
FE: Forward Euler (first order explicit)

INPUTS: 
w0: the natural frequency of the system
z: the damping ratio of the system
m: the mass of the system
w: the forcing frequency
x0: the initial condition
T: the final time to solve
N: the number of timesteps to use while solving

OUTPUTS:
x: a list of the displacement values at each timestep
t: a list of the times corresponding to the values of x 
"""
def FE(w0, z, m, w, x0, T, N):
    # Initialize A, dt
    # Preallocate x, t
    A = np.matrix([[0, 1], [-w0**2, -2*z*w0]])
    dt = T/N
    x = [0]*(N+1)
    t = [0]*(N+1)
    
    # Fill in time values
    # Set initial displacement as x[0]
    for i in range(1, N+1):
        t[i] = t[i-1] + dt
    x[0] = x0[0]
    
    # Loop through each timestep
    xN_Vec = x0
    for n in range(0, N):
        # Calculate b Vector using w, t, m
        b_Vec = np.matrix([[0], [math.cos(w*t[n])/m]])
        # Perform update step
        xN_Vec = xN_Vec + dt * np.matmul(A, xN_Vec) + dt * b_Vec
        # Update displacement
        x[n+1] = xN_Vec[0]
        
    return (x,t)

"""
BE: Backward Euler (first order implicit)

INPUTS:
w0: the natural frequency of the system
z: the damping ratio of the system
m: the mass of the system
w: the forcing frequency
x0: the initial condition
T: the final time to solve
N: the number of timesteps to use while solving

OUTPUTS: 
x: a list of the displacement values at each timestep
t: a list of the times corresponding to the values of x 
"""
def BE(w0, z, m, w, x0, T, N):
    # Initialize A, dt
    # Preallocate x, t
    A = np.matrix([[0, 1], [-w0**2, -2*z*w0]])
    dt = T/N
    x = [0]*(N+1)
    t = [0]*(N+1)
    
    # Fill in time values
    # Set initial displacement as x[0]
    for i in range(1, N+1):
        t[i] = t[i-1] + dt
    x[0] = x0[0]
    
    # Obtain "BE multiplicative" factor
    inversionMultiplier = np.linalg.inv(np.identity(2) - dt * A)
    
    # Loop through each timestep
    xN_Vec = x0
    for n in range(0, N):
        # Calculate b Vector using w, t+1 (use t+1 for BE), m
        b_Vec = np.matrix([[0], [math.cos(w*t[n+1])/m]])
        # Perform update step
        xN_Vec = np.matmul(inversionMultiplier, xN_Vec + dt * b_Vec)
        # Update displacement
        x[n+1] = xN_Vec[0]   
        
    return (x,t)

"""
CN: Crank-Nicolson (second order implicit)

INPUTS: 
w0: the natural frequency of the system
z: the damping ratio of the system
m: the mass of the system
w: the forcing frequency
x0: the initial condition
T: the final time to solve
N: the number of timesteps to use while solving

OUTPUTS:
x: a list of the displacement values at each timestep
t: a list of the times corresponding to the values of x 
"""
def CN(w0, z, m, w, x0, T, N):
    # Initialize A, dt
    # Preallocate x, t
    A = np.matrix([[0, 1], [-w0**2, -2*z*w0]])
    dt = T/N
    x = [0]*(N+1)
    t = [0]*(N+1)
    
    # Fill in time values
    # Set initial displacement as x[0]
    for i in range(1, N+1):
        t[i] = t[i-1] + dt
    x[0] = x0[0]
    
    # Obtain "multiplicative" factor
    multiplier = np.identity(2) + dt/2 * A
    inversionMultiplier = np.linalg.inv(np.identity(2) - dt/2 * A)
    
    # Loop through each timestep
    xN_Vec = x0
    for n in range(0, N):
        # Calculate b vector for t and t+1 using w and 
        b_Vec_0 = np.matrix([[0], [math.cos(w*t[n])/m]])
        b_Vec_1 = np.matrix([[0], [math.cos(w*t[n+1])/m]])
        # Perform update step
        xN_Vec = np.matmul(inversionMultiplier, np.matmul(multiplier, xN_Vec) + dt * ((b_Vec_1 + b_Vec_0)/2))
        # Update displacement
        x[n+1] = xN_Vec[0]
        
    return (x,t)

"""
RK4: Fourth-Order Runge-Kutta

INPUTS:
w0: the natural frequency of the system
z: the damping ratio of the system
m: the mass of the system
w: the forcing frequency
x0: the initial condition
T: the final time to solve
N: the number of timesteps to use while solving

OUTPUTS: 
x: a list of the displacement values at each timestep
t: a list of the times corresponding to the values of x 
"""
def RK4(w0, z, m, w, x0, T, N):
    # Initialize A, dt
    # Preallocate x, t
    A = np.matrix([[0, 1], [-w0**2, -2*z*w0]])
    dt = T/N
    x = [0]*(N+1)
    t = [0]*(N+1)
    
    # Fill in time values
    # Set initial displacement as x[0]
    for i in range(1, N+1):
        t[i] = t[i-1] + dt
    x[0] = x0[0]
    
    # Loop through each timestep
    xN_Vec = x0
    for n in range(0, N):
        # Calculate b vector for t, t + dt/2, and t + dt
        b_Vec_k1 = np.matrix([[0], [math.cos(w*t[n])/m]])
        b_Vec_k2k3 = np.matrix([[0], [math.cos(w*(t[n]+dt/2))/m]])
        b_Vec_k4 = np.matrix([[0], [math.cos(w*(t[n]+dt))/m]])
        # Calculate k1, k2, k3, k4 for each iteration
        k1 = dt * np.matmul(A, xN_Vec) + dt * b_Vec_k1
        k2 = dt * np.matmul(A, xN_Vec + k1/2) + dt * b_Vec_k2k3
        k3 = dt * np.matmul(A, xN_Vec + k2/2) + dt * b_Vec_k2k3
        k4 = dt * np.matmul(A, xN_Vec + k3) + dt * b_Vec_k4
        # Perform update step
        xN_Vec = xN_Vec + k1/6 + k2/3 + k3/3 + k4/6
        # Update displacement
        x[n+1] = xN_Vec[0]
    
    return (x,t)

"""
Function of Exact Solution for Error

INPUTS:
t: time 

OUTPUTS:
float: exact solution
"""
def exact(t):
    return (1/2) * (math.sin(t) - t * math.exp(-t))

'''
Function for Best-Fit Slope

INPUTS:
xs: x values
ys: y values 

OUTPUTS:
float: slope
'''

def best_fit_slope(xs,ys):
    m = (((mean(xs)*mean(ys)) - mean(xs*ys)) /
         ((mean(xs)**2) - mean(xs**2)))
    return m

"""
main function
"""
if __name__ == '__main__':
    """
    Part 3: Testing the methods
    """
    # Initialize array N, that contains the number of timesteps.
    # Values of N should be: 100, 1000, 10^4, 10^5, 10^6, 5 different N vals
    # Set total time T=10
    N = [10**(i+2) for i in range(0, 5)]
    T = 10
    
    ### For each method, calculate error at time T = 10 for a range of N values ###
    # Preallocate arrays for error
    err_FE = [0]*5
    err_BE = [0]*5
    err_CN = [0]*5
    err_RK4 = [0]*5
    
    # Set variables used for testing
    w0 = 1
    z = 1
    m = 1
    w = 1
    x0 = np.matrix([[0], [0]])
        
    # Set errors
    for i in range(0, 5):
        err_FE[i] = abs(exact(T) - FE(w0, z, m, w, x0, T, N[i])[0][-1])
        err_BE[i] = abs(exact(T) - BE(w0, z, m, w, x0, T, N[i])[0][-1])
        err_CN[i] = abs(exact(T) - CN(w0, z, m , w, x0, T, N[i])[0][-1])
        err_RK4[i] = abs(exact(T) - RK4(w0, z, m, w, x0, T, N[i])[0][-1])
    
    # Plot FE
    plt.figure()
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.plot(N, np.array(err_BE).flatten(), label="Absolute Error Per Timestep Count (log-log)")
    legend = ax.legend(loc="upper right")
    plt.title("Absolute Error of Backwards Euler with Increasing Timesteps (log-log)")
    plt.xlabel("log(# of Timesteps)")
    plt.ylabel("log(Error)")
    plt.savefig("BE_error.png", bbox_inches="tight") # Save as a png file.
    print("FE Slope: ", best_fit_slope(np.log10(N), np.log10(np.array(err_FE).flatten())))

    # Plot BE
    plt.figure()
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.plot(N, np.array(err_FE).flatten(), label="Absolute Error Per Timestep Count (log-log)")
    legend = ax.legend(loc="upper right")
    plt.title("Absolute Error of Forward Euler with Increasing Timesteps (log-log)")
    plt.xlabel("log(# of Timesteps)")
    plt.ylabel("log(Error)")
    plt.savefig("FE_error.png", bbox_inches="tight") # Save as a png file.
    print("BE Slope: ", best_fit_slope(np.log10(N), np.log10(np.array(err_BE).flatten())))
    
    # Plot CN
    plt.figure()
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.plot(N, np.array(err_CN).flatten(), label="Absolute Error Per Timestep Count (log-log)")
    legend = ax.legend(loc="upper right")
    plt.title("Absolute Error of Crank Nicolson with Increasing Timesteps (log-log)")
    plt.xlabel("log(# of Timesteps)")
    plt.ylabel("log(Error)")
    plt.savefig("CN_error.png", bbox_inches="tight") # Save as a png file.
    print("CN Slope: ", best_fit_slope(np.log10(N), np.log10(np.array(err_CN).flatten())))
    
    # Plot RK4
    plt.figure()
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.plot(N, np.array(err_RK4).flatten(), label="Absolute Error Per Timestep Count (log-log)")
    legend = ax.legend(loc="upper right")
    plt.title("Absolute Error of Fourth Order Runge-Kutta with Increasing Timesteps (log-log)")
    plt.xlabel("log(# of Timesteps)")
    plt.ylabel("log(Error)")
    plt.savefig("RK4_error.png", bbox_inches="tight") # Save as a png file.
    
    '''
    Part 4: Beats and Resonance (FE)
    '''
    
    # Set variables used for beats and resonance
    w0 = 1
    z = 0
    m = 1
    x0 = np.matrix([[0], [0]])
    # set total time T = 100
    T = 100
    # set N = 10^5 to let dt = 10^-3 given T = 100
    N = 10**5
    # Case 1: w = 0.8
    # Plot x(t) versus time
    plt.figure()
    fig, ax = plt.subplots()
    ax.plot(np.array(FE(w0,z,m,0.8,x0,T,N)[1]).flatten(), np.array(FE(w0,z,m,0.8,x0,T,N)[0]).flatten())
    plt.title("x(t) v.s. t when w = 0.8")
    plt.xlabel("t")
    plt.ylabel("x(t)")
    plt.savefig("P4-1.png", bbox_inches="tight") # Save as a png file.
    
    # Case 2: w = 0.9
    # Plot x(t) versus time
    plt.figure()
    fig, ax = plt.subplots()
    ax.plot(np.array(FE(w0,z,m,0.9,x0,T,N)[1]).flatten(), np.array(FE(w0,z,m,0.9,x0,T,N)[0]).flatten())
    plt.title("x(t) v.s. t when w = 0.9")
    plt.xlabel("t")
    plt.ylabel("x(t)")
    plt.savefig("P4-2.png", bbox_inches="tight") # Save as a png file.
    
    # Case 3: w = 1.0
    # Plot x(t) versus time
    plt.figure()
    fig, ax = plt.subplots()
    ax.plot(np.array(FE(w0,z,m,1.0,x0,T,N)[1]).flatten(), np.array(FE(w0,z,m,1.0,x0,T,N)[0]).flatten())
    plt.title("x(t) v.s. t when w = 1.0")
    plt.xlabel("t")
    plt.ylabel("x(t)")
    plt.savefig("P4-3.png", bbox_inches="tight") # Save as a png file.
    
    '''
    Part 5: Frequency Response Function (FE)
    '''
    # set variables used for frequency response function
    w0 = 1
    z = 1/10
    m = 1
    x0 = np.matrix([[0], [0]])
    # set total time T = 100
    T = 100
    # set N = 10^5 to let dt = 10^-3 given T = 100
    N = 10**5
    # initialize frequency list
    freqList = []
    # initialize maximum displacement list 
    maxDispList = []
    # initialize w = 0.1
    w = 0.1
    while w <= 10:
        freqList.append(w)
        # ValueError: The truth value of an array with more than one element is ambiguous. Use a.any() or a.all()
        maxDisp = max(np.abs(np.array(FE(w0,z,m,w,x0,T,N)[0]).flatten()))
        maxDispList.append(maxDisp)
        w += 0.1
    # plot the maximum displacement versus frequency on a log-log scale
    plt.figure()
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.plot(np.array(freqList),np.array(maxDispList))
    plt.xlabel("Frequency (w)")
    plt.ylabel('Maximum Displacement (max(x))')
    plt.title("Maximum Displacement v.s. Frequency")
    plt.savefig("P5", bbox_inches="tight") # Save as a png file.