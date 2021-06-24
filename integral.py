"""
integral.py
"""

import math
import numpy as np
import matplotlib.pyplot as plt

"""
This is the function f(x) to be integrated.    
"""
def f(x):
    return math.cos(x/10)

"""
Right Point Rule

INPUTS:
N: number of subintervals 
a: integral lower bound
b: integral upper bound

OUTPUTS:
integral evaluation via the right point rule
"""
def right(N, a, b):
    # N is the number of subintervals, so there are N+1 collocation points
    n = N+1
    # subinterval width 
    # (total length of interval divided by the number of subintervals)
    deltaX = (b-a)/(N)
    # initialize summation counter
    summation = 0
    
    # we want to apply a summation from 1 to n-1 (inclusive)
    for i in range (1, n):
        summation += f(a + i*deltaX)
    return deltaX * summation

"""
Left Point Rule

INPUTS:
N: number of subintervals (there are N+1 collocation points)
a: integral lower bound
b: integral upper bound

OUTPUTS:
integral evaluation via the left point rule
"""
def left(N, a, b):
    # N is the number of subintervals, so there are N+1 collocation points
    n = N+1
    # subinterval width 
    # (total length of interval divided by the number of subintervals)
    deltaX = (b-a)/(N)
    # initialize summation counter
    summation = 0
    
    # we want to apply a summation from 0 to n-2 (inclusive)
    for i in range (0, n-1):
        summation += f(a + i*deltaX)
    return deltaX * summation

"""
Trapezoidal Rule

INPUTS:
N: number of subintervals (there are N+1 collocation points)
a: integral lower bound
b: integral upper bound

OUTPUTS:
integral evaluation via the trapezoidal rule
"""
def trap(N, a, b):
    # N is the number of subintervals, so there are N+1 collocation points
    n = N + 1
    # subinterval width 
    # (total length of interval divided by the number of subintervals)
    deltaX = (b-a)/(N)
    # initialize summation counter
    summation = 0
    
    # we appply the method to minimize operations from 1 to n-2 (inclusive)
    for i in range (1, n-1): 
        summation += f(a + i*deltaX)
    return (deltaX/2) * (f(a) + f(b) + 2*summation)

"""
main function
"""
if __name__ == '__main__':

    # Set the values of a and b.
    a = 0
    b = 1/2

    # Set the correct solution to compare against.
    correct = 10*math.sin(1/20)

    # Create the range of N values.
    N = [x for p in range(2,8) for x in (10**p,5*(10**p))]
    logN = [np.log10(n) for n in N]
    
    # initialize lists right point, left point, and trapezoidal rules
    valuesRight = []
    errorsRight = []
    
    valuesLeft = []
    errorsLeft = []
    
    valuesTrap = []
    errorsTrap = []
    
    # right point rule
    for n in N:
        rightApprox = right(n, a, b)
        valuesRight.append(rightApprox)
        errorsRight.append(abs(correct - rightApprox))
        
    logErrorsRight = [np.log10(err) for err in errorsRight]
    
    # left point rule
    for n in N:
        leftApprox = left(n, a, b)
        valuesLeft.append(leftApprox)
        errorsLeft.append(abs(correct - leftApprox))
    
    logErrorsLeft = [np.log10(err) for err in errorsLeft]
    
    # trapezoidal rule
    for n in N:
        trapApprox = trap(n, a, b)
        valuesTrap.append(trapApprox)
        errorsTrap.append(abs(correct - trapApprox))
        
    logErrorsTrap = [np.log10(err) for err in errorsTrap]
    print(errorsTrap)
    
    # right error plot
    plt.figure()
    fig, ax = plt.subplots()
    ax.plot(logN, logErrorsRight, label="Absolute Error Per Subinterval Count (log-log)")
    legend = ax.legend(loc="upper right")
    plt.title("Absolute Error of Right Point Rule with Increasing Subinterval Counts (log-log)")
    plt.xlabel("log(# of Subintervals)")
    plt.ylabel("log(Error)")
    plt.savefig("right.png", bbox_inches="tight") # Save as a png file.
    plt.show()
    
    # left error plot
    plt.figure()
    fig, ax = plt.subplots()
    ax.plot(logN, logErrorsLeft, label="Absolute Error Per Subinterval Count (log-log)")
    legend = ax.legend(loc="upper right")
    plt.title("Absolute Error of Left Point Rule with Increasing Subinterval Counts (log-log)")
    plt.xlabel("log(# of Subintervals)")
    plt.ylabel("log(Error)")
    plt.savefig("left.png", bbox_inches="tight") # Save as a png file.
    plt.show()
    
    # trap error plot
    plt.figure()
    fig, ax = plt.subplots()
    ax.plot(logN, logErrorsTrap, label="Absolute Error Per Subinterval Count (log-log)")
    legend = ax.legend(loc="upper right")
    plt.title("Absolute Error of Trapezoid Rule with Increasing Subinterval Counts (log-log)")
    plt.xlabel("log(# of Subintervals)")
    plt.ylabel("log(Error)")
    plt.savefig("trap.png", bbox_inches="tight") # Save as a png file.
    plt.show()