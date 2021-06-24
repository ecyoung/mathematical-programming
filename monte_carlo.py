"""
monte_carlo.py
"""

import math
import random
import numpy as np

def buffonNeedle(trials): # input a number of trials in the range of 100,000 to 10,000,000
    counter  = 0 # initialize while loop counter to 0
    crossings = 0 # initialize the number of crossings to 0
    while (counter < trials): 
        D = random.random() # randomly select a real number between 0 and 1
        # unit circle simulation
        X = (random.random() - 0.5) * 2 # randomly select a real number between -1 and 1 
        Y = random.random() # # randomly select a real number between 0 and 1
        thresh = math.sqrt(X**2+Y**2)
        if (thresh > 1): continue # go back to start of loop if out of unit circle
        else: d = Y / thresh # this is sine theta
        if d > D: crossings += 1 # crossing occurs when d > D
        counter += 1 # increment counter by 1
    if crossings == 0: pi = math.inf # to avoid the divide by 0 edge case
    pi = (2 * trials) / crossings 
    return pi # returns the approximate value of pi 
    
def func(x): # define function with input as x
    return x**2 # returns the function output

def erfunc(t): # define error function with input as t
    return (2*math.exp(-t**2)) / math.sqrt(math.pi) # returns the function output

def MCpoly(a, b, N): # define integral function with inputs a (lower bound), b(upper bound), and N(number of trials)
    total = 0 # initalize summation to 0
    lst = np.random.uniform(a, b, N) # this is a list of N random numbers between a and b 
    for i in lst: # add every item in the list applied with the function to the summation variable
        total += func(i) 
    apprx = ((b-a) / N) * total # calculate approximate integral
    return apprx # returns the integral approximation

def MCerf(x, N):
    total = 0 # initalize summation to 0
    lst = np.random.uniform(0, x, N) # this is a list of N random numbers between 0 and x 
    for i in lst: # add every item in the list applied with the error function to the summation variable
        total += erfunc(i)
    apprx = ((x-0) / N) * total # calculate approximate error integral
    return apprx # returns the error integral approximation