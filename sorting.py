"""
sorting.py
"""

"""
merge: the helper function for MergeSort that merges two sorted lists.

INPUTS:
L: original input list from MergeSort(L)
L1: sorted list 1 
L2: sorted list 2 

OUTPUTS:
L: merged sorted lists
"""

def merge(L, L1, L2):
    # initialize counter variables a, b, c
    a = 0 # iterates through left half
    b = 0 # iterates through right half
    c = 0 # iterates through the main list
    # make appropriate swaps
    while a < len(L1) and b < len(L2):
        if L1[a] < L2[b]:
            L[c] = L1[a]
            a += 1
        else:
            L[c] = L2[b]
            b += 1
        c += 1
    # assign the rest (one of the halves might be incompletely traversed)
    while a < len(L1):
        L[c] = L1[a]
        a += 1
        c += 1
    while b < len(L2):
        L[c] = L2[b]
        b += 1
        c += 1
    return L
"""
MergeSort: A divide-and-conquer recursive sorting algorithm

INPUTS: 
L: an unsorted list 

OUTPUTS:
L: a sorted list
"""
def MergeSort(L):
    # if the input list has more than 1 element
    if len(L) > 1: 
        # evaluate mid entry and assign left and right half lists
        mid = len(L) // 2
        left = L[:mid]
        right = L[mid:]
        # recursively call each half
        MergeSort(left)
        MergeSort(right)
        # call merge 
        merge(L, left, right)
    return L

"""
partition: the helper function for QuickSort that partitions the input
           list L based on the value of the pivot at L[p].
           
INPUTS:
L: unsorted list 
p: starting (pivot) index
q: ending index

OUTPUTS:
newL: new list
i: new pivot index     
"""

def partition(L, p, q):
    # define initial pivot
    pivot = L[p] 
    lower = p + 1 # index of lower bound 
    upper = q # index of upper bound 
    while True:
        # traverse if order is correct
        while lower <= upper and L[upper] >= pivot:
            upper -= 1 # decrement
        while lower <= upper and L[lower] <= pivot:
            lower += 1 # increment
        # swap if order is incorrect
        if lower <= upper:
            L[lower], L[upper] = L[upper], L[lower]
        else:
            break
    # define new pivot
    L[p], L[upper] = L[upper], L[p]
    newL = L 
    i = upper
    return (newL, i)

"""
QuickSort: A divide-and-conquer recursive sorting algorithm

INPUTS:
L: unsorted list 
p: starting (pivot) index
q: ending index

OUTPUTS:
L: sorted list
"""
def QuickSort(L, p, q):
    # base case
    if p >= q:
        return
    # call partition function
    part = partition(L, p, q)
    # call quick sort on left half
    QuickSort(L, p, part[1]-1)
    # call quick sort on right half
    QuickSort(L, part[1]+1, q)
    # return sorted list 
    return L