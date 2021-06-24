"""
matrix.py
"""
import copy
import numpy as np
"""
printMat

This function will print a matrix in a readable format. You will not need to
alter this function.

INPUTS
mat: the matrix represented as a list of lists in row major form.

OUTPUTS
s: a string representing the matrix formatted nicely.
"""
def printMat(mat):
    s = ''
    for row in mat:
        s = s + ''
        for col in row:
            # Display 2 digits after the decimal, using 9 chars.
            s = s + ('%9.2e' % col) + ' '
        s = s + '\n'
    return s

"""
Matrix Class

This class will have code that implements:
- matrix multiplication
- LU factorization
- backward substitution
- forward substitution
- Gaussian elimination
"""
class Matrix:

    """
    Class attributes:
    mat:     the matrix itself, represented as a list of lists.
    numRows: the number of rows in the matrix.
    numCols: the number of columns in the matrix.
    L:       the lower triangular matrix from the LU Factorization.
    U:       the upper triangular matrix from the LU Factorization.
    P:       the permutation matrix from the LU Factorization.
    """

    # Constructor method.
    def __init__(self, mat):
        self.mat = mat
        self.numRows = len(mat)
        self.numCols = len(mat[0])
        self.L = None
        self.U = None
        self.P = None

    # Special method used for printing this Matrix.
    # You will not have to alter this function.
    def __repr__(self):
        s = ''
        s += 'The %dx%d Matrix itself:\n\n' % (self.numRows, self.numCols)
        s += printMat(self.mat)
        s += '\n'
        if self.L != None:
            s += 'The lower triangular matrix L:\n\n'
            s += printMat(self.L.mat)
            s += '\n'
        if self.U != None:
            s += 'The upper triangular matrix U:\n\n'
            s += printMat(self.U.mat)
            s += '\n'
        if self.P != None:
            s += 'The permutation matrix P:\n\n'
            s += printMat(self.P.mat)
        return s
    '''
    Matrix Multiplication Function

    @input B    a Matrix object B which is the right operand of multiplication
    
    @output C   resulting Matrix multiplication object
    
    '''
    
    def matMult(self,B): # inputs: matrix to be multiplied with 
        # raise exception if the inner dimensions do not agree
        if self.numCols != B.numRows:
            raise Exception ("Inner dimensions do not agree!")
        # preallocate matrix C with zeros in place of all entries
        C = [[0 for x in range (B.numCols)] for y in range (self.numRows)]
        # iterate through the rows of self
        for i in range (self.numRows):
            # iterate through the columns of B
            for j in range (B.numCols):
                # iterate through the rows of B
                for k in range (B.numRows):
                    # compute the dot product of the row entries of self 
                    # and the column entries of B 
                    C[i][j] += self.mat[i][k] * B.mat[k][j] 
        return Matrix(C) # output: resulting matrix

    '''
    LU Factorization function that calculates the Lower triangular Matrix object,
    the upper triangular Matrix Object, and the permutation Matrix object
    
    No inputs or outputs
    '''
    def LUfact(self):
        # check singularity
        for j in range (0, len(self.mat) - 1):
            checkPivotRow = j
            for i in range (j, len(self.mat)):
                if abs(self.mat[i][j]) > abs(self.mat[checkPivotRow][j]):
                    checkPivotRow = i
                    
        if checkPivotRow == j:
            raise Exception("Input matrix is singular!")
            
        # preallocate matrix L with zeros in place of all entries
        self.L = [[0 for x in range (self.numCols)] for y in range (self.numRows)]
        for r in range (0, self.numRows):
            for c in range (0, self.numCols):
                if r == c:
                    self.L[r][c] += 1
        
        # preallocate matrix P with zeros in place of all entries
        self.P = [[0 for x in range (self.numCols)] for y in range (self.numRows)]
        for r in range (0, self.numRows):
            for c in range (0, self.numCols):
                if r == c:
                    self.P[r][c] += 1
                    
        # copy input matrix to U in order to prevent changes to self
        self.U = copy.deepcopy(self.mat)
        
        # forward elimination
        for j in range (0, len(self.U) - 1):
            pivotRow = j
            for i in range (j, len(self.U)):
                if abs(self.U[i][j]) > abs(self.U[pivotRow][j]):
                    pivotRow = i 
                    
            # row swaps for U
            tempRow = self.U[j]
            self.U[j] = self.U[pivotRow]
            self.U[pivotRow] = tempRow
            # row swaps for P
            tempRow = self.P[j]
            self.P[j] = self.P[pivotRow]
            self.P[pivotRow] = tempRow
            # row swaps for L 
            if j < i:
                # swap rows that are only below the diagonal if j < i
                for c in range (0,j):
                    self.L[j][c] = self.L[i][c]
            else:
                # otherwise swap normally
                tempRow = self.L[j]
                self.L[j] = self.L[pivotRow]
                self.L[pivotRow] = tempRow
            
            # update entries in U
            for k in range (j+1, len(self.U)):
                l = self.U[k][j] / self.U[j][j]
                # set L[k][j] to the value of lkj
                self.L[k][j] = l
                for c in range (j, len(self.U)):
                    self.U[k][c] += -l * self.U[j][c]
                    
        # Make the L, U, and P lists into Matrix Objects 
        self.L = Matrix(self.L)
        self.U = Matrix(self.U)
        self.P = Matrix(self.P)
        return

    '''
    Performs back substitution on the Matrix object C
    
    @param c    Matrix object that is result of forward substitution on b
    
    @output x   final solution to Gaussian elimination
    
    '''
    def backSub(self,c):
        # Initialize empty x column vector
        x = [[0] for i in range (self.numRows)]
            
        # Use upper triangular matrix to find Xn
        x[self.numRows-1][0] += c.mat[self.numRows-1][0]/self.U.mat[self.numRows-1][self.numCols-1]
        
        for i in range (self.numRows-1):
            j = self.numRows - 2 - i
            knownVals = 0
            for k in range (j+1, self.numRows):
                knownVals += self.U.mat[j][k] * x[k][0]
            x[j][0] += (1/self.U.mat[j][j])*(c.mat[j][0] - knownVals)
        return Matrix(x)
    '''
    Performs forward substitution on the image vector, which is the Matrix object
    that is b
    
    @param b    Matrix object that is our initial image vector 
    
    @output c   Resulting Matrix object after forward substitution on b
    
    '''
    def forSub(self,b):
        c = [[0] for i in range (self.numRows)]
        c[0][0] = b.mat[0][0]
        knownVals = 0
        for j in range (1, self.numRows):
            knownVals = 0
            for k in range (1, j):
                knownVals += self.L.mat[j][k] * c[k][0]
            c[j][0] += b.mat[j][0] - knownVals
        return Matrix(c)

    '''
    Wrapper functions that calls helper functions to perform Gaussian elimination
    on any Matrix object with a given image vector, b
    
    @param b    Matrix object that is our given image vector
    
    @output     Matrix object that is the solution to our Gaussian elimination
    '''
    def gaussElim(self,b):
        # perform LUfact if it has not been done
        if self.L == None or self.U == None or self.P == None:
            self.LUfact() 
        bHat = self.P.matMult(b)
        c = self.forSub(bHat)
        x = self.backSub(c)
        return x