from scipy.integrate.quadpack import dblquad
from scipy.spatial.distance import euclidean
from math import sqrt
import numpy as np




def main():
    res = 0    
    
    X0, X1, X2 = calcLagrange((0, 0), (1, 0), (0, 1))
    lagrange = [X0, X1, X2]
    for i in range(3): # Laufvariable für i = 0 bis 5
        ans, err = dblquad(integrand, a=0, b=1, gfun=0, hfun=1/2, args=[i, lagrange]) # Integralfunktion, welche fkt, äußere schranken, innere schranken, Laufvariable (für benutzte fkt)
        res += abs(ans) # aufzsummieren
        print("i=" + str(i) + ": " + '{:.4f}'.format(abs(ans))) # jedes Integral in Laufvariable ausgegeben
        #print("Err: " + str(err))

    print("Ergebnis: " + '{:.4f}'.format(res)) # aufsummiertes Endergebnis


def integrand(y, x, i, lagrange):   # Inneres vom Integral
    fxy = 0
    li = 0
    if i == 0: 
        li = lagrange[0][0] + lagrange[0][1] * x + lagrange[0][2] * y
        fxy = euclideanNorm(0,0)
    elif i == 1: 
        li = lagrange[1][0] + lagrange[1][1] * x + lagrange[1][2] * y
        fxy = euclideanNorm(1,0)
    elif i == 2: 
        li = lagrange[2][0] + lagrange[2][1] * x + lagrange[2][2] * y
        fxy = euclideanNorm(0,1)
    return fxy * li


def calcLagrange(xy0, xy1, xy2):
    matrix = np.array([[1, xy0[0], xy0[1]], [1, xy1[0], xy1[1]], [1, xy2[0], xy2[1]]])
    
    l0 = np.array([1,0,0])
    l1 = np.array([0,1,0])
    l2 = np.array([0,0,1])
    
    X0 = np.linalg.inv(matrix).dot(l0)
    X1 = np.linalg.inv(matrix).dot(l1)
    X2 = np.linalg.inv(matrix).dot(l2)
    return X0, X1, X2


def euclideanNorm(x, y):
    return euclidean(x, y)
    return sqrt(x*x + y*y)


main()