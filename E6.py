import numpy as np
from scipy.integrate.quadpack import dblquad
from scipy.spatial.distance import euclidean
from math import sqrt


def main():
    X0, X1, X2, X3, X4, X5 = calcLagrange((0, 0), (1, 0), (0, 1), (0.5, 0), (0, 0.5), (0.5, 0.5)) # Lagrange Basen berechnen
    lagrange = [X0, X1, X2, X3, X4, X5]
    
    res = 0
    for i in range(6): # Laufvariable für i = 0 bis 5
        ans, err = dblquad(integrand, a=0, b=1, gfun=0, hfun=1/2, args=[i, lagrange]) # Integralfunktion, welche fkt, äußere schranken, innere schranken, Laufvariable (für benutzte fkt)
        res += abs(ans) # aufzsummieren
        print("i=" + str(i) + ": " + '{:.4f}'.format(abs(ans))) # jedes Integral in Laufvariable ausgegeben

    print("Ergebnis: " + '{:.4f}'.format(res)) # aufsummiertes Endergebnis


def integrand(y, x, i, lagrange): # inneres vom Integral (6 Stützstellen)
    fxy = 0
    li = 0
    if i == 0: 
        li = lagrange[0][0] + lagrange[0][1]*x + lagrange[0][2]*x*x - lagrange[0][3]*y + lagrange[0][4]*y*y + lagrange[0][5]*x*y
        fxy = euclideanNorm(0,0)
    elif i == 1: 
        li = lagrange[1][0] + lagrange[1][1]*x + lagrange[1][2]*x*x - lagrange[1][3]*y + lagrange[1][4]*y*y + lagrange[1][5]*x*y
        fxy = euclideanNorm(1,0)
    elif i == 2: 
        li = lagrange[2][0] + lagrange[2][1]*x + lagrange[2][2]*x*x - lagrange[2][3]*y + lagrange[2][4]*y*y + lagrange[2][5]*x*y
        fxy = euclideanNorm(0,1)
    elif i == 3: 
        li = lagrange[3][0] + lagrange[3][1]*x + lagrange[3][2]*x*x - lagrange[3][3]*y + lagrange[3][4]*y*y + lagrange[3][5]*x*y
        fxy = euclideanNorm(1/2,0)
    elif i == 4: 
        li = lagrange[4][0] + lagrange[4][1]*x + lagrange[4][2]*x*x - lagrange[4][3]*y + lagrange[4][4]*y*y + lagrange[4][5]*x*y
        fxy = euclideanNorm(0,1/2)
    elif i == 5: 
        li = lagrange[5][0] + lagrange[5][1]*x + lagrange[5][2]*x*x - lagrange[5][3]*y + lagrange[5][4]*y*y + lagrange[5][5]*x*y
        fxy = euclideanNorm(1/2,1/2)
    return fxy * li


def calcLagrange(xy0, xy1, xy2, xy3, xy4, xy5):
    matrix = np.array([
        [1, xy0[0], xy0[0]*xy0[0], xy0[1], xy0[1]*xy0[1], xy0[0]*xy0[1]], 
        [1, xy1[0], xy1[0]*xy1[0], xy1[1], xy1[1]*xy1[1], xy1[0]*xy1[1]], 
        [1, xy2[0], xy2[0]*xy2[0], xy2[1], xy2[1]*xy2[1], xy2[0]*xy2[1]], 
        [1, xy3[0], xy3[0]*xy3[0], xy3[1], xy3[1]*xy3[1], xy3[0]*xy3[1]], 
        [1, xy4[0], xy4[0]*xy4[0], xy4[1], xy4[1]*xy4[1], xy4[0]*xy4[1]], 
        [1, xy5[0], xy5[0]*xy5[0], xy5[1], xy5[1]*xy5[1], xy5[0]*xy5[1]]])
    
    l0 = np.array([1,0,0,0,0,0])
    l1 = np.array([0,1,0,0,0,0])
    l2 = np.array([0,0,1,0,0,0])
    l3 = np.array([0,0,0,1,0,0])
    l4 = np.array([0,0,0,0,1,0])
    l5 = np.array([0,0,0,0,0,1])
    
    X0 = np.linalg.inv(matrix).dot(l0)
    X1 = np.linalg.inv(matrix).dot(l1)
    X2 = np.linalg.inv(matrix).dot(l2)
    X3 = np.linalg.inv(matrix).dot(l3)
    X4 = np.linalg.inv(matrix).dot(l4)
    X5 = np.linalg.inv(matrix).dot(l5)
    return X0, X1, X2, X3, X4, X5


def euclideanNorm(x, y):
    return euclidean(x, y)
    return sqrt(x*x + y*y)


main()