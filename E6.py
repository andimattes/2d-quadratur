import numpy as np
from scipy.integrate import quad
from math import cos, sin, sqrt

pts = (0, 0), (1, 0), (0, 1), (0.5, 0), (0, 0.5), (0.5, 0.5)

def main():
    res1 = 0
    res2 = 0
    res3 = 0
    
    lagrange = calcLagrange(pts) # Lagrange Basen berechnen       
    for i in range(6): # Laufvariable für i = 0 bis 5
        alpha, err = quad(outerIntegral, a=0, b=1, args=(i, lagrange)) # Integralfunktion, welche fkt, äußere schranken, innere schranken, Laufvariable (für benutzte fkt)
        
        print("alpha #" + str(i) + ": " + '{:.4f}'.format(abs(alpha))) # jedes Integral in Laufvariable ausgegeben
        
        res1 += alpha * fxy1(pts[i][0], pts[i][1]) # aufsummieren
        res2 += alpha * fxy2(pts[i][0], pts[i][1]) # aufsummieren
        res3 += alpha * fxy3(pts[i][0], pts[i][1]) # aufsummieren
        print()
        
    print("Ergebnis1: " + '{:.4f}'.format(res1)) # aufsummiertes Endergebnis
    print("Ergebnis2: " + '{:.4f}'.format(res2)) # aufsummiertes Endergebnis
    print("Ergebnis3: " + '{:.4f}'.format(res3)) # aufsummiertes Endergebnis


def outerIntegral(x, i, lagrange):
    ans, abserr = quad(innerIntegral, a=0, b=1-x, args=(x, i, lagrange))
    return ans


def fxy1(x, y):
    fxy = pow(x, 2) + pow(y, 2) 
    print('fxy1: ' + str(fxy))
    return fxy
    
def fxy2(x, y):
    fxy = sin(x) + cos(y) - sin(x*x) + y*y - x*y - 0.17
    print('fxy2: ' + str(fxy))
    return fxy

def fxy3(x, y):
    fxy = sqrt(pow((x+1)*(y+1),2)-1)
    print('fxy3: ' + str(fxy))
    return fxy


def innerIntegral(y, x, i, lagrange): # inneres vom Integral (6 Stützstellen)
    li = 0
    if i == 0: 
        li = lagrange[0][0] + lagrange[0][1]*x + lagrange[0][2]*x*x - lagrange[0][3]*y + lagrange[0][4]*y*y + lagrange[0][5]*x*y
    elif i == 1: 
        li = lagrange[1][0] + lagrange[1][1]*x + lagrange[1][2]*x*x - lagrange[1][3]*y + lagrange[1][4]*y*y + lagrange[1][5]*x*y
    elif i == 2: 
        li = lagrange[2][0] + lagrange[2][1]*x + lagrange[2][2]*x*x - lagrange[2][3]*y + lagrange[2][4]*y*y + lagrange[2][5]*x*y
    elif i == 3: 
        li = lagrange[3][0] + lagrange[3][1]*x + lagrange[3][2]*x*x - lagrange[3][3]*y + lagrange[3][4]*y*y + lagrange[3][5]*x*y
    elif i == 4: 
        li = lagrange[4][0] + lagrange[4][1]*x + lagrange[4][2]*x*x - lagrange[4][3]*y + lagrange[4][4]*y*y + lagrange[4][5]*x*y
    elif i == 5: 
        li = lagrange[5][0] + lagrange[5][1]*x + lagrange[5][2]*x*x - lagrange[5][3]*y + lagrange[5][4]*y*y + lagrange[5][5]*x*y
    return li


def calcLagrange(pts):
    xy0, xy1, xy2, xy3, xy4, xy5 = pts[0], pts[1], pts[2], pts[3], pts[4], pts[5]
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


main()