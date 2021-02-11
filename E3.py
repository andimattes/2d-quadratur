from scipy.integrate import quad
import numpy as np
from math import sin, cos, sqrt

pts = (0, 0), (1, 0), (0, 1)

def main():
    res0 = 0
    res1 = 0
    res2 = 0
        
    lagrange = calcLagrange(pts)
    
    for i in range(3):
        alpha, abserr = quad(outerIntegral, a=0, b=1, args=(i, lagrange))
        
        print("alpha #" + str(i) + ": " + '{:.4f}'.format(alpha))
                
        res0 += alpha * fxy0(pts[i][0], pts[i][1]) 
        res1 += alpha * fxy1(pts[i][0], pts[i][1]) 
        res2 += alpha * fxy2(pts[i][0], pts[i][1]) 
        print()

    print("Ergebnis0: " + '{:.4f}'.format(res0)) 
    print("Ergebnis1: " + '{:.4f}'.format(res1)) 
    print("Ergebnis2: " + '{:.4f}'.format(res2)) 


def outerIntegral(x, i, lagrange): 
    ans, abserr = quad(innerIntegral, a=0, b=1-x, args=(x, i, lagrange))
    return ans


def fxy0(x, y):
    fxy = pow(x, 2) + pow(y, 2) 
    print('fxy0: ' + str(fxy))
    return fxy
    
def fxy1(x, y):
    fxy = sin(x) + cos(y) 
    print('fxy1: ' + str(fxy))
    return fxy

def fxy2(x, y):
    fxy = - 4*x - 2*y + 3 
    print('fxy2: ' + str(fxy))
    return fxy


def innerIntegral(y, x, i, lagrange):  
    li = 0
    if i == 0: 
        li = lagrange[0][0] + lagrange[0][1] * x + lagrange[0][2] * y
    elif i == 1: 
        li = lagrange[1][0] + lagrange[1][1] * x + lagrange[1][2] * y
    elif i == 2: 
        li = lagrange[2][0] + lagrange[2][1] * x + lagrange[2][2] * y
    return li


def calcLagrange(pts): 
    matrix = np.array([[1, pts[0][0], pts[0][1]], 
                       [1, pts[1][0], pts[1][1]], 
                       [1, pts[2][0], pts[2][1]]])
    
    l0 = np.array([1,0,0])
    l1 = np.array([0,1,0])
    l2 = np.array([0,0,1])
    
    X0 = np.linalg.inv(matrix).dot(l0)
    X1 = np.linalg.inv(matrix).dot(l1)
    X2 = np.linalg.inv(matrix).dot(l2)
    return X0, X1, X2


main()