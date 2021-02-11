from scipy.integrate.quadpack import quad
from math import cos, sin, sqrt
import numpy as np


x0, y0 = 2, 4 
x1, y1 = 4, 5
x2, y2 = 5, 3

matrix =  np.array([[x1 - x0, x2 - x0], [y1 - y0, y2 - y0]]) 
det = 0

def main(x0, y0, x1, y1, x2, y2):
    arr = [[0, 0], [0, 1], [1, 0]] 
    
    pts = [] 
    for e, n in arr:
        pnt = generalT(e, n) 
        pts.append([pnt[0][0], pnt[1][0]]) 
    
    global det 
    det = np.linalg.det(matrix) 

    res1 = 0
    res2 = 0
    res3 = 0
    
    lagrange = calcLagrange((x0, y0), (x1, y1), (x2, y2))
    
    for i in range(3):
        alpha, err = quad(outerIntegral, a=0, b=1, args=(i, lagrange))
        
        print("alpha #" + str(i) + ": " + '{:.2f}'.format(alpha))
        
        res1 += alpha * fxy1(pts[i][0], pts[i][1])
        res2 += alpha * fxy2(pts[i][0], pts[i][1]) 
        res3 += alpha * fxy3(pts[i][0], pts[i][1]) 
        print()
        
    print("Ergebnis0: " + '{:.4f}'.format(res1)) 
    print("Ergebnis1: " + '{:.4f}'.format(res2)) 
    print("Ergebnis2: " + '{:.4f}'.format(res3)) 
    return res2      #für g


def outerIntegral(x, i, lagrange):
    ans, abserr = quad(innerIntegral, a=0, b=1-x, args=(x, i, lagrange))
    return ans


def fxy1(x, y):
    fxy = pow(x, 2) + pow(y, 2) 
    print('fxy0: ' + str(fxy))
    return fxy
    
def fxy2(x, y):
    fxy = sin(x) + cos(y)
    print('fxy1: ' + str(fxy))
    return fxy

def fxy3(x, y):
    fxy = sqrt(pow((x+1)*(y+1),2)-1)
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
    
    return li * abs(det)


def calcLagrange(xy0, xy1, xy2): 
    mtrx = np.array([[1, xy0[0], xy0[1]], [1, xy1[0], xy1[1]], [1, xy2[0], xy2[1]]])
    
    l0 = np.array([1,0,0])
    l1 = np.array([0,1,0])
    l2 = np.array([0,0,1])
    
    X0 = np.linalg.inv(mtrx).dot(l0)
    X1 = np.linalg.inv(mtrx).dot(l1)
    X2 = np.linalg.inv(mtrx).dot(l2)
    return X0, X1, X2
    

def generalT(e, n):    
    t1 = np.array([[x0], [y0]])
    t2 = matrix                
    t3 = np.array([[e], [n]])  
    
    x = t1 + np.dot(t2, t3) 
    return x


main(x0, y0, x1, y1, x2, y2)