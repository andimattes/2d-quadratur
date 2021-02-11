from scipy.integrate import quadpack
from scipy.integrate.quadpack import quad  
from math import sqrt, cos, sin
import numpy as np

x0, y0 = 2, 4
x1, y1 = 4, 5
x2, y2 = 5, 3

# Streckenmittelpunkte
x01, y01 = (x0+x1)/2, (y0+y1)/2
x02, y02 = (x0+x2)/2, (y0+y2)/2
x12, y12 = (x1+x2)/2, (y1+y2)/2

matrix =  np.array([[x1 - x0, x2 - x0], [y1 - y0, y2 - y0]])
det = 0

def main():
    arr = [[0, 0], [0, 1], [1, 0], [1/2, 0], [1/2, 1/2], [0, 1/2]] 
    
    pts = [] 
    for e, n in arr: 
        pnt = generalT(e, n)
        pts.append([pnt[0][0], pnt[1][0]]) 
       
    global det
    det = np.linalg.det(matrix)
        
    res1 = 0
    res2 = 0
    res3 = 0
    
    lagrange = calcLagrange((x0, y0), (x1, y1), (x2, y2), (x01, y01), (x02, y02), (x12, y12))
    for i in range(6):
        alpha, err = quad(outerIntegral, a=0, b=1, args=(i, lagrange))
        
        print("alpha #" + str(i) + ": " + '{:.2f}'.format(alpha))
        
        res1 += alpha * fxy1(pts[i][0], pts[i][1]) 
        res2 += alpha * fxy2(pts[i][0], pts[i][1]) 
        res3 += alpha * fxy3(pts[i][0], pts[i][1]) 
        print()

    print("Ergebnis0: " + '{:.4f}'.format(res1)) 
    print("Ergebnis1: " + '{:.4f}'.format(res2)) 
    print("Ergebnis2: " + '{:.4f}'.format(res3)) 
        
        
def outerIntegral(x, i, lagrange):
    ans, abserr = quad(innerIntegral, a=0, b=1-x, args=(x, i, lagrange))
    return ans


def fxy1(x, y):
    fxy = pow(x, 2) + pow(y, 2)
    print('fxy0: ' + str(fxy))
    return fxy
    
def fxy2(x, y):
    fxy = - 4*x - 2*y + 3 + 2*x*x - 5*x*y
    print('fxy1: ' + str(fxy))
    return fxy

def fxy3(x, y):
    fxy = sqrt(pow((x+1)*(y+1),2)-1)
    print('fxy2: ' + str(fxy))
    return fxy


def innerIntegral(y, x, i, lagrange):    
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
    return li * abs(det)
    

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
    

def generalT(e, n):    
    t1 = np.array([[x0], [y0]])
    t2 = matrix
    t3 = np.array([[e], [n]])
    
    x = t1 + np.dot(t2, t3)
    return x


main()