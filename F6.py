from scipy.integrate.quadpack import dblquad  
from scipy.spatial.distance import euclidean
from math import sqrt
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
    arr = [[0, 0], [0, 1], [1, 0], [1/2, 0], [1/2, 1/2], [0, 1/2]] # Einheitsdreieck mit 6 Stützstellen
    
    pts = [] # leeres Array
    for e, n in arr: # e und n sind Punkte vom Einheitsdreieck
        pnt = generalT(e, n)
        pts.append([pnt[0][0], pnt[1][0]]) 
       
    global det
    det = np.linalg.det(matrix)
    print("Det: " + '{:.2f}'.format(det))
    print("Points: ", end="")
    print(pts)
    print([x01, y01], [x02, y02], [x12, y12])
    # newpts = [[int(x[0]*abs(det)), int(x[1]*abs(det))] for x in pts]
    #print("New Points: ", end="")
    #print(newpts)   
    
    # Integrationsgrenzen bestimmen
    a = min([x0, x1, x2]) # äußere grenze unten
    b = max([x0, x1, x2]) # äußere grenze oben
    
    c = min([y0, y1, y2]) # innere grenze unten
    d = max([y0, y1, y2]) # innere grenze oben
        
    res = 0
    
    X0, X1, X2, X3, X4, X5 = calcLagrange((x0, y0), (x1, y1), (x2, y2), (x01, y01), (x02, y02), (x12, y12))
    lagrange = [X0, X1, X2, X3, X4, X5]
    for i in range(6):
        ans, err = dblquad(integrand, a=a, b=b, gfun=c, hfun=d, args=[i, pts, lagrange])
        res += ans
        print("i=" + str(i) + ": " + '{:.2f}'.format(ans))

    print("Ergebnis: " + '{:.2f}'.format(res))
        

def integrand(y, x, i, pts, lagrange):    
    li = 0 # lagrange basis mit 6 Punkten von allgemeinen Dreieck
    fxy = euclideanNorm(pts[i][0], pts[i][1])
    
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
    print(det)
    return fxy * li * abs(det)
    

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


def euclideanNorm(x, y):
    return euclidean(x, y)
    return sqrt(x*x + y*y)


main()