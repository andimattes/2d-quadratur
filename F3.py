from scipy.integrate.quadpack import dblquad
from scipy.spatial.distance import euclidean  
from math import sqrt
import numpy as np


x0, y0 = 2, 4 # Pkt (X0,Y0) von bel Dreieck
x1, y1 = 4, 5
x2, y2 = 5, 3

matrix =  np.array([[x1 - x0, x2 - x0], [y1 - y0, y2 - y0]]) # Matrix aus affiner Abbildung
det = 0

def main(x0, y0, x1, y1, x2, y2):
    arr = [[0, 0], [0, 1], [1, 0]] # Punkte des Einheitsdreiecks
    
    pts = [] # leeres Array definieren
    for e, n in arr:
        pnt = generalT(e, n) # Punkte des Einheitsdreiecks
        pts.append([pnt[0][0], pnt[1][0]]) # append: Array an ein Array anhängen
    
    global det
    det = np.linalg.det(matrix) # Determinante berechnen (fertige Fkt von numpy)
    #print("Det: " + '{:.2f}'.format(det)) # Determinante ausgeben (2 Kommastellen genau)
    #print("Points: ", end="") # neuen Punkte des Dreiecks ausgeben
    #print(pts)
    newpts = [[int(x[0]*abs(det)), int(x[1]*abs(det))] for x in pts] # Betrag der Determinante mit den neuen Punkten multipliziern
    #print("New Points: ", end="") #neu multiplizierte Punkte ausgeben
    #print(newpts)   
    
    # Integrationsgrenzen bestimmen
    a = min([x0, x1, x2]) # äußere grenze unten
    b = max([x0, x1, x2]) # äußere grenze oben
    
    c = min([y0, y1, y2]) # innere grenze unten
    d = max([y0, y1, y2]) # innere grenze oben
    
    res = 0
    
    X0, X1, X2 = calcLagrange((x0, y0), (x1, y1), (x2, y2))
    lagrange = [X0, X1, X2]
    for i in range(3):
        ans, err = dblquad(integrand, a=a, b=b, gfun=c, hfun=d, args=[i, pts, lagrange])
        res += ans
        #print("i=" + str(i) + ": " + '{:.2f}'.format(ans))

    #print("Ergebnis: " + '{:.2f}'.format(res))     
    return res   


def integrand(y, x, i, pts, lagrange):
    fxy = euclideanNorm(pts[i][0], pts[i][1]) # x und y an der Position i "herausfiltern"
    li = 0
    if i == 0:
        li = lagrange[0][0] + lagrange[0][1] * x + lagrange[0][2] * y
    elif i == 1:
        li = lagrange[1][0] + lagrange[1][1] * x + lagrange[1][2] * y
    elif i == 2:
        li = lagrange[2][0] + lagrange[2][1] * x + lagrange[2][2] * y
    
    return li * fxy * abs(det)


def calcLagrange(xy0, xy1, xy2):
    mtrx = np.array([[1, xy0[0], xy0[1]], [1, xy1[0], xy1[1]], [1, xy2[0], xy2[1]]])
    
    l0 = np.array([1,0,0])
    l1 = np.array([0,1,0])
    l2 = np.array([0,0,1])
    
    X0 = np.linalg.inv(mtrx).dot(l0)
    X1 = np.linalg.inv(mtrx).dot(l1)
    X2 = np.linalg.inv(mtrx).dot(l2)
    return X0, X1, X2
    

def generalT(e, n):    # affine Abbildung
    t1 = np.array([[x0], [y0]]) # 1. Term
    t2 = matrix                 # 2. term
    t3 = np.array([[e], [n]])   # 3. Term
    
    x = t1 + np.dot(t2, t3) # dot ist Vektor-Matrix Multiplikation und dann wird anderer vektor addiert
    return x


def euclideanNorm(x, y):
    return euclidean(x, y)
    return sqrt(x*x + y*y)


main(x0, y0, x1, y1, x2, y2)