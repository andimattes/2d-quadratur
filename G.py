import dmsh
import matplotlib.pyplot as plt
import F3


def main():
    generate([1, 3, 5, 10, 30, 50])
    

def generate(iterations):
    fig = plt.figure() # Leeres pyplot sheet
    index = 231 # 231 --> 2 spalten, 3 zeilen, 1 = links oben
    
    geo = dmsh.Union( #fertige Fkt von dmsh: 2 Quadrate werden 'übereinander gelegt'
        [dmsh.Rectangle(-1.0, +0.5, -1.0, +0.5), dmsh.Rectangle(-0.5, +1.0, -0.5, +1.0)] #fertige Fkt von dmsh: erstellt Rechteck mit Längen als Parameter
    )
        
    for it in iterations: #für jedes Element in iterations (1,10,25,50)
        ax = fig.add_subplot(index) # Subplot hinzufügen
        drawPlot(it, ax, geo)
        index += 1 # = nächstes feld im plot
    
    plt.show() # Anzeigen
    
    
def drawPlot(it, ax, geo):
    edge_size = 0.2 #0.25 #Größe der Dreiecke
    X, cells = dmsh.generate(geo, edge_size, show=False, max_steps=it, verbose=False) #X mit gesetzten Pkte, cells mit Dreiecken, erstellt quasi das gamze
    #X, cells = optimesh.optimize_points_cells(X, cells, "CVT (full)", 1.0e-10, 100) # high quality

    ''' PLOT ''' 
    eps = 1.0e-10 #epsilon-umgebung 
    is_inside = geo.dist(X.T) < eps #stetigkeit? ist sowieso gefordert durch stetig diffbar

    x, y, z = geo._get_xyz() #pkte von geo 
    ax.contour(x, y, z, levels=[0.0], colors="k") # Original shape (schwarzer Rand)
    
    ax.plot(X[is_inside, 0], X[is_inside, 1], ".")  #
    ax.triplot(X[:, 0], X[:, 1], cells)
    
    pts = []
    for x, y in zip(X[:, 0], X[:, 1]):
        pts.append((ff(x), ff(y)))    
    
    triangles = []
    for cell in cells:        
        obj = [pts[cell[0]], pts[cell[1]], pts[cell[2]]]
        triangles.append(obj)

    res = 0
    for tria in triangles:
        erg = F3.main(tria[0][0], tria[0][1], tria[1][0], tria[1][1], tria[2][0], tria[2][1])
        res += erg
    
    print("Iteration #" + str(it) + ": " + '{:.4f}'.format(res))
    print("Number of triangles: " + str(len(triangles)))
    print()
    print("------")

    # Format
    ax.axis("square") # Macht die Plots quadratisch
    ax.title.set_text("Steps: " + str(it) + " | " + str(len(triangles)) + " triangles | " + " size: " + '{:.2f}'.format(res))


def ff(x):
    return float("{:.5f}".format(x))

main()