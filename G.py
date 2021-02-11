import dmsh
import matplotlib.pyplot as plt
import F3


def main():
    fig = plt.figure()
    position = 231 
    
    geo = dmsh.Union( 
        [dmsh.Rectangle(-1.0, +0.5, -1.0, +0.5), dmsh.Rectangle(-0.5, +1.0, -0.5, +1.0)] 
    )
        
    for it in [1, 3, 5, 10, 30, 50]: 
        ax = fig.add_subplot(position) 
        drawPlot(it, ax, geo)
        position += 1
    
    plt.show() 
    
    
def drawPlot(it, ax, geo):
    edge_size = 0.25
    X, cells = dmsh.generate(geo, edge_size, show=False, max_steps=it, verbose=False) 

    eps = 1.0e-10
    is_inside = geo.dist(X.T) < eps 

    x, y, z = geo._get_xyz()
    ax.contour(x, y, z, levels=[0.0], colors="k")
    
    ax.plot(X[is_inside, 0], X[is_inside, 1], ".") 
    ax.triplot(X[:, 0], X[:, 1], cells)
    
    pts = []
    for x, y in zip(X[:, 0], X[:, 1]): 
        pts.append((float("{:.5f}".format(x)), float("{:.5f}".format(y))))
    
    triangles = []
    for cell in cells:         
        tria = [pts[cell[0]], pts[cell[1]], pts[cell[2]]] 
        triangles.append(tria) 

    res = 0
    for tria in triangles:
        erg = F3.main(tria[0][0], tria[0][1], tria[1][0], tria[1][1], tria[2][0], tria[2][1])  
        res += erg
    
    print("\n\n\nIteration #" + str(it) + ": " + '{:.4f}'.format(res))
    print("Number of triangles: " + str(len(triangles)))
    
    # Format
    ax.axis("square") 
    ax.title.set_text("Steps: " + str(it) + " | " + str(len(triangles)) + " triangles | " + " size: " + '{:.2f}'.format(res))
 

main()