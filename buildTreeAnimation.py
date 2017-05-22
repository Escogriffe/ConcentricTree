
import sys
from time import time 
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def readEdges(fileEdges):
    ledges = []

    IN = open(fileEdges , "r")
    l = IN.readline()

    while l != "":
        sl = l.strip().split()

        p1 = [int(i) for i in sl[0].split(",")]
        p2 = [int(i) for i in sl[1].split(",")]

        ledges.append(( p1[0], p1[1] , p2[0], p2[1] ))

        l = IN.readline()

    IN.close()
    return ledges

def readPoints(filePoints):
    Points = []

    IN = open(filePoints , "r")
    l = IN.readline()

    while l != "":
        sl = l.strip().split()

        gen,index = sl[0].split(",")

        gen = int(gen)
        index = int(index)

        while len(Points)<= gen:
            Points.append([])

        while len(Points[gen])<= index:
            Points[gen].append( None )


        x = float(sl[1])
        y = float(sl[2])
        

        Points[gen][index] = (x,y)
        
        l = IN.readline()

    IN.close()
    return Points





if __name__ == "__main__":


    help = """
    This script plots a concentric tree
    It takes 2 arguments:
        - file containing the edges of the tree
        - file containing the points of the tree

    """
    if len(sys.argv) <3:
        print help
        exit(1)

    fileEdge = sys.argv[1]
    filePoint = sys.argv[2]
    
    
    t1 = time()
    print "reading..."
    
    
    Edges = readEdges(fileEdge)
    Points = readPoints(filePoint)
    
    t2 = time()
    print "reading ok (",t2 - t1,"sec)"
    print "generating frames..."

    maxGen = 0
    for i,e in enumerate(Edges):
        if maxGen < e[2] : 
            maxGen = e[2]
    

    fig = plt.figure(figsize = (100/2.53 , 100/2.53 ), facecolor='white' ,dpi=300)
    ax = plt.axes(frameon=False)
    ax.axes.get_yaxis().set_visible(False)
    ax.axes.get_xaxis().set_visible(False)


    DRAWS = []

    genPerFrame = 50 ## should take about 30sec

    for GENERATION in range( (maxGen/genPerFrame) + 1 ):

        Xs = []
        Ys = []
        
        for i,e in enumerate(Edges):
        

            p1 = Points[e[0]][e[1]]
            p2 = Points[e[2]][e[3]]
    
            if e[0] < (GENERATION * genPerFrame):

                Xs.append( p1[0] )
                Xs.append( p2[0] )
                Xs.append( None )
            
                Ys.append( p1[1] )
                Ys.append( p2[1] )
                Ys.append( None )
                
            print int( GENERATION*1000./(maxGen/genPerFrame) )/10.,"%\t(", int( i*1000./len(Edges) )/10. ,"%)\r",
        
    
    
    #for s in SEUIL_XTINCT:
    #    circle1 = plt.Circle((0, 0), s, color='r' , fill=False)
    #    ax.add_artist(circle1)
    
        DRAWS.append( plt.plot(Xs,Ys, color="black") )
    



    duration = 3 * 1000 ## milliseconds of video

    ani = animation.ArtistAnimation(fig, DRAWS, interval=duration / len(DRAWS), blit=True)##making the animation

    writer = animation.ImageMagickFileWriter()
    ani.save('BigTree.gif', writer=writer)

    t3 = time()
    print "took",t3 - t2,"s"
    print "written animation"
    #plt.show()