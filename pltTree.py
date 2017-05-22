
import sys
from time import time 
import matplotlib.pyplot as plt


def readEdges(fileEdges):
    ledges = []

    IN = open(fileEdges , "r")
    l = IN.readline()

    PROP = []

    while l != "":
        sl = l.strip().split()

        p1 = [int(i) for i in sl[0].split(",")]
        p2 = [int(i) for i in sl[1].split(",")]

        if len(sl)>2:
            PROP.append(int(sl[2]))
        else:
            PROP.append(1)

        ledges.append(( p1[0], p1[1] , p2[0], p2[1] ))

        l = IN.readline()

    IN.close()
    return ledges,PROP

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
    
    
    Edges,EdgeProp = readEdges(fileEdge)
    Points = readPoints(filePoint)
    
    t2 = time()
    print "reading ok (",t2 - t1,"sec)"
    print "plot..."
    
    
    

    Xs = [[],[]]
    Ys = [[],[]]
    Cs = ["lightgrey" , "black"]
    
    for i,e in enumerate(Edges):
    
        try:
    
            p1 = Points[e[0]][e[1]]
            p2 = Points[e[2]][e[3]]
    
            index = EdgeProp[i]

            Xs[index].append( p1[0] )
            Xs[index].append( p2[0] )
            Xs[index].append( None )
        
            Ys[index].append( p1[1] )
            Ys[index].append( p2[1] )
            Ys[index].append( None )
    
        except:
    #       continue
            print e
            print len(Points[e[0]]) , len(Points[e[2]])
    
        print i*100./len(Edges),"%","\r",
    
    fig = plt.figure(facecolor='white')
    ax = plt.axes(frameon=False)
    ax.axes.get_yaxis().set_visible(False)
    ax.axes.get_xaxis().set_visible(False)
    
    
    #for s in SEUIL_XTINCT:
    #    circle1 = plt.Circle((0, 0), s, color='r' , fill=False)
    #    ax.add_artist(circle1)
    
    ax.plot(Xs[0],Ys[0], color=Cs[0])
    ax.plot(Xs[1],Ys[1], color=Cs[1])
    
    
    t3 = time()
    print "took",t3 - t2,"s"
    plt.show()