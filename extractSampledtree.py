import networkx as nx


def computeAgeLineage(G,n):
    """ compute the number of generation to the last (previous) speciation event """
    age = 0
    S = G.successors(n)
    P = G.predecessors(n)

    if len(S)>1:
        return 0

    elif len(P)==0: ## root of the tree
        return 0

    return 1 + computeAgeLineage(G,P[0])

def computePathLineage(G,n):
    """ compute the Path to the last (previous) speciation event (excluding this speciation) """

    S = G.successors(n)
    P = G.predecessors(n)

    if len(S)>1:
        return []

    elif len(P)==0: ## root of the tree
        return []

    return [n] + computePathLineage(G,P[0])



def readEdges(fileEdges):
    ledges = []

    IN = open(fileEdges , "r")
    l = IN.readline()

    while l != "":
        sl = l.strip().split()
        ledges.append((sl[0],sl[1]))
        l = IN.readline()

    IN.close()
    return ledges

def getNodeGeneration(n):
    return int( n.split(",")[0] )

def getMaxGeneration(G):
    maxGen = 0
    for n in G.nodes():
        g = getNodeGeneration(n)
        if g > maxGen:
            maxGen = g
    return maxGen

def writeEdges(filename,G , keep ):
    OUT = open(filename,"w")

    #print "keep" , keep

    for e in G.edges():

        sampled = int( G.edge[ e[0] ][ e[1] ]['sampled'] ) 

        if not keep:
            if sampled:
                OUT.write(e[0] + " " + e[1] + "\n")
        else:
            OUT.write(e[0] + " " + e[1] + " " + str(sampled) + "\n")

    OUT.close()


def yieldPredessors(G , n ):
    current = n
    pre = G.predecessors(current)
    while len(pre)==1:
        yield pre[0]
        current = pre[0]
        pre = G.predecessors(current)




### very crude implemntation

import sys

help = """ 
This script crop the branches of a concentricTree that have become extinct before present so that only the tree of extant lineages remain
2 arguments:
  - file containing the edges of the tree
  - 1 ( or anything other than 0 , actually ) if the edges should all remain but a third column should be added to the edge file ( 0 for unsampled, 1 for sampled ) 

exemple : python extractSampledTree.py BigTree.edges.txt 
The cropped edges are written in a new file which has the name of the input file + ".sampled.txt"
"""

if len(sys.argv) <3:
    print help
    exit(1)

try:

    fileEdges = sys.argv[1]
    keep = bool(int(sys.argv[2])) # all lineage that aren't older than that will be removed
    
    ledges = readEdges(fileEdges)
    G = nx.DiGraph()

    for e in ledges:
        G.add_edge( e[0] , e[1] , sampled = False)
    
    SpecialGen = [] ## generation we want to exclude from this cropping (the last generation for instance ... )
    
    
    N1 = G.number_of_edges()
    
    maxGen = getMaxGeneration(G)
    #print "maxGen" , maxGen
    
    for e in G.edges():
        if G.edge[ e[0] ][ e[1] ]['sampled']:
            continue

        if getNodeGeneration(e[1]) == maxGen:
            current = e[1]
            #print "found", e[1]

            for pre in yieldPredessors(G , e[1] ):

                if G.edge[ pre ][ current ]['sampled'] :
                    break

                G.edge[ pre ][ current ]['sampled'] = True
                #print pre , current
                current = pre

    outPutName = fileEdges + ".sampled.txt"
 
    writeEdges(outPutName,G , keep)
    print "written new edges in" , outPutName

except Exception as e:
    print help
    print e