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

def writeEdges(filename,G):
    OUT = open(filename,"w")

    for e in G.edges():

        OUT.write(e[0] + " " + e[1] + "\n")
    OUT.close()


### very crude implemntation

import sys

help = """ 
This script crop the branches of a concentricTree that become extinct before a given cut-off.
2 arguments:
  - file containing the edges of the tree
  - cut-off (integer , >0)

exemple : python cropTree.py BigTree.edges.txt 3
The cropped edges are written in a new file which has the name of the input file + ".cropped.cut-off.txt"
"""

if len(sys.argv) <3:
    print help
    exit(1)

try:

    fileEdges = sys.argv[1]
    AgeLimite = int(sys.argv[2]) # all lineage that aren't older than that will be removed
    
    ledges = readEdges(fileEdges)
    G = nx.DiGraph(ledges)
    
    SpecialGen = [] ## generation we want to exclude from this cropping (the last generation for instance ... )
    
    
    N1 = G.number_of_edges()
    
    
    CONTINUE = True

    while CONTINUE:
    
        toRemove = []
        
        for n in G.nodes():
            if len(G.successors(n))==0:
        
                gen = getNodeGeneration(n)
                if not gen in SpecialGen:
                    path = computePathLineage(G,n)
                    age = len(path)
        
                    if age <= AgeLimite:
                        toRemove += path 
        
        G.remove_nodes_from(toRemove)
    
        CONTINUE = (len(toRemove)!=0)
    
    
    print "removed " ,  N1 - G.number_of_edges() , "out of" , N1 , "edges (",(100*(N1 - G.number_of_edges() ) )/ N1,"%)"
    outPutName = fileEdges + ".cropped." + str(AgeLimite) + ".txt" 
    writeEdges(outPutName,G)
    print "written new edges in" , outPutName

except Exception as e:
    print help
    print e