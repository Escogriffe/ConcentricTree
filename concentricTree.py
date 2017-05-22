#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt




def angleDistToXY(angle,distance):
	x = np.cos(angle) * distance
	y = np.sin(angle) * distance
	return x,y


def formatAngle(angle):
	formatedAngle = angle

	if angle < 0:
		formatedAngle = np.pi * 2 + angle

	while formatedAngle > np.pi*2:
		formatedAngle -= np.pi * 2

	return formatedAngle

def getAngleDist(a1 , a2):
	""" always a positive angle """

	if a2 > a1:
		return a2 - a1

	return a2 + (np.pi*2 - a1)



def drawNewAngleOffset(angle,var, bounds = [] , offset = 0.):

	newAngle =  offset + np.random.randn()*var

	#print "drawNewAngleOffset", offset , var, bounds,"->" , newAngle  

	if len(bounds)==2:

		lowerB =  bounds[0]
		upperB =  bounds[1] 

		##drawing sign
		newAngle = abs(newAngle - offset)
		space = upperB - lowerB
		r = np.random.rand()
		if r > ( upperB / space):
			newAngle *= -1.

		newAngle += offset

		ok = False

		while not ok:
			ok = True
			if newAngle >= upperB:
				ok=False
				newAngle =  upperB - (newAngle-upperB) 

			if newAngle <= lowerB:
				ok=False
				newAngle =  lowerB + (lowerB-newAngle) 

		#print newAngle , lowerB,upperB

	return newAngle

def drawNewAngle(angle,var, bounds = []):

	newAngle =  drawNewAngleOffset(angle,var, bounds = [])

	return formatAngle( newAngle + angle )

def makeXChildren(p, lowerB,upperB , VAR , nbCh, diffToFather = 0.):
	""" returns new angles and difference to Father """

	ch = [ drawNewAngleOffset(p, VAR  , ( lowerB , upperB ) ,diffToFather )  for i in range(nbCh)]
	ch.sort()

	return [ formatAngle( a+p ) for a in ch]  , ch

def insertXchildren(Edges, newPoints, fatherIndex, fatherGen , children):
	""" return list offset """

	for ch in children:
		newPoints.append( ch )
		Edges.append( (fatherGen, fatherIndex, fatherGen+1,len(newPoints)-1  ) )

	return len(children)-1


def shiftPoints(Edges, ToFather, newPoints, currentGen , offset = 1):
	""" shift points in new points """

	newToFather = {}

	for i in range(offset):## point offset
		newPoints.insert(0, newPoints.pop(-1))

	def getNewIndex(i,o,l):
		""" get new index given index i and offset o and size of newPoints l """
		return (i+o)%l

	for i,e in enumerate(Edges):
		
		if e[2] == currentGen:
			newIndex = getNewIndex( e[3] , offset , len(newPoints) )
			newToFather[ newIndex ] = ToFather[ e[3] ]
			
			#print e[3], "-(",offset,",",len(newPoints),")->" ,newIndex
			Edges[i] = ( e[0] , e[1] , e[2] , newIndex )


	ToFather = newToFather.copy()

	return



	
def write_tree(filename , Points, Edges):
	"""
	Takes:
		- filename (str): prefix of the name of the files were the points and edges will be written
	"""
	
	OUT_POINT = open(filename + ".points.txt","w")
	OUT_EDGE = open(filename + ".edges.txt","w")
   
	for generation in range(len(Points)):
		for i,p in enumerate( Points[generation] ):
			OUT_POINT.write(str(generation)  + "," + str(i) + " " + str(p[0]) + " " + str(p[1]) + "\n")
		
	for e in Edges:
		OUT_EDGE.write(str(e[0]) + "," + str(e[1]) + " " + str(e[2]) + "," + str(e[3]) + "\n")
	
	OUT_POINT.close()
	OUT_EDGE.close()

	return
	

## Bref en % du rayon total et en distance en cm depuis la périphérie d'un cercle de 50cm, les 5 extinctions nous donnent: 
## 1 - 444MA - 12.7% - 6.34cm 
## 2 - 367MA - 10.5% - 5.24cm
## 3 - 250MA - 7.1% - 3.57cm
## 4 - 210MA - 6% - 3cm
## 5 - 65MA - 1.86% - 0.93cm




##############################################
###### PARAMETERS ############################

mmPerGen = 0.5  ### number of millimeters per generations. 
nbGen = 500 / mmPerGen ### number of generation. Here set so that the final tree has a total radius of 500mm (1 meter of diameter)

EXPECTEDNBPTPERMM = 0.05 * mmPerGen ## determine how many "biological niche" there is per millimeter of circle perimeter (1 circle = 1 generation)
### currently 0.05 niche per mm --> should create about 1 point every 2 cm.


VAR=0.25 ## variance of the angle between a father and its children
GRANDFATHERLAG=0.5 ##DO NOT SET ABOVE 1 importance of the position of the grand-father -> determine if the father and children will follow a coherent curve 

if GRANDFATHERLAG>1.:
	print "you have set the parameter GRANDFATHERLAG above 1. Beware weird behavior."

## the angle of a children, relative to the angle of its father 
## is determined by a normal law 
## of mean = GRANDFATHERLAG * ( difference of angle between the father and the grandfather )
## of variance = VAR * chidlrenFacor
##                     where  chidlrenFacor = ( number of children )^CHDISPERSIONPOWER  so that the more children you have, the more dispersed they are
##
## and bounded by the right and the left neighbors 

reductionFactor = 0.8 ## factor of reduction of the space allowed between the left and the right neighbor. The smaller the parameter, the farther apart the lineages should be

CHDISPERSIONPOWER= 2. ## affect the dispersion of different children round their father.

VARNBCHILDREN = 1. ## variance of the number of children per lineage.
MAXNBCHILDREN = 5 ##maximum number of children

TREENAME = "BigTree" ## prefix of the output files



SEUIL_XTINCT = [ #]
		int( ( (1. - (12.7/100.) ) *nbGen  ) ),
		int( ( (1. - (10.5/100.) ) *nbGen  ) ),
		int( ( (1. - (7.1 /100.) ) *nbGen  ) ),
		int( ( (1. - ( 6  /100.) ) *nbGen  ) ),
		int( ( (1. - (1.86/100.) ) *nbGen  ) )] ## extinction thresholds


THIN = 10 ## extinction thinning. here it means that only 1 in 10 lineage are left alive (at least 2)


##############################################


currentPTS = [ 0.0 ]

Points = []
Edges = []



SURVIVE = [True for i in currentPTS] 

d = 0

from time import time
t1 = time()


ToFather = {}


Points.append([])
i = 0
while i<len(currentPTS):
	p = currentPTS[i]
	Points[-1].append( angleDistToXY( p,0 ) )
	i+=1



while d < nbGen:


	if d%100 == 0:
		print "generation", d, "(",len(currentPTS),"points)"

	newPoints = []
	newToFather = {}


	i = 0
	while i<len(currentPTS):

		if not SURVIVE[i]:
			i += 1
			continue


		pointIndex = i

		p = currentPTS[pointIndex]



		lowerB = -np.pi 
		upperB = np.pi
		if len(currentPTS)>1:
			U = i+1
			if U ==len(currentPTS):
				upperB =  getAngleDist(p ,newPoints[0] )
			else:
				upperB =  getAngleDist(p ,currentPTS[U] )
	
			if len(newPoints) == 0:
				lowerB = - getAngleDist( currentPTS[pointIndex-1], p) 
			else:
				lowerB = - getAngleDist( newPoints[-1], p)


		lowerB *= reductionFactor
		upperB *= reductionFactor

		space =  ( upperB - lowerB ) * (d+1)
		#print lowerB , upperB ,"->", space
		XpectedNbCh = space * EXPECTEDNBPTPERMM
		#print "->",XpectedNbCh,

		nbCh =  min(MAXNBCHILDREN, max(0 , int(XpectedNbCh + np.random.randn() * VARNBCHILDREN ) ) ) ##np.random.poisson(XpectedNbCh) )

		##avoid total extinction
		if len(newPoints) == 0:
			nbCh = max(1 , nbCh)

		##looking at my father's position
		#if ToFather.has_key(pointIndex):
		#	FatherAngle = Points [ d-1 ][ ToFather[pointIndex] ]
		#	print "father angle" , FatherAngle , "<>", pointIndex

		FatherDiff = ToFather.get(i,0.) * GRANDFATHERLAG

		children , childrenDiffToFather = makeXChildren(p, lowerB,upperB , (nbCh)**CHDISPERSIONPOWER * VAR/((d+1)) , nbCh , diffToFather = FatherDiff )
		for j,a in enumerate(childrenDiffToFather):
			newToFather[ len(newPoints) + j ] = a

		#print FatherDiff,"->" , childrenDiffToFather

		offset  = insertXchildren(Edges, newPoints, pointIndex, d , children)



		i += 1#(offset + 1)
	
	## shifting points index to avoid spiralling too much
	shiftPoints(Edges, newToFather, newPoints, d+1 , offset = len(newPoints)/2)

	ToFather = newToFather.copy()

	#print ToFather

	currentPTS = newPoints[:]

	#print 'generation',d,Points[-1]

	d +=1

	Points.append([])
	for p in currentPTS:
		Points[-1].append( angleDistToXY( p,d ) )

	SURVIVE = [True for i in currentPTS] 

	if d in SEUIL_XTINCT:
		## extinction
		nbSurvivor = max( 2 , len(currentPTS)/THIN )

		print "extinction at",d," nb of surviving lineages : ",nbSurvivor

		survivors = np.random.choice(range(len(currentPTS)), size=nbSurvivor, replace=False)


		SURVIVE = [False for i in currentPTS] 
		for i in survivors:
			SURVIVE[i] = True


t2 = time()
print "done (",t2-t1,"s)"

print sum([len(g) for g in Points]),"points" , len(Edges),"edges"

print "plot..."
t2 = time()

Xs = []
Ys = []

for i,e in enumerate(Edges):

	try:

		p1 = Points[e[0]][e[1]]
		p2 = Points[e[2]][e[3]]

		Xs.append( p1[0] )
		Xs.append( p2[0] )
		Xs.append( None )
	
		Ys.append( p1[1] )
		Ys.append( p2[1] )
		Ys.append( None )

	except:
#		continue
		print e
		print len(Points[e[0]]) , len(Points[e[2]])

	print i*100./len(Edges),"%","\r",

fig = plt.figure(facecolor='white')
ax = plt.axes(frameon=False)
ax.axes.get_yaxis().set_visible(False)
ax.axes.get_xaxis().set_visible(False)


for s in SEUIL_XTINCT:
	circle1 = plt.Circle((0, 0), s, color='r' , fill=False)
	ax.add_artist(circle1)

ax.plot(Xs,Ys, color="black")


t3 = time()
print "took",t3 - t2,"s"
plt.show()




write_tree( TREENAME , Points, Edges)

print "results written in ",TREENAME + ".*"