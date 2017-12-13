
python concentricTree2.py test 0
time python cropTree.py test.edges.txt 5
time python extractSampledtree.py test.edges.txt.cropped.5.txt 1
python pltTree.py test.edges.txt.cropped.5.txt.sampled.txt test.points.txt 
