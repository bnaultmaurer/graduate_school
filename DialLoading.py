# Code to find shortest path using Dijkstra's algorithm
# Script created for CEGE8217 FALL 2016 at the University of Minnesota
# Author: Benjamin Nault-Maurer

def ShortestPath(start,stop,nodeList,outLinks,link):
    Q = {}           # holds current shortest distance for node and previous node
    visited = []     # list of visited nodes to avoid repetitious calculations
    SEL = {}        # scan eligible list
    dist = {}         # holds shortest distance from start node to destination
    path = {}       # holds shortest path from start node to destination
    
    a = 0
    visited = []
    print "nodeList is:",nodeList
    print "outLinks is:", outLinks
    print "link is",link
    for a in nodeList:
        dist[start] = []
        path[start,a] = [i]
        Q[a] = [float("inf"),a]
        if a == start:
            Q[start] = [0,start]
            SEL[start] = [0,start]
            dist[start,a] = 0
            a = a + 1
    
    # searching for shortest path
    while len(Q) != 0:
        inode = min(SEL)
        for onode in outLinks[inode]:
            if onode not in visited:
                SEL[onode] = [Q[onode][0],onode]
                alt = Q[inode][0] + link[inode,onode][2]
                if alt < Q[onode][0]:
                    Q[onode][0] = alt
                    SEL[onode] = [alt,inode]
                    dist[start,onode] = Q[onode][0]
                    print "dist[",start,onode,"]:",dist[start,onode]
                    path[start,onode] = path[start,inode] + [onode]
                
        # prep for next iteration
        visited.append(inode)
        del Q[inode]
        del SEL[inode]
    
        if min(SEL) == [] or onode == stop:
            break
    print "shortest distances are:", dist

###############################################################################
#Start Code
###############################################################################

nodeList = []  # list of nodes in network
outLinks = {} # list of outLinks for each node
link = {}        # holds upstream node, downstream node, and travel time (in that order)

# opening files
netFile = open("C:/Users/nault021/Google Drive/Classes/1TransportationNetworks/HW3/DialNetwork.txt",'r')  
outFile = open("C:/Users/nault021/Google Drive/Classes/1TransportationNetworks/HW3/DialOutput.txt",'w')
     
# Initalize nodeList, outLinks, and link
i = 0
for b in netFile:
    net = b.split("\t")
    if i > 0:
        inode = int(net[1])         # in node, or upstream node
        onode = int(net[2])         # out node, or downstream node
        net[4] = net[4].split("/n")
        travelTime = int(net[4][0])
        link[inode,onode] = [inode,onode,travelTime]
        if inode not in nodeList:
            nodeList.append(inode)
        if onode not in nodeList:
            nodeList.append(onode)
        if inode in outLinks:
            outLinks[inode] = outLinks[inode] + [onode]
        else:
            outLinks[inode] = [onode]
    i += 1 

ShortestPath(1,6,nodeList,outLinks,link)
ShortestPath(2,9,nodeList,outLinks,link)

'''
# write to output files
outFile.write("Output for Problem 1, HW3, created for CEGE8217 Fall 2016" + "\n")
outFile.write("Created by Benjamin Nault-Maurer" + "\n\n")

# write header line for shortest paths
outFile.write("\n\nShortest Paths" + "\n")
outFile.write("start node" + "\t" + "end node" + "\t" + "path" + "\n\n")
  
# write shortest paths to output file
for i in nodeList:
    for a in nodeList:
        outFile.write(str(i) + "\t" + str(a) + "\t" + str(path[i,a]) + "\n")
'''
# close output file and import file
outFile.close()
netFile.close()
print "run completed"


