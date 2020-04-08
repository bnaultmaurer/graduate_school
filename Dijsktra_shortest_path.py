# Code to find shortest path using Dijkstra's algorithm
# Script created for CEGE8217 FALL 2016 at the University of Minnesota
# Author: Benjamin Nault-Maurer

Q = {}           # holds current shortest distance for node and previous node
visited = []     # list of visited nodes to avoid repetitious calculations
SEL = {}         # scan eligible list
nodeList = []    # list of nodes in network
netList = {}     # temporary storage for network data from file
outLinks = {}    # list of outLinks for each node
link = {}        # holds upstream node, downstream node, and travel time (in that order)
OD = {}          # holds shortest distance from start node to destination
path = {}        # holds shortest path from start node to destination

# opening files
netFile = open("C:/Users/nault021/Desktop/SiouxFalls_net.tntp",'r')  
outFile = open("C:/Users/nault021/Desktop/ODMatrix_BenjaminNaultMaurer2.txt",'w')
     
# Initalize nodeList, outLinks, and link
i = 0
for b in netFile:
    net = b.split("\t")
    if i > 7:
        netList[i] = net
		
		# link[in node, out node] = [in node, out node, length (or travel time)]
        link[int(netList[i][1]),int(netList[i][2])] = [int(netList[i][1]),int(netList[i][2]),int(netList[i][4])]
		
		# if the in node isn't in the node list, add it
        if int(netList[i][1]) not in nodeList:
            nodeList.append(int(netList[i][1]))
			
		# if the in node isn't in the list of outgoing links, add it
		# this makes up a list of outgoing nodes from a given node (therefore you can call link[in node, out node] later)
        if int(netList[i][1]) in outLinks:
            outLinks[int(netList[i][1])] = outLinks[int(netList[i][1])] + [int(netList[i][2])]
        else:
            outLinks[int(netList[i][1])] = [int(netList[i][2])]
    i += 1

a = 0
for i in nodeList:
	# initialize values
    start = i
    visited = []
	
	# for each node, initialize OD matrix (for output), path values, and Q
    for a in nodeList:
        OD[i] = []
        path[i,a] = [i]
        Q[a] = [float("inf"),a]
		
		# if start node, set distance to 0
        if a == start:
            Q[start] = [0,start]
            SEL[start] = [0,start]
            OD[i,a] = 0
        a = a + 1

    # Dijkstra's algorithm, searching for shortest path
    while len(Q) != 0:
	
		# find the node with the shortest distance
        inode = min(SEL)
		
		# for each node in the outgoing nodes that haven't already been visited, 
		# check the distance from origin
        for onode in outLinks[inode]:
            if onode not in visited:
                SEL[onode] = [Q[onode][0],onode]
                alt = Q[inode][0] + link[inode,onode][2]
				
				# if the distance is shorter, reassign the distance
                if alt < Q[onode][0]:
                    Q[onode][0] = alt
                    SEL[onode] = [alt,inode]
                    OD[i,onode] = Q[onode][0]
                    path[i,onode] = path[i,inode] + [onode]
                    
        # prep for next iteration
        visited.append(inode)
        del Q[inode]
        del SEL[inode]

# write to output files
outFile.write("OD matrix for the Souix Falls network, created for CEGE8217 Fall 2016" + "\n")
outFile.write("Created by Benjamin Nault-Maurer" + "\n\n")
outFile.write("Travel Time:" + "\n")
outFile.write("(First column and row describe the origin and destination, respectively)" + "\n\n")

# write header line for OD matrix
outFile.write(str(0) + "\t")
for i in nodeList:
    outFile.write(str(i) + "\t")
outFile.write("\n")

# write out OD matrix
for i in nodeList:
    outFile.write(str(i) + "\t")
    for a in nodeList:
        outFile.write(str(OD[i,a]) + "\t")
    outFile.write("\n")

# write header line for shortest paths
outFile.write("\n\nShortest Paths" + "\n")
outFile.write("start node" + "\t" + "end node" + "\t" + "path" + "\n\n")
  
# write shortest paths to output file
for i in nodeList:
    for a in nodeList:
        outFile.write(str(i) + "\t" + str(a) + "\t" + str(path[i,a]) + "\n")

# close output file and import file
outFile.close()
netFile.close()
print "run completed"


