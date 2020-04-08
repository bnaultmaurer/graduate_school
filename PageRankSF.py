# Page Rank Algorithm for Transit Network
# To implement on Sioux Falls Network
# Author: Benjamin Nault-Maurer
# Date: Winter 2016-2017

import time

# initialize network from file
def initNetwork(netFile):
# Initalize nodeList, outLinks, and link
    i = 0
    for b in netFile:
        net = b.split("\t")
        if i > 7:
            netList[i] = net
            upnode = float(netList[i][1])
            downnode = float(netList[i][2])

            # if the in node isn't in the node list, add it
            if upnode not in nodeList:
                nodeList.append(upnode)
            if downnode not in nodeList:
                nodeList.append(downnode)

            # if the in node isn't in the list of outgoing links, add it
            # this makes up a list of outgoing nodes from a given node (therefore you can call link[upnode, downnnode] later)
            if upnode in outNode:
                outNode[upnode] = outNode[upnode] + [downnode]
            else:
                outNode[upnode] = [downnode]
        i += 1

def PR(tmpPageRank, t):
    #print "\ntmpPageRank", tmpPageRank
    #print "t", t
    pageRank = {i: tmpPageRank[i] for i in tmpPageRank}
    if t == 0:
        pageRank = {i: 1/float(len(outNode[i])) for i in outNode}
        return pageRank
    else:
        for i in nodeList:
            pageRank[i] = 0
            #print "\noutNode", i , outNode[i]
            for j in outNode[i]:
                #print "len", tmpPageRank[j] / len(outNode[j])
                pageRank[i] += round(tmpPageRank[j]/float(len(outNode[j])),3)
                #print "pageRank[",i ,"]", pageRank[i]
    return pageRank


netFile = open("C:/Users/nault021/Google Drive/Research/PageRank/SiouxFalls/SiouxFalls_net.tntp",'r')
logFile = open('C:/Users/nault021/Google Drive/Research/PageRank/log.txt','w')
outFile = open('C:/Users/nault021/Google Drive/Research/PageRank/output.txt','w')

nodeList = []   # list of nodes in network
netList = {}    # temporary storage for network data from file
outNode = {}    # list of outLinks for each node
link = {}       # holds upstream node, downstream node, and travel time (in that order)
pageRank = {}   # hold the value of page rank for iteration
tmpPageRank = {}    # temp holding for pageRank

# start timer
start_time_1 = time.time()

initNetwork(netFile)

t = 0
while t < 100:
    pageRank = PR(pageRank, t)
    t += 1
'''
outFile.write("Node Importance of Sioux Falls Network after " + str(t) + " iterations\n\n" + "node\tpageRank\n")
for key, value in sorted(pageRank.iteritems(), key=lambda (k,v):(v,k)):
        outFile.write(str(key) + " \t" + str(value) + "\n")
'''

logFile.close()
outFile.close()
netFile.close()