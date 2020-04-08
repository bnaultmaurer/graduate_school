# MSA Algorithm to solve AV Routing problem
# Written for Final Course Project of CEGE 8217 at the University of Minnesota
# Last Work Date: 12/15/2016 (ongoing version)
# Programmer: Benjamin Nault-Maurer

import time

# calculate relative gap for MSA assignment
def gap(link, pathTT, demand, linkPath):
    TSTT = float(0)
    SPTT = float(0)
    for linkit in link:
        upnode = int(linkit[0])
        downnode = int(linkit[1])
        if (upnode, downnode) in linkPath:
            SPTT += link[upnode,downnode][5] * demand[upnode,downnode]
        TSTT += link[upnode,downnode][5]*link[upnode,downnode][6]
    gap = abs(TSTT - SPTT) / SPTT
    retGap = [gap, TSTT, SPTT]
    return retGap


# initialize demand from OD file
def initDemand(demandFile):
    i = 0
    for line in demandFile:
        if i > 4:
            tmp = line.split("\t")
            last = int(len(tmp)) - 1
            tmp[last] = tmp[last].split("\n")
            tmp[last] = tmp[last][0]
            node = int(tmp[0])
            demandWork[node] = []
            demandHome[node] = []
            for a in tmp:
                demandWork[node].append(int(a))
                demandHome[node].append(int(a))
        i += 1
    return [demandWork, demandHome]


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
            capacity = float(netList[i][3])
            length = float(netList[i][4])
            ffTime = float(netList[i][5])
            trTime = float(netList[i][5])
            flow = float(0)

            link[upnode, downnode] = [upnode, downnode, length, capacity, ffTime, trTime, flow]
            totLink[upnode, downnode] = [upnode, downnode, length, capacity, ffTime, trTime, flow]
            homeLink[upnode, downnode] = [upnode, downnode, length, capacity, ffTime, trTime, flow]
            workLink[upnode, downnode] = [upnode, downnode, length, capacity, ffTime, trTime, flow]
            auxFlow[upnode,downnode] = 0

            # if the in node isn't in the node list, add it
            if upnode not in nodeList:
                nodeList.append(upnode)
            if downnode not in nodeList:
                nodeList.append(downnode)

            # if the in node isn't in the list of outgoing links, add it
            # this makes up a list of outgoing nodes from a given node (therefore you can call link[upnode, downnnode] later)
            if upnode in outLinks:
                outLinks[upnode] = outLinks[upnode] + [downnode]
            else:
                outLinks[upnode] = [downnode]
        i += 1


# calculate travel time using BPR function
def travelTime(flow, capacity, ffTime, totLinkFlow, length, flag, fuelCost):
    flow = float(flow)
    capacity = float(capacity)
    ffTime = float(ffTime)
    totLinkFlow = float(totLinkFlow)
    
    # to use the typical BPR function
    if flag == "human":
        time = ffTime * (1 + 0.15*((totLinkFlow + flow)/capacity)**(4))
        
    # to calculate travel time based on fuel and parking costs
    elif flag == "AV":
        travelTime = 60 * ffTime * (1 + 0.15*((totLinkFlow + flow)/capacity)**(4))
        fuel = fuelCons(length, travelTime, flow, fuelCost)/12.30  # VOT: $12.30/hr
        if fuel == 0:
            time = ffTime
        elif fuel != 0:
            time = fuel
        else:
            print "something wrong with travelTime on AV trip"
            time = 0
    return time


# calculate the shortest path using Dijkstra's algorithm
def shortestPath(start, end, link, outLinks, nodeList, linkPath, parkingNode):
    Q = {}              # holds current shortest distance for node and previous node
    visited = []        # list of visited nodes to avoid repetitious calculations
    SEL = {}            # scan eligible list for nodes to recalculate
    nodePath = []       # list of nodes from origin to destination
    tmplinkPath = []    # temporary path of links from origin to destination

    for i in nodeList:
        visited = []
        OD[i] = []
        Q[i] = [float("inf"),i]
        prev[i] = 0
        if i == start:
            Q[start] = [0,start]
            SEL[start] = [0,start]
            OD[i,start] = 0   
        if i == parkingNode:
            continue
    
    while len(SEL) != 0:
        b = float("inf")
        if flag == "human":
            for n in SEL:
                if b > SEL[n][0] and n != parkingNode:
                    b = SEL[n][0]
                    upnode = int(n)
        elif flag == "AV":
           for n in SEL:
                if b > SEL[n][0]:
                    b = SEL[n][0]
                    upnode = int(n)
        else:
            print "Something wrong in shortestPath"
            break
        if upnode == end:
            break
        for downnode in outLinks[upnode]:
            if downnode not in visited:
                alt = Q[upnode][0] + link[upnode,downnode][5]
            if alt < Q[downnode][0]:
                Q[downnode][0] = alt
                Q[downnode][1] = upnode
                SEL[downnode] = [alt,upnode]
        visited.append(upnode)
        del SEL[upnode]

    # create nodePath
    i = end
    nodePath = [Q[end][1]] + [end]
    while i != start:
        i = Q[i][1]
        if i not in nodePath:
            nodePath = [Q[i][1]] + nodePath

    # print out list of nodes and links of shortest path from origin to destination
    a = 0
    j = 0
    while a+1 < len(nodePath):
        j = nodePath[a]
        jnxt = nodePath[a+1]
        if (j, jnxt) not in linkPath:
            linkPath.append((j, jnxt))
        if (j, jnxt) not in tmplinkPath:
            tmplinkPath.append((j, jnxt))
        a += 1
    pathCost = [nodePath, linkPath, Q[end][0], tmplinkPath]
    return pathCost
    
    
# calculate fuel cost based on travel time
def fuelCons(length, trTime, flow, fuelCost):
    cost = flow*(length/36.44)*(14.58*(length/trTime)**(-0.6253))
    return cost
    
    
# MSA traffic assignment
def MSA(link, tolerance, demandFile, totLink, logFile, parkingNode, flag):
    
    k = 0
    gapk = 500
    
    # Perform MSA Algorithm
    while gapk > tolerance:
        if k > 100000:
            break
        k += 1
        if k % 5000 == 0:
            print flag, "iteration", k ,"started at time", round((time.time() - start_time_1),2), "sec", "with gap", gapk
        logFile.write("iteration: " + str(k) + "\n")
    
        # initialize auxFlow
        for linkit in link:
            upnode = int(linkit[0])
            downnode = int(linkit[1])
            auxFlow[upnode,downnode] = 0
            demand[upnode,downnode] = 0
    
        # Step 2: find shortest path and assign auxilary flow to shortest path
        linkPath = []
        for origin in demandFile:
            for dest in demandFile:
                if demandFile[origin][dest] != 0:
                    tmpDemand = demandFile[origin][dest]
                    pathTT = shortestPath(origin, dest, link, outLinks, nodeList, linkPath, parkingNode)            
                    logFile.write("shortest Path for OD: " + str((origin,dest)) + " is " + str(pathTT[0]) + "\ntravel time " + str(pathTT[2]) + "\n")
                    for linkit in link:
                        upnode = int(linkit[0])
                        downnode = int(linkit[1])
                        if (upnode, downnode) in pathTT[3]:
                            demand[upnode,downnode] += tmpDemand
                            auxFlow[upnode,downnode] += tmpDemand
                            
        for linkit in link:
            upnode = int(linkit[0])
            downnode = int(linkit[1])
            
            # Step 3: update link flows: x[i,j,k+1] = x[i,j,k] + (1/k)*(y[i,j,k] - x[i,j,k])
            link[upnode,downnode][6] = ((float(k)-1)/float(k))*(float(link[upnode,downnode][6])) + (1/float(k))*(float(auxFlow[upnode,downnode]))
                        
            # Step 4: Update link travel time: BPR Function
            # if link to parking node, set parking cost (length)
            if downnode == parkingNode:
                link[upnode,downnode][5] = link[upnode,downnode][2]
            else:
                link[upnode,downnode][5] = travelTime(link[upnode,downnode][6], link[upnode,downnode][3], link[upnode,downnode][4], totLink[upnode,downnode][6], link[upnode,downnode][2], flag, fuelCost)
            logFile.write("link " + str(upnode) + " " + str(downnode) + " flow: " + str(link[upnode,downnode][6]) + " trTime: " + str(link[upnode,downnode][5]) + "\n")
        
        # Step 5: Termination: if gap < tolerance, terminate.
        if k != 1:
            gapk = gap(link, pathTT, demand, linkPath)
            logFile.write("TSTT: " + str(round(gapk[1],2)) + " SPTT: " + str(round(gapk[2],2)) + "\n")
            gapk = gapk[0]
        logFile.write("Gap: " + str(gapk) + "\n\n") 
        gapFile.write("iteration: " + str(k) + " Gap: " + str(gapk) + "\n")
    return link
            

###############################################################################
########################### Main Program ######################################
###############################################################################

demandFileWork = open("C:/Users/nault021/Google Drive/Classes/1TransportationNetworks/Final Project/SiouxFalls/Sioux_Falls_Demand_Matrix_Work_5Times.txt",'r')
demandFileHome = open("C:/Users/nault021/Google Drive/Classes/1TransportationNetworks/Final Project/SiouxFalls/Sioux_Falls_Demand_Matrix_Home_AVParking_5Times.txt",'r')
netFile = open("C:/Users/nault021/Google Drive/Classes/1TransportationNetworks/Final Project/SiouxFalls/Sioux_Falls_net_3TimesParking.txt",'r')
logFileWork = open('C:/Users/nault021/Google Drive/Classes/1TransportationNetworks/Final Project/SiouxFalls/logWork.txt','w')
logFileHome = open('C:/Users/nault021/Google Drive/Classes/1TransportationNetworks/Final Project/SiouxFalls/logHome.txt','w')
results = open('C:/Users/nault021/Google Drive/Research/CubeAssignment/SiouxFalls/results/results_5Demand_3PCost.txt','w')
gapFile = open('C:/Users/nault021/Google Drive/Research/CubeAssignment/SiouxFalls/gap/gap_5Demand_3PCost.txt','w')

nodeList = []   # list of nodes in network
nodePath = []   # hold path from origin to destination (by node)
linkPath = []   # path from origin to destination (by link)
pathCost = []   # hold information about cost of path
netList = {}    # temporary storage for network data from file
outLinks = {}   # list of outLinks for each node
link = {}       # holds upstream node, downstream node, and travel time (in that order)
OD = {}         # DO I USE THIS?
demandWork = {} # demand to CBD for initial run
demandHome = {} # demand to home for final run
prev = {}       # holds shortest path from start node to destination
auxFlow = {}    # hold auxiliary flow for MSA assignment
demand = {}     # to hold demand for assignment
totLink = {}    # to hold total link flows
workLink = {}   # to hold link flows from work assignment
homeLink = {}   # to hold link flows from home assignment

# Step 1: Initialize
initNetwork(netFile)
initDemand(demandFileWork)
parkingNode = 25
fuelCost = 2
allFlow = 0

k = 0
tolerance = 10**(-3)

# start timer
start_time_1 = time.time()

for linkit in link:
    upnode = int(linkit[0])
    downnode = int(linkit[1])
    #print "Work pretotLink[",upnode,downnode,"][6]", totLink[upnode,downnode][6]
    allFlow += totLink[upnode,downnode][6]
    #totLink[upnode,downnode][5] += workLink[upnode,downnode][5]
    #print "Work posttotLink[",upnode,downnode,"][6]", totLink[upnode,downnode][6]
    #print "work travelTime[",upnode,downnode,"]", travelTime(workLink[upnode,downnode][6], workLink[upnode,downnode][2], workLink[upnode,downnode][4], 0)

print "Before work assignment allFlow =", allFlow

# MSA for to work trips
flag = "human"
gapFile.write(str(flag) + "\n\n")
workLink = MSA(workLink, tolerance, demandWork, totLink, logFileWork, parkingNode, flag)

for linkit in workLink:
    upnode = int(linkit[0])
    downnode = int(linkit[1])
    totLink[upnode,downnode][6] += workLink[upnode,downnode][6]

print "After work assignment allFLow =", allFlow
# MSA for home/parking trips
initDemand(demandFileHome)
# set flag as "AV" to account for parking, set as "human" to not allow parking node
flag = "AV"
gapFile.write(str(flag) + "\n\n")
#print "\n\ntotlink", totLink
homeLink = MSA(homeLink, tolerance, demandHome, totLink, logFileHome, parkingNode, flag)

# sum it all together
allFlow = 0
for linkit in link:
    upnode = int(linkit[0])
    downnode = int(linkit[1])
    #print "\naddition", totLink[upnode,downnode][6], " + ", homeLink[upnode,downnode][6], " = ", totLink[upnode,downnode][6] + homeLink[upnode,downnode][6]
    totLink[upnode,downnode][4] = homeLink[upnode,downnode][4]
    totLink[upnode,downnode][5] = homeLink[upnode,downnode][5]
    totLink[upnode,downnode][6] += homeLink[upnode,downnode][6]
    allFlow += totLink[upnode,downnode][6]
    #print "final totLink[",upnode,downnode,"][6]", totLink[upnode,downnode][6]
    #print "home travelTime[",upnode,downnode,"]", travelTime(totLink[upnode,downnode][6], totLink[upnode,downnode][2], totLink[upnode,downnode][4], 0)

print "After work assignment allFLow =", allFlow
#print "\n\nfinal totLink", totLink



###############################################################################
###################### Output, Write to Files #################################
###############################################################################

# initialize final output measures
VMT = 0  
VHT = 0
ffVHT = 0
cost = 0
parkingCost = 0
totalCost = 0
allFlow = 0
parkingFlow = 0
sumFuelCost = 0
z1 = 0
z2 = 0
z3 = 0
z4 = 0
VOT = 0.21   # $0.21/min ($12.30/hr) as standard VOT (USDOT 2014)
for linkit in totLink:
        upnode = int(linkit[0])
        downnode = int(linkit[1])
        capacity = float(totLink[upnode,downnode][3])
        length = float(totLink[upnode,downnode][2])
        ffTime = float(totLink[upnode,downnode][4])
        trTime = float(totLink[upnode,downnode][5])
        flow = float(totLink[upnode,downnode][6])
        
        allFlow += flow
        VMT += flow * length
        VHT += flow * trTime
        ffVHT += flow * ffTime
        sumFuelCost += fuelCons(length, trTime, flow, fuelCost)
        if downnode == parkingNode:
            parkingCost += VOT*(flow * length)
            parkingFlow += flow
            if upnode == 10 or upnode == 16 or upnode == 17:
                z1 += flow
            if upnode == 9 or upnode == 11 or upnode == 15 or upnode == 18 or upnode == 19 or upnode == 20:
                z2 += flow
            if upnode == 5 or upnode == 6 or upnode == 7 or upnode == 8 or upnode == 12 or upnode == 14 or upnode == 21 or upnode == 22 or upnode == 23 or upnode == 24:
                z3 += flow
            if upnode == 1 or upnode == 2 or upnode == 3 or upnode == 4 or upnode == 13:
                z4 += flow
        if parkingFlow == 0:
            parkingFlow = 1
            
totalCost = sumFuelCost + parkingCost

results.write("Overall Parameters\n")

results.write("Flow: " + str(round(allFlow,2)) + "\n")
results.write("VMT: " + str(round(VMT,2)) + "\n")
results.write("VMT per vehicle: " + str(round(VMT/allFlow,2)) + "\n\n")

results.write("VHT (hr): " + str(round(VHT/60,2)) + "\n")
results.write("VHT per vehicle (min): " + str(round(((VHT)/allFlow),2)) + "\n")
results.write("Free Flow VHT (min): " + str(round((ffVHT),2)) + "\n")
results.write("VHT/Free Flow (min): " + str(round(((VHT-ffVHT)/ffVHT),2)) + "\n")
results.write("VHT/Free Flow per vehicle (min): " + str(round(((VHT-ffVHT)/allFlow),2)) + "\n")
results.write("Avg Travel Time: Avg FFtime: " + str(round( (VHT/allFlow)/(ffVHT/allFlow) ,2)) + " : " + str(round((ffVHT/ffVHT),2)) + "\n\n")

results.write("Fuel Cost ($): " + str(round(sumFuelCost,2)) + "\n")
results.write("Avg Fuel Cost per vehicle ($): " + str(round(sumFuelCost/allFlow,2)) + "\n")
results.write("% vehicles parking outside CBD: " + str(round((parkingFlow/allFlow)*100,2)) + "\n")
results.write("Parking Cost ($): " + str(round(parkingCost,2)) + "\n")
results.write("Avg Parking Cost ($): " + str(round(parkingCost/parkingFlow,2)) + "\n")
results.write("Cost($): " + str(totalCost) + "\n")
results.write("Avg Cost per vehicle ($): " + str(round(totalCost/allFlow,2)) + "\n\n")

results.write("% Parking in Zone 1:" + str(round(100*(z1/parkingFlow),2)) + "\n")
results.write("% Parking in Zone 2:" + str(round(100*(z2/parkingFlow),2)) + "\n")
results.write("% Parking in Zone 3:" + str(round(100*(z3/parkingFlow),2)) + "\n")
results.write("% Parking in Zone 4:" + str(round(100*(z4/parkingFlow),2)) + "\n\n")


results.write("Link Flows\n")
for linkit in link:
        upnode = int(linkit[0])
        downnode = int(linkit[1])
        capacity = float(totLink[upnode,downnode][3])
        length = float(totLink[upnode,downnode][2])
        ffTime = float(totLink[upnode,downnode][4])
        trTime = float(totLink[upnode,downnode][5])
        flow = float(totLink[upnode,downnode][6])
        results.write("link[ " + str(upnode) + " " + str(downnode) + " ]: " + str(flow) + "\n")

results.write("\nTravel Times\n")
for linkit in link:
        upnode = int(linkit[0])
        downnode = int(linkit[1])
        capacity = float(totLink[upnode,downnode][3])
        length = float(totLink[upnode,downnode][2])
        ffTime = float(totLink[upnode,downnode][4])
        trTime = float(totLink[upnode,downnode][5])
        flow = float(totLink[upnode,downnode][6])
        results.write("link[ " + str(upnode) + " " + str(downnode) + " ]: " + str(round(trTime,2)) + "\n")
        
results.write("\nTravel Costs\n")
for linkit in link:
        upnode = int(linkit[0])
        downnode = int(linkit[1])
        capacity = float(totLink[upnode,downnode][3])
        length = float(totLink[upnode,downnode][2])
        ffTime = float(totLink[upnode,downnode][4])
        trTime = float(totLink[upnode,downnode][5])
        flow = float(totLink[upnode,downnode][6])
        sumFuelCost += fuelCons(length, trTime, flow, fuelCost)
        results.write("link[ " + str(upnode) + " " + str(downnode) + " ]: " + str(round(fuelCost,2)) + "\n")
        
results.write("\nParking Costs\n")
for linkit in link:
        upnode = int(linkit[0])
        downnode = int(linkit[1])
        capacity = float(totLink[upnode,downnode][3])
        length = float(totLink[upnode,downnode][2])
        ffTime = float(totLink[upnode,downnode][4])
        trTime = float(totLink[upnode,downnode][5])
        flow = float(totLink[upnode,downnode][6])
        if downnode == parkingNode:
            parkingCost = VOT * (flow * length)
            results.write("link[ " + str(upnode) + " " + str(downnode) + " ]: " + str(round(parkingCost,2)) + "\n")

sumFuelCost = 0
results.write("\nTotal Costs\n")
for linkit in link:
        upnode = int(linkit[0])
        downnode = int(linkit[1])
        capacity = float(totLink[upnode,downnode][3])
        length = float(totLink[upnode,downnode][2])
        ffTime = float(totLink[upnode,downnode][4])
        trTime = float(totLink[upnode,downnode][5])
        flow = float(totLink[upnode,downnode][6])
        sumFuelCost += fuelCons(length, trTime, flow, fuelCost)
        if downnode == parkingNode:
            parkingCost = VOT *(flow * length)
            totalCost = sumFuelCost + parkingCost
            results.write("link[ " + str(upnode) + " " + str(downnode) + " ]: " + str(round(totalCost,2)) + "\n")
        else:
            totalCost = sumFuelCost + parkingCost
            results.write("link[ " + str(upnode) + " " + str(downnode) + " ]: " + str(round(totalCost,2)) + "\n")

###############################################################################
###############################################################################
###############################################################################
###############################################################################

demandFileWork.close()
demandFileHome.close()
gapFile.close()
logFileWork.close()
logFileHome.close()
results.close()
print "\nrun completed"
print "total run time",round((time.time() - start_time_1),2), "sec", round((time.time() - start_time_1)/60,2), "min"