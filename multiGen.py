# coding=utf-8

import time
import os
import sys
import requests
from io import StringIO
from termcolor import colored # For printing coloured text in terminal
from cdk_pywrapper.cdk_pywrapper import Compound
import math
import copy
import json

# Demo
# http://service.iris.edu/fdsnws/event/1/

###################################################################
###################################################################
###################################################################
###################################################################

class mainAnalysis():

###################################################################
###################################################################

    def __init__(self, obsPathwayQuery, predPathwayQuery):

        __author__ = "Jason Tam"
        __copyright__ = "Copyright 2020"
        __credits__ = ["Jason Tam"]
        __license__ = "GPL"
        __version__ = "1.0.1"
        __maintainer__ = "Jason Tam"
        __email__ = "tam@envipath.com"
        __status__ = "Development"

        self.obsPathwayQuery = obsPathwayQuery
        self.predPathwayQuery = predPathwayQuery
        self.queryParams = {'accept':'application/json'}

        return

###################################################################
###################################################################

    def runAnalysis(self):

        print(colored("Run analysis", "cyan"))
        
        obsPathwayJson, predPathwayJson = self.requestPathways()

        delta = "https://envipath.org/package/32de3cf4-e3e6-4168-956e-32fa5ddb0ce1/pathway/9b38552a-5c55-4032-9483-96615fa9c4b7/node/b603c102-386f-424a-97b9-66d00f032e9b"
        catechol = "https://envipath.org/package/32de3cf4-e3e6-4168-956e-32fa5ddb0ce1/pathway/9b38552a-5c55-4032-9483-96615fa9c4b7/node/02fcee4b-b9a0-4865-bb27-86bb5a8223ed"

        score = self.comparePathways(obsPathwayJson, predPathwayJson)
        print("Score: " + str(score))

        print("Success")

        return

###################################################################
###################################################################

    def comparePathways(self, obsPathwayJson, predPathwayJson):

        # Get Pathways with corrected Depth
        adjObsPathwayJson = self.getPathwayWithDepth(obsPathwayJson)
        adjPredPathwayJson = self.getPathwayWithDepth(predPathwayJson)

        # Get weights
        obsPathwayWeightDict = self.setPathwayEvalWeight(obsPathwayJson)
        predPathwayWeightDict = self.setPathwayEvalWeight(predPathwayJson)

        # Get Common Nodes
        predToObsCommonDict = self.getCommonNodes(obsPathwayJson, predPathwayJson)
        predToObsCommonDict_Loose = self.getCommonNodes(obsPathwayJson, predPathwayJson, fullInchiKey=False)

        # Get Unique Nodes
        uniquePredNodesIDList = self.getUniqueNodes(predPathwayJson, predToObsCommonDict)
        uniqueObsNodesIDList = self.getUniqueNodes(obsPathwayJson, predToObsCommonDict)

        startAndEndNodesIDs = self.findIntermediates(obsPathwayJson, predPathwayJson) # {intNode, {startNode, endNode}}

        print("predToObsCommonDict: " + str(len(predToObsCommonDict)))
        print("startAndEndNodesIDs: " + str(len(startAndEndNodesIDs)))
        print("uniquePredNodesIDList: " + str(len(uniquePredNodesIDList)))
        print("uniqueObsNodesIDList: " + str(len(uniqueObsNodesIDList)))

        intermediatesIDs = []
        for intNodeID in startAndEndNodesIDs:
            intermediatesIDs.append(intNodeID)

        if(len(intermediatesIDs)>0):
            predPathwayJson = self.getDepthAdjustedPathway(predPathwayJson, intermediatesIDs)

        # Initialize values to determine comparison score
        score_test_TP = 0.0
        score_pred_TP = 0.0
        score_FP = 0.0
        score_FN = 0.0
        finalScore = 0.0

        for commonNodeID in predToObsCommonDict_Loose:
            predNodeID = commonNodeID
            if (self.getDepthInPathway(predPathwayJson, predNodeID)>0):
                score_test_TP = score_test_TP + predPathwayWeightDict.get(predNodeID)

        for nodeID in uniqueObsNodesIDList:
            if (self.getDepthInPathway(obsPathwayJson, nodeID)>0):
                score_FN = score_FN + obsPathwayWeightDict.get(nodeID)

        for predNodeID in uniquePredNodesIDList:
            if (self.getDepthInPathway(predPathwayJson, predNodeID)>0):
                score_FP = score_FP + predPathwayWeightDict.get(predNodeID)

        if (score_test_TP + score_FP + score_FN) > 0:
            finalScore = score_test_TP / (score_test_TP + score_FP + score_FN)

        return finalScore

###################################################################
###################################################################

    def getDepthInPathway(self, pathwayJson, nodeID):

        depth = -99

        for node in pathwayJson['nodes']:
            if nodeID==node['id']:
                depth = node['depth']
                break

        return depth

###################################################################
###################################################################

    def getPathwayWithDepth(self, pathwayJson):

        outputPathwayJson = copy.deepcopy(pathwayJson)

        # Find root nodes
        rootNodes = []
        allNodes = []
        for x in pathwayJson['nodes']:
            if x['depth']==0:
                rootNodes.append(x['id'])
            allNodes.append(x['id'])

        currentDepth = 0

        nodeList = outputPathwayJson['nodes']

        upstreamNodes = rootNodes

        # Mark the root nodes as processed
        for upstreamNode in upstreamNodes:
            allNodes.remove(upstreamNode)

        while(len(allNodes)>0):
            newUpstreamNodes = []
            # loop through the upstreamNodes
            for upstreamNode in upstreamNodes:

                # find edge where upstreamNode is startNode
                for edge in pathwayJson['links']:
                    eductid = ""
                    productid = ""
                    for i in range(0,len(nodeList)):
                        # Match the indes from the edge json object
                        if i==edge['source']:
                            eductid=nodeList[i]['id']
                        elif i==edge['target']:
                            productid=nodeList[i]['id']
                    
                    # print("eductid: " + eductid)
                    # print("productid: " + productid)

                    if eductid==upstreamNode:
                        # Check if depth of product node is already processed
                        if productid in allNodes:
                            # Find the product in output pathway
                            for node in outputPathwayJson['nodes']:
                                if node['id']==productid:
                                    # Set new depth
                                    node['depth'] = currentDepth + 1
                                    # Mark processed
                                    allNodes.remove(productid)
                                    # Set product to be upstream nodes of next while iteration
                                    newUpstreamNodes.append(productid)

            # Set new depth for next while iteration
            currentDepth = currentDepth + 1
            # Set upstream nodes for next while iteration
            upstreamNodes = newUpstreamNodes

        return outputPathwayJson

###################################################################
###################################################################

    def getDepthAdjustedPathway(self, predPathwayJson, intermediatesIDs):

        outputPathwayJson = copy.deepcopy(predPathwayJson)

        rootNodesID = []

        for predNode in predPathwayJson['nodes']:
            if predNode['depth']==0:
                rootNodesID.append(predNode['id'])

        for predNode in outputPathwayJson['nodes']:
            if predNode['depth']<1:
                continue

            if predNode in intermediatesIDs:
                predNode['depth']==-99
            else:
                shortestPathsList = []
                maxSize = 0
                for rootNodeID in rootNodesID:
                    shortestPathIDList = self.getShortestPath(predPathwayJson, rootNodeID, predNode['id'])
                    if len(shortestPathIDList)>0:
                        shortestPathsList.append(shortestPathIDList)
                    if len(shortestPathIDList)>maxSize:
                        maxSize = len(shortestPathIDList)

                if len(shortestPathsList)>0:
                    # Get the shortest path from each of the paths to a root node
                    shortestIndex = 0
                    for i in range(0, len(shortestPathsList)):
                        shortestPathIDList = shortestPathsList[i]
                        if len(shortestPathIDList)<maxSize:
                            maxSize = len(shortestPathIDList)
                            shortestIndex = i

                    #  Determine the number of intermediates in this set of in between nodes
                    shortestPathIDList = shortestPathsList[shortestIndex]
                    numInts = 0
                    for nodeID in shortestPathIDList:
                        if nodeID in intermediatesIDs:
                            numInts = numInts + 1

                    # Adjust node depth
                    if numInts > 0:
                        diff = predNode['depth']-numInts
                        predNode['depth'] = diff

        return outputPathwayJson

###################################################################
###################################################################

    def getCommonNodes(self, obsPathwayJson, predPathwayJson, fullInchiKey=True):

        predToObsCommonDict = {}

        for x in obsPathwayJson['nodes']:
            for y in predPathwayJson['nodes']:
                if ((x['depth']>=0) and (y['depth']>=0)):
                    x_comp = Compound(compound_string=x['smiles'], identifier_type='smiles')
                    y_comp = Compound(compound_string=y['smiles'], identifier_type='smiles')
                    x_ikey = x_comp.get_inchi_key()
                    y_ikey = y_comp.get_inchi_key()

                    if (fullInchiKey==True):
                        if x_ikey==y_ikey:
                            predToObsCommonDict[y['id']] = x['id']
                    else:
                        x_ikey_front = x_ikey.split("-")[0]
                        y_ikey_front = y_ikey.split("-")[0]
                        if x_ikey_front==y_ikey_front:
                            predToObsCommonDict[y['id']] = x['id']

        return predToObsCommonDict

###################################################################
###################################################################

    def getUniqueNodes(self, pathwayJson, predToObsCommonDict):

        uniqueNodesIDList = []

        for x in pathwayJson['nodes']:
            if (x['id'] not in predToObsCommonDict.keys()) and (x['id'] not in predToObsCommonDict.values()):
                if x['depth']>=0:
                    uniqueNodesIDList.append(x['id'])

        return uniqueNodesIDList

###################################################################
###################################################################

    def setPathwayEvalWeight(self, pathwayJson):

        evalWeightDict = {}

        for x in pathwayJson['nodes']:
            if x['depth']>0:
                nodeID = x['id']
                depth = x['depth']
                evalWeightDict[nodeID] = 1 / math.pow(2,depth)

        return evalWeightDict


###################################################################
###################################################################

    def getDownStreamNodes(self, pathwayJson, nodeID):

        listOfDSNodesID = []

        # Get list of nodes from the pathway
        nodeList = pathwayJson['nodes']

        # find edge where upstreamNode is startNode
        for edge in pathwayJson['links']:
            eductid = ""
            productid = ""
            for i in range(0,len(nodeList)):
                # Match the indes from the edge json object
                if i==edge['source']:
                    eductid=nodeList[i]['id']
                elif i==edge['target']:
                    productid=nodeList[i]['id']

            # if educt ID matches node ID, downstream Node found
            if eductid==nodeID:
                listOfDSNodesID.append(productid)


        return listOfDSNodesID

###################################################################
###################################################################

    def getShortestPath(self, pathwayJson, startNodeID, endNodeID):
        
        # Test if endNodeID is reachable from startNodeID
        # Get set of downStreamNodes for each node
        allDownStreamNodes = {}
        for node in pathwayJson['nodes']:
            downstreamNodesIDList = self.getDownStreamNodes(pathwayJson, node['id'])
            allDownStreamNodes[node['id']] = downstreamNodesIDList

        # Get list of nodes from the pathway
        allNodeIDs = []
        for node in pathwayJson['nodes']:
            allNodeIDs.append(node['id'])

        pred = {}
        dist = {}

        # vector path stores the shortest path
        reversePath = []

        if (self.isReachable(allDownStreamNodes, allNodeIDs, startNodeID, endNodeID, pred, dist) == False):
            return reversePath

        # construct path
        crawl = endNodeID
        reversePath.append(crawl)
        
        while (len(pred[crawl])>0):
            reversePath.append(pred[crawl])
            crawl = pred[crawl]

        shortestPathIDList = []

        for i in reversed(range(1, (len(reversePath)-1))):
            if (self.getDepthInPathway(pathwayJson, reversePath[i])>=0):
                shortestPathIDList.append(reversePath[i])

        return shortestPathIDList

###################################################################
###################################################################

    def isReachable(self, allDownStreamNodes, allNodeIDs, startNodeID, endNodeID, pred, dist):

        # a queue to maintain queue of vertices whose
        # adjacency list is to be scanned as per normal
        # DFS algorithm
        queue = []
    
        # boolean array visited[] which stores the
        # information whether ith vertex is reached
        # at least once in the Breadth first search
        # Also set the default values for dist and pred
        visited = {}
        for nodeID in allNodeIDs:
            visited[nodeID] = False
            dist[nodeID] = 1000000
            pred[nodeID] = ""
    
        # now source is first to be visited and
        # distance from source to itself should be 0
        visited[startNodeID] = True
        dist[startNodeID] = 0
        queue.append(startNodeID)
    
        # standard BFS algorithm
        while (len(queue) != 0):
            u = queue[0]
            queue.pop(0)
            for i in range(len(allDownStreamNodes.get(u))):
            
                if (visited[allDownStreamNodes.get(u)[i]] == False):
                    visited[allDownStreamNodes.get(u)[i]] = True
                    dist[allDownStreamNodes.get(u)[i]] = dist[u] + 1
                    pred[allDownStreamNodes.get(u)[i]] = u
                    queue.append(allDownStreamNodes.get(u)[i])
    
                    # We stop BFS when we find
                    # destination.
                    if (allDownStreamNodes.get(u)[i] == endNodeID):
                        return True

        return

###################################################################
###################################################################

    def findIntermediates(self, obsPathwayJson, predPathwayJson):

        predToObsCommonDict = self.getCommonNodes(obsPathwayJson, predPathwayJson)
        startAndEndNodesIDs = {}


        for x in predToObsCommonDict:
            predNodeID = x
            obsNodeID = predToObsCommonDict.get(x)

            downstreamNodesIDList = self.getDownStreamNodes(obsPathwayJson, obsNodeID)

            for downstreamNodeID in downstreamNodesIDList:
                # Check if downstream node is a common node
                for predToObs in predToObsCommonDict:
                    if (downstreamNodeID==predToObsCommonDict.get(predToObs)):
                        predDSNodeID = predToObs
                        listOfAllIntIDs = self.getShortestPath(predPathwayJson, predNodeID, predDSNodeID)

                        for intNodeID in listOfAllIntIDs:
                            startAndEndNodeIDs = {}
                            startAndEndNodeIDs[obsNodeID] = downstreamNodeID
                            startAndEndNodesIDs[intNodeID] = startAndEndNodeIDs
                            intName = ""
                            startNodeName = ""
                            endNodeName = ""
                            for node in predPathwayJson['nodes']:
                                if node['id']==intNodeID:
                                    intName = node['name']
                                if node['id']==obsNodeID:
                                    startNodeName = node['name']
                                if node['id']==downstreamNodeID:
                                    starendNodeNametNodeName = node['name']
                            print("Int: " + intName + " startNode: " + startNodeName + " endNode: " + endNodeName)

        return startAndEndNodesIDs

###################################################################
###################################################################

    def requestPathways(self):

        # ### Observed Pathway

        # # Send HTTMP request to query observed pathway data
        obsPathwayReq = self.requestData(self.obsPathwayQuery,self.queryParams)

        # Convert String to json
        obsPathwayJson = obsPathwayReq.json()

        ### Predicted Pathway

        # # # Send HTTMP request to query predicted pathway data
        predPathwayReq = self.requestData(self.predPathwayQuery,self.queryParams)

        # Convert String to json
        predPathwayJson = predPathwayReq.json()

        return obsPathwayJson, predPathwayJson

###################################################################
###################################################################

    def requestData(self, path, queryParams):

        print("Sending HTTP request for Pathway...")
        r = requests.get(path, headers=queryParams)

        if (r.status_code != 200):
            print("Status Code: " + str(r.status_code))
            print(colored("HTML request failed", "red"))
            sys.exit()
        else :
            print(colored("HTML request succeeded", "green"))
        
        print("")

        return r

###################################################################
###################################################################

if __name__ == "__main__":

    startTime = time.time()

    ### Observed Pathway (Swap URL for desired pathway)
    baseLink = "https://envipath.org/"
    # BBD
    packageID = "32de3cf4-e3e6-4168-956e-32fa5ddb0ce1"
    # Pathway "beta-1,2,3,4,5,6-Hexachlorocyclohexane (an/aerobic)"
    pathwayID = "9b38552a-5c55-4032-9483-96615fa9c4b7"

    # Construct link for query
    obsPathwayQuery = baseLink+"/package/"+packageID+"/pathway/"+pathwayID

    ### Predicted Pathway (Swap URL for desired pathway)
    predPathwayQuery = "http://localhost:8080/package/67a42bea-db61-404d-871b-0d06ba3eb858/pathway/48415fe8-0d00-4d01-bbc3-d85613673372"

    Analysis_Object = mainAnalysis(obsPathwayQuery, predPathwayQuery)

    Analysis_Object.runAnalysis()

    endTime = time.time()

    print("")
    print ("Time elapsed: " + repr(endTime-startTime))