from hivDist import *
from njHelper import *

def nodeSep(node,nodeL,distD):
    '''A measure of the separation between a node of interest and
    the other nodes.'''
    distSum=0
    for iternode in nodeL:
        if node != iternode:
            distSum+=distD[(node,iternode)]
    return(float(distSum)/(len(nodeL)-2))

def njMetric(node1,node2,nodeL,distD):
    '''Calculates the neighbor joining metric between two nodes.'''
    d=distD[(node1,node2)]
    s1=nodeSep(node1,nodeL,distD)
    s2=nodeSep(node2,nodeL,distD)
    return(d-s1-s2)

def branchLength(nodeA,nodeB,nodeL,distD):
    '''Takes two nodes we're planning to merge, nodeA and nodeB. Calculates the
    branch lengths from their common ancestor to each.'''
    dist=distD[(nodeA,nodeB)]
    sepA=nodeSep(nodeA,nodeL,distD)
    sepB=nodeSep(nodeB,nodeL,distD)
    branchA=0.5*(dist+(sepA-sepB))
    branchB=0.5*(dist+(sepB-sepA))
    return(branchA,branchB)

def bestPair(nodeL,distD):
    """ inputs: nodeL, a list of nodes 
    and distD, a dictionary of distances
    outputs: a tuple, with a pair of nodes from node list 
    with the minimum value for the neighbor-joining metric"""
    bestIndices = [0, 1]
    bestDist = njMetric(nodeL[0], nodeL[1], nodeL, distD)
    for i in range(len(nodeL)):
        for j in range(len(nodeL)):
            
            if i != j:
                newDist = njMetric(nodeL[i], nodeL[j], nodeL, distD)
                if newDist < bestDist:
                    bestIndices = [i, j]
                    bestDist = newDist
    return nodeL[bestIndices[0]], nodeL[bestIndices[1]]

def mergeNodes(nodeA,nodeB,branchLenA,branchLenB):
    """ inputs: nodeA and nodeB, tree nodes
    branchLen A and branchLen B,  lengths from this newly created node to each of its daughters
    outputs:  newly created node, string 'anc' in the name position
    0 in the branch length position
    """
    newnodeA = (nodeA[0], nodeA[1], nodeA[2], branchLenA)
    newnodeB = (nodeB[0], nodeB[1], nodeB[2], branchLenB)
    return ('anc', newnodeA, newnodeB, 0)


def updateDistances(nodeA,nodeB,newNode,nodeL,distD):
    """inputs: nodeA and nodeB, two nodes that were merged to make newNode
    newNode, a new node whose distance needs to be calculated
    nodeL, a list of current nodes without newNode
    distD, a dictionary of distances
    outputs: does not return anything, but updates distD"""

    newDist = 0
    for n in range(len(nodeL)):
        if nodeL[n] != nodeA and nodeL[n] != nodeB:
            newDist = 0.5* (distD[(nodeA, nodeL[n])] + distD[(nodeB, nodeL[n])] - distD[nodeA, nodeB])
            distD[(nodeL[n], newNode)] = newDist
            distD[(newNode, nodeL[n])] = newDist
    
    
def nj(nodeL,distD):
    """inputs: nodeL, a list of nodes
    distD, a distance dictionary
    outputs: a phylogenetic tree in tuple tree format
    """
    #nextPair = tuple()
    #nextLength = tuple()
    #newNode = tuple()
    while len(nodeL) > 2:
        node1, node2 = bestPair(nodeL, distD)
        nextLength = branchLength(node1, node2, nodeL, distD)
        newNode = mergeNodes(node1, node2, nextLength[0], nextLength[1])
        updateDistances(node1, node2, newNode,nodeL,distD)
        nodeL.remove(node1)
        nodeL.remove(node2)
        nodeL += [newNode]
    return terminate(nodeL, distD)


def terminate(nodeL,distD):
    """inputs: nodeL, a list of nodes
    distD, a dictionary of distances
    outputs: a final merged node made from the two nodes in nodeL 
    """
    dist=distD[(nodeL[0],nodeL[1])]
    branchA=0.5*dist
    branchB=0.5*dist
    return mergeNodes(nodeL[0], nodeL[1], branchA,branchB)