import unittest,random
import nj as hw


class SimpleTests(unittest.TestCase):

    def testB(self):
        # bestPair
        node1,node2=hw.bestPair(hw.tc1L,hw.tc1D)
        abCase=(node1==('A', (), (), 0) and node2==('B', (), (), 0)) or (node2==('A', (), (), 0) and node1==('B', (), (), 0))
        cdCase=(node1==('C', (), (), 0) and node2==('D', (), (), 0)) or (node2==('C', (), (), 0) and node1==('D', (), (), 0))
        self.assertTrue(abCase or cdCase,msg="Problem in bestPair.")
        node1,node2=hw.bestPair(hw.tc2L,hw.tc2D)
        self.assertTrue((node1==('B', (), (), 0) and node2==('D', (), (), 0)) or (node2==('B', (), (), 0) and node1==('D', (), (), 0)),msg="Problem in bestPair.")

        
    def testC(self):
        # bestPair
        node1,node2=hw.bestPair(hw.tc3L,hw.tc3D)
        self.assertTrue((node1==('F', (), (), 0) and node2==('G', (), (), 0)) or (node2==('F', (), (), 0) and node1==('G', (), (), 0)),msg="Problem in bestPair.")
        node1,node2=hw.bestPair(hw.tc4L,hw.tc4D)
        self.assertTrue((node1==('A', (), (), 0) and node2==('B', (), (), 0)) or (node2==('A', (), (), 0) and node1==('B', (), (), 0)),msg="Problem in bestPair.")

    def testD(self):
        # mergeNodes
        nodeA=('anc', ('anc', ('A', (), (), 3), ('B', (), (), 4), 5), ('anc', ('C', (), (), 1), ('D', (), (), 4), 3), 0)
        nodeB=('anc', ('anc', ('D', (), (), 9), ('F', (), (), 2), 12), ('anc', ('G', (), (), 13), ('H', (), (), 17), 33), 0)
        sol=('anc', ('anc', ('anc', ('A', (), (), 3), ('B', (), (), 4), 5), ('anc', ('C', (), (), 1), ('D', (), (), 4), 3), 20), ('anc', ('anc', ('D', (), (), 9), ('F', (), (), 2), 12), ('anc', ('G', (), (), 13), ('H', (), (), 17), 33), 30), 0)
        self.assertTrue(hw.mergeNodes(nodeA,nodeB,20,30)==sol,msg="Problem in mergeNodes.")

    def testE(self):
        # updateDistances
        L=[('X',(),(),0), ('Y',(),(),0), ('Z',(),(),0)] 
        D={(('Z', (), (), 0), ('Y', (), (), 0)): 13, (('Z', (), (), 0), ('X', (), (), 0)): 12, (('X', (), (), 0), ('Z', (), (), 0)): 12, (('X', (), (), 0), ('Y', (), (), 0)): 5, (('Z', (), (), 0), ('Z', (), (), 0)): 0, (('Y', (), (), 0), ('Z', (), (), 0)): 13, (('X', (), (), 0), ('X', (), (), 0)): 0, (('Y', (), (), 0), ('X', (), (), 0)): 5, (('Y', (), (), 0), ('Y', (), (), 0)): 0}
        node1=('X', (), (), 0)
        node2=('Y', (), (), 0)
        newNode=('anc', ('X', (), (), 2), ('Y', (), (), 3), 0)
        hw.updateDistances(node1,node2,newNode,L,D)
        solD={(('anc', ('X', (), (), 2), ('Y', (), (), 3), 0), ('Z', (), (), 0)): 10.0, (('X', (), (), 0), ('X', (), (), 0)): 0, (('Y', (), (), 0), ('X', (), (), 0)): 5, (('Z', (), (), 0), ('X', (), (), 0)): 12, (('Y', (), (), 0), ('Z', (), (), 0)): 13, (('Y', (), (), 0), ('Y', (), (), 0)): 0, (('X', (), (), 0), ('Z', (), (), 0)): 12, (('X', (), (), 0), ('Y', (), (), 0)): 5, (('Z', (), (), 0), ('Y', (), (), 0)): 13, (('Z', (), (), 0), ('anc', ('X', (), (), 2), ('Y', (), (), 3), 0)): 10.0, (('Z', (), (), 0), ('Z', (), (), 0)): 0}
        ansL=sorted(D.items())
        solL=sorted(solD.items())
        self.assertTrue(len(ansL)==len(solL),msg="Problem in updateDistances.")
        for i in range(len(ansL)):
            ansK,ansV=ansL[i]
            solK,solV=solL[i]
            self.assertTrue(ansK==solK,msg="Problem in updateDistances.")
            self.assertTrue(round(ansV,3)==round(solV,3),msg="Problem in updateDistances.")
        
    def testF(self):
        # updateDistances
        L=[('A', (), (), 0), ('B', (), (), 0), ('C', (), (), 0), ('D', (), (), 0), ('E', (), (), 0)]
        D={(('A', (), (), 0), ('A', (), (), 0)): 0.0, (('B', (), (), 0), ('B', (), (), 0)): 0.0, (('E', (), (), 0), ('C', (), (), 0)): 1.4, (('B', (), (), 0), ('E', (), (), 0)): 1.55, (('C', (), (), 0), ('C', (), (), 0)): 0.0, (('A', (), (), 0), ('D', (), (), 0)): 0.85, (('E', (), (), 0), ('A', (), (), 0)): 1.55, (('C', (), (), 0), ('E', (), (), 0)): 1.4, (('C', (), (), 0), ('D', (), (), 0)): 0.7, (('D', (), (), 0), ('A', (), (), 0)): 0.85, (('C', (), (), 0), ('B', (), (), 0)): 0.55, (('A', (), (), 0), ('E', (), (), 0)): 1.55, (('A', (), (), 0), ('C', (), (), 0)): 0.55, (('B', (), (), 0), ('D', (), (), 0)): 0.85, (('A', (), (), 0), ('B', (), (), 0)): 0.1, (('D', (), (), 0), ('C', (), (), 0)): 0.7, (('D', (), (), 0), ('B', (), (), 0)): 0.85, (('C', (), (), 0), ('A', (), (), 0)): 0.55, (('B', (), (), 0), ('A', (), (), 0)): 0.1, (('D', (), (), 0), ('E', (), (), 0)): 1.5, (('E', (), (), 0), ('D', (), (), 0)): 1.5, (('E', (), (), 0), ('E', (), (), 0)): 0.0, (('B', (), (), 0), ('C', (), (), 0)): 0.55, (('D', (), (), 0), ('D', (), (), 0)): 0.0, (('E', (), (), 0), ('B', (), (), 0)): 1.55}
        hw.updateDistances(('A', (), (), 0),('B', (), (), 0),('anc', ('A', (), (), 0.05), ('B', (), (), 0.05), 0),L,D)
        solD={(('A', (), (), 0), ('A', (), (), 0)): 0.0, (('B', (), (), 0), ('B', (), (), 0)): 0.0, (('E', (), (), 0), ('C', (), (), 0)): 1.4, (('B', (), (), 0), ('E', (), (), 0)): 1.55, (('C', (), (), 0), ('C', (), (), 0)): 0.0, (('A', (), (), 0), ('D', (), (), 0)): 0.85, (('E', (), (), 0), ('A', (), (), 0)): 1.55, (('anc', ('A', (), (), 0.05), ('B', (), (), 0.05), 0), ('E', (), (), 0)): 1.5, (('C', (), (), 0), ('E', (), (), 0)): 1.4, (('C', (), (), 0), ('D', (), (), 0)): 0.7, (('anc', ('A', (), (), 0.05), ('B', (), (), 0.05), 0), ('D', (), (), 0)): 0.7999999999999999, (('D', (), (), 0), ('A', (), (), 0)): 0.85, (('C', (), (), 0), ('B', (), (), 0)): 0.55, (('A', (), (), 0), ('E', (), (), 0)): 1.55, (('A', (), (), 0), ('C', (), (), 0)): 0.55, (('B', (), (), 0), ('D', (), (), 0)): 0.85, (('A', (), (), 0), ('B', (), (), 0)): 0.1, (('D', (), (), 0), ('C', (), (), 0)): 0.7, (('anc', ('A', (), (), 0.05), ('B', (), (), 0.05), 0), ('C', (), (), 0)): 0.5, (('D', (), (), 0), ('B', (), (), 0)): 0.85, (('C', (), (), 0), ('A', (), (), 0)): 0.55, (('B', (), (), 0), ('A', (), (), 0)): 0.1, (('D', (), (), 0), ('anc', ('A', (), (), 0.05), ('B', (), (), 0.05), 0)): 0.7999999999999999, (('D', (), (), 0), ('E', (), (), 0)): 1.5, (('E', (), (), 0), ('D', (), (), 0)): 1.5, (('E', (), (), 0), ('E', (), (), 0)): 0.0, (('B', (), (), 0), ('C', (), (), 0)): 0.55, (('C', (), (), 0), ('anc', ('A', (), (), 0.05), ('B', (), (), 0.05), 0)): 0.5, (('E', (), (), 0), ('anc', ('A', (), (), 0.05), ('B', (), (), 0.05), 0)): 1.5, (('D', (), (), 0), ('D', (), (), 0)): 0.0, (('E', (), (), 0), ('B', (), (), 0)): 1.55}
        ansL=sorted(D.items())
        solL=sorted(solD.items())
        self.assertTrue(len(ansL)==len(solL),msg="Problem in updateDistances.")
        for i in range(len(ansL)):
            ansK,ansV=ansL[i]
            solK,solV=solL[i]
            self.assertTrue(ansK==solK,msg="Problem in updateDistances.")
            self.assertTrue(round(ansV,3)==round(solV,3),msg="Problem in updateDistances.")

    
    def testG(self):
        # nj, tc1
        sol=('anc', ('anc', ('C', (), (), 3.0), ('D', (), (), 7.0), 0.5), ('anc', ('E', (), (), 10.0), ('anc', ('A', (), (), 6.0), ('B', (), (), 2.0), 1.0), 0.5), 0)
        tcTree=hw.nj(hw.tc1L,hw.tc1D)
        self.assertTrue(sorted(leafList(sol))==sorted(leafList(tcTree)),msg="nj fails on test case 1. Tree doesn't have correct leaves.")
        self.assertTrue(tt(sol,tcTree),msg="nj fails on test case 1. Tree doesn't have correct topology.")
        self.assertTrue(round(branchSum(sol),4)==round(branchSum(tcTree),4),msg="nj fails on test case 1. Tree doesn't have correct branch lengths.")
        
    def testH(self):
        # nj, tc2
        sol=('anc', ('anc', ('A', (), (), 1.5), ('E', (), (), 1.0), 0.75), ('anc', ('C', (), (), 2.0), ('anc', ('B', (), (), 3.0000000000000004), ('D', (), (), 0.9999999999999996), 2.0), 0.75), 0)
        tcTree=hw.nj(hw.tc2L,hw.tc2D)
        self.assertTrue(sorted(leafList(sol))==sorted(leafList(tcTree)),msg="nj fails on test case 1. Tree doesn't have correct leaves.")
        self.assertTrue(tt(sol,tcTree),msg="nj fails on test case 1. Tree doesn't have correct topology.")
        self.assertTrue(round(branchSum(sol),4)==round(branchSum(tcTree),4),msg="nj fails on test case 1. Tree doesn't have correct branch lengths.")

    def testI(self):
        # nj, tc3
        sol=('anc', ('anc', ('C', (), (), 4.0), ('anc', ('A', (), (), 2.0), ('B', (), (), 2.0), 3.0), 0.5), ('anc', ('D', (), (), 6.0), ('anc', ('E', (), (), 8.0), ('anc', ('F', (), (), 4.0), ('G', (), (), 1.0), 7.0), 2.0), 0.5), 0)
        tcTree=hw.nj(hw.tc3L,hw.tc3D)
        self.assertTrue(sorted(leafList(sol))==sorted(leafList(tcTree)),msg="nj fails on test case 1. Tree doesn't have correct leaves.")
        self.assertTrue(tt(sol,tcTree),msg="nj fails on test case 1. Tree doesn't have correct topology.")
        self.assertTrue(round(branchSum(sol),4)==round(branchSum(tcTree),4),msg="nj fails on test case 1. Tree doesn't have correct branch lengths.")

            
    def testJ(self):
        # nj, tc4
        sol=('anc', ('anc', ('C', (), (), 1.0), ('anc', ('A', (), (), 1.833333333333333), ('B', (), (), 3.166666666666667), 7.5), 0.25), ('anc', ('D', (), (), 3.75), ('E', (), (), 3.25), 0.25), 0)
        tcTree=hw.nj(hw.tc4L,hw.tc4D)
        self.assertTrue(sorted(leafList(sol))==sorted(leafList(tcTree)),msg="nj fails on test case 1. Tree doesn't have correct leaves.")
        self.assertTrue(tt(sol,tcTree),msg="nj fails on test case 1. Tree doesn't have correct topology.")
        self.assertTrue(round(branchSum(sol),4)==round(branchSum(tcTree),4),msg="nj fails on test case 1. Tree doesn't have correct branch lengths.")

    
## support stuff

BITS = 32
def labelLeaves(numLeaves, labelLength):
    leafLabels = [] 
    for internal in range(numLeaves):
        leafLabels.append(randomLabel(labelLength))
    return leafLabels

def randomLabel(labelLength):
    label = ""
    for i in range(labelLength):
        label = label + random.choice(["0", "1"])
    return label

def labelInternalNodes(Tree, leafLabels):
    labelD = {}
    labelInternal(Tree, leafLabels, labelD, Tree[0][0], [])  # recursively label internal nodes beginning at node 0
    return labelD

def labelInternal(Tree, leafLabels, labelD, node, visited):
    visited = visited + [node]
    N = len(Tree)
    leaves = int((N+2)/2)
    if node < leaves:  # it's a leaf
        return leafLabels[node]
    else:
        neighborLabels = []
        for neighbor in Tree[node]:
            if neighbor not in visited:
                neighborLabels.append(labelInternal(Tree, leafLabels, labelD, neighbor, visited))
        newLabel = xor(neighborLabels)
        labelD[newLabel] = node
        return newLabel

def xor(stringList):
    output = ""
    length = len(stringList[0])
    for pos in range(length):
        count = 0
        for string in stringList:
            if string[pos] == "1": count += 1
        if count % 2 == 0: output = output + "0"
        else: output = output + "1"
    return output

def distance(Tree1, Tree2):
    '''Robinson-Foulds distance, expects trees in adjacency format.'''
    N = len(Tree1) # nodes in tree
    numLeaves = int((N+2)/2)
    leafLabels = labelLeaves(numLeaves, BITS)
    labelT1 = labelInternalNodes(Tree1, leafLabels)
    labelT2 = labelInternalNodes(Tree2, leafLabels)
    # counter1 counts T1 - T2
    counter = 0
    for T1key in labelT1.keys():
        if not T1key in labelT2:
            counter += 1
    return counter

def leafCount(Tree):
    '''Count leaves (not internal nodes).'''
    if Tree[1]==():
        return 1
    else:
        return leafCount(Tree[1]) + leafCount(Tree[2])

def leafList(Tree):
    '''Return the names of all the leaves as a list.'''
    if Tree[1]==():
        return [Tree[0]]
    else:
        return leafList(Tree[1]) + leafList(Tree[2])

def branchSum(Tree):
    '''Return sum of lengths of all branches.'''
    if Tree[1]==():
        return Tree[3]
    else:
        return Tree[3] + branchSum(Tree[1]) + branchSum(Tree[2])
    
    
def nodeNumTree(Tree,leafNamesD,currentNodeNum):
    '''Convert the node names in Tree to numbers. For leaves, use mapping in leafNamesD. For internal node, use currentNodeNum.'''
    if Tree[1]==():
        return (leafNamesD[Tree[0]],(),(),Tree[3]),currentNodeNum
    else:
        newLeft,currentNodeNum=nodeNumTree(Tree[1],leafNamesD,currentNodeNum)
        newRight,currentNodeNum=nodeNumTree(Tree[2],leafNamesD,currentNodeNum)
        return (currentNodeNum,newLeft,newRight,Tree[3]),currentNodeNum+1

def makeAdj(nodeTree,adjL):
    '''Takes tree with numbered nodes and converts it to adjacency format, putting it in adjL.'''
    if nodeTree[1]==():
        return
    else:
        makeAdj(nodeTree[1],adjL) # left
        makeAdj(nodeTree[2],adjL) # right

        # now add info for branch to left subtree
        proximal=nodeTree[0]
        distal=nodeTree[1][0]
        adjL[proximal].append(distal)
        adjL[distal].append(proximal)
        # now add info for branch to right subtree
        distal=nodeTree[2][0]
        adjL[proximal].append(distal)
        adjL[distal].append(proximal)

        
def tup2Adj(Tree,leafNamesD):
    '''Convert a tree in bio 52 4-tuple format to an adjacency tree.'''
    firstInternal=max(leafNamesD.values()) + 1
    nodeTree,maxNodeNum=nodeNumTree(Tree,leafNamesD,firstInternal)
    adjL=[[] for i in range(maxNodeNum)]
    makeAdj(nodeTree,adjL)
    return adjL

def rfdist(Tree1, Tree2):
    '''Robinson-Foulds distance, expects trees in bio 52 4-tuple format.'''
    leafNamesL=leafList(Tree1)
    leafNamesL2=leafList(Tree2)
    if sorted(leafNamesL2)!=sorted(leafNamesL): return float('inf')
    
    leafNamesD={} # mapping between names and node numbers for leaves
    for i in range(len(leafNamesL)):
        leafNamesD[leafNamesL[i]]=i

    Tree1Adj=tup2Adj(Tree1,leafNamesD)
    Tree2Adj=tup2Adj(Tree2,leafNamesD)
    #return Tree1Adj
    return distance(Tree1Adj,Tree2Adj)

def tt(Tree1,Tree2):
    '''Return True if these unrooted trees are equivalent (have RF
       distance of 0).'''
    return rfdist(Tree1, Tree2)==0
    
if __name__ == '__main__':
    unittest.main()
