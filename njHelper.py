
def tupleTree2Newick(tree):
    '''Convert a four tuple based tree (root,left,right,branchLen) into a
newick formated string.'''
    if tree[1]==():
        return tree[0]+":"+str(tree[3])
    else:
        leftStr=tupleTree2Newick(tree[1])
        rightStr=tupleTree2Newick(tree[2])
        return "("+leftStr+","+rightStr+"):"+str(tree[3])

def writeTree(tree,fileName):
    '''Write tree to fileName (in newick format).'''
    f=open(fileName,"w")
    f.write("("+tupleTree2Newick(tree)+");")
    f.close()
