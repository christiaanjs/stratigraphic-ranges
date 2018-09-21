package beast.evolution.tree;

import beast.core.StateNode;
import sranges.StratigraphicRange;

import java.util.List;

public class SRMixedTree extends SRTree {

    public SRMixedTree(){}

    public SRMixedTree(Node root, List<StratigraphicRange> sRanges, List<Node> symmetricNodes){
        super(root, sRanges);
        initSymmetric();
        for(Node symmetricNode: symmetricNodes){
            setNodeIsSymmetric(symmetricNode, true);
        }
    }

    private boolean[] nodeIsSymmetric;
    private boolean[] storedNodeIsSymmetric;

    @Override
    public void initAndValidate(){
        super.initAndValidate();
        initSymmetric();
    }

    protected void initSymmetric(){
        nodeIsSymmetric = new boolean[getNodeCount()];
        storedNodeIsSymmetric = new boolean[getNodeCount()];
    }

    private void storeSymmetric(){
        System.arraycopy(nodeIsSymmetric, 0, storedNodeIsSymmetric, 0, nodeIsSymmetric.length);
    }

    private void assignSymmetric(boolean[] otherSymmetric){
        System.arraycopy(otherSymmetric, 0, nodeIsSymmetric, 0, nodeIsSymmetric.length);
    }

    @Override
    public void assignFrom(StateNode other){
        super.assignFrom(other);
        final SRMixedTree tree = (SRMixedTree) other;
        initSymmetric();
        assignSymmetric(tree.nodeIsSymmetric);
    }

    @Override
    public void assignFromFragile(StateNode other){
        super.assignFromFragile(other);
        final SRMixedTree tree = (SRMixedTree) other;
        initSymmetric();
        assignSymmetric(tree.nodeIsSymmetric);

    }

    @Override
    public void store(){
        super.store();
        storeSymmetric();
    }

    @Override
    public void restore() {

        super.restore();

        boolean[] tmp = storedNodeIsSymmetric;
        storedNodeIsSymmetric = nodeIsSymmetric;
        nodeIsSymmetric = tmp;
    }

    public boolean getNodeIsSymmetric(Node node){
        if(node.isDirectAncestor() | node.isLeaf()){
            throw new RuntimeException("Only speciation event nodes have symmetry");
        }
        return nodeIsSymmetric[node.getNr()];
    }

    public void setNodeIsSymmetric(Node node, boolean symmetric){
        if(node.isFake() | node.isLeaf()){
            throw new RuntimeException("Only speciation event nodes have symmetry");
        }
        nodeIsSymmetric[node.getNr()] = symmetric;
    }

}
