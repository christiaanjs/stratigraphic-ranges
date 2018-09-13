package beast.evolution.tree;

import java.util.TreeMap;

public class SRMixedNode extends Node {

    @Override
    public int sort()  {
        SRMixedTree tree = (SRMixedTree) getTree();
        if(tree.getNodeIsSymmetric(this)){
            return super.sort();
        } else {
            throw new RuntimeException("Do not sort asymmetric nodes.");
        }

    }

    /**
     * @return (deep) copy of node
     */
    public SRMixedNode copy() {
        final SRMixedNode node = new SRMixedNode();
        node.height = height;
        node.labelNr = labelNr;
        node.metaDataString = metaDataString;
        node.metaData = new TreeMap<>(metaData);
        node.parent = null;
        node.setID(getID());

        for (final Node child : getChildren()) {
            node.addChild(child.copy());
        }
        return node;
    }

    @Override
    public void setLeft(final Node leftChild) {
        super.setLeft(leftChild);
        leftChild.setParent(this, false);
    }

    @Override
    public void setRight(final Node rightChild) {
        super.setRight(rightChild);
        rightChild.setParent(this, false);
    }

}
