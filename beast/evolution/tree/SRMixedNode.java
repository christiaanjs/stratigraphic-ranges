package beast.evolution.tree;

import java.util.TreeMap;

public class SRMixedNode extends Node {

    @Override
    public int sort()  {
        throw new RuntimeException("Do not sort potentially asymmetric nodes.");
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

    public void setLeft(final Node leftChild, boolean inOperator) {
        super.setLeft(leftChild);
        leftChild.setParent(this, inOperator);
    }

    @Override
    public void setRight(final Node rightChild) {
        super.setRight(rightChild);
        rightChild.setParent(this, false);
    }

    public void setRight(final Node rightChild, boolean inOperator) {
        super.setRight(rightChild);
        rightChild.setParent(this, inOperator);
    }


}
