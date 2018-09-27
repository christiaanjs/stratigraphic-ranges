package beast.evolution.tree;

import sranges.StratigraphicRange;

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

    public String toShortNewick(final boolean printInternalNodeNumbers) {
        final StringBuilder buf = new StringBuilder();

        if (!isLeaf()) {
            buf.append("(");
            boolean isFirst = true;
            for (Node child : getChildren()) {
                if (isFirst)
                    isFirst = false;
                else
                    buf.append(",");
                buf.append(child.toShortNewick(printInternalNodeNumbers));
            }
            buf.append(")");
        }

        if (isLeaf() || printInternalNodeNumbers) {
            buf.append(getNr());
        }
        if(getID() != null){
            buf.append("/").append(getID());
        }
        SRTree tree = (SRTree) getTree();
        StratigraphicRange range = tree.getRangeOfNode(this);
        if(range != null){
            buf.append("/").append(range.getID());
        }

        buf.append(getNewickMetaData());
        buf.append(":").append(getNewickLengthMetaData()).append(getLength());
        return buf.toString();
    }


}
