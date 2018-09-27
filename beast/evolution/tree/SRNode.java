package beast.evolution.tree;

import sranges.StratigraphicRange;

import java.util.TreeMap;

/**
 * Created by gavryusa on 04/05/17.
 */
public class SRNode extends Node {

    /**
     * @return (deep) copy of node
     */
    public SRNode copy() {
        final SRNode node = new SRNode();
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
    } // copy

    @Override
    public int sort()  {
        throw new RuntimeException("Do not sort ordered trees. Calculation stopped.");
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
