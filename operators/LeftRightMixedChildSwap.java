package operators;

import beast.evolution.tree.Node;
import beast.evolution.tree.SRMixedTree;
import beast.evolution.tree.SRTree;
import beast.util.Randomizer;

import java.util.ArrayList;

public class LeftRightMixedChildSwap extends SRMixedTreeOperator {
    @Override
    public void initAndValidate() {
    }

    /**
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {

        SRMixedTree tree = treeInput.get(this);

        // choose a random node avoiding root and leaves that are direct ancestors
        int nodeCount = tree.getNodeCount();

        Node node = null;

        for (int i=0; i<5; i++) {
            node = tree.getNode(Randomizer.nextInt(nodeCount));
            if (!node.isLeaf() && !node.isFake() && !tree.belongToSameSRange(node.getNr(), node.getLeft().getNr())) {
                break;
            }
            node = null;
        }


        if (node == null) {
            ArrayList<Integer> allowableNodeIndices = new ArrayList<Integer>();

            for (int index=0; index<nodeCount; index++) {
                Node candidateNode = tree.getNode(index);
                //the node is not a leaf or sampled ancestor, the node is not fake, non of its children belongs to the same srange as node, the node is asymmetric
                if (!candidateNode.isLeaf() &&
                        !candidateNode.isFake() &&
                        !tree.belongToSameSRange(index, candidateNode.getLeft().getNr()) &&
                        !tree.getNodeIsSymmetric(candidateNode))
                    allowableNodeIndices.add(index);
            }

            int allowableNodeCount = allowableNodeIndices.size();

            if (allowableNodeCount == 0) {
                return Double.NEGATIVE_INFINITY;
            }

            node=tree.getNode(allowableNodeIndices.get(Randomizer.nextInt(allowableNodeCount)));
        }


        Node left = node.getLeft();
        Node right = node.getRight();

        node.setLeft(right);
        node.setRight(left);

        return 0.0;
    }
}
