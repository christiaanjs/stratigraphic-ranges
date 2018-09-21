package operators;

import beast.core.Input;
import beast.core.Operator;
import beast.evolution.tree.Node;
import beast.evolution.tree.SRMixedTree;
import beast.evolution.tree.Tree;

abstract public class SRMixedTreeOperator extends Operator {

    final public Input<SRMixedTree> treeInput = new Input<>("tree", "beast.tree on which this operation is performed", Input.Validate.REQUIRED);
    final public Input<Boolean> markCladesInput = new Input<>("markclades", "Mark all ancestors of nodes changed by the operator as changed," +
            " up to the MRCA of all nodes changed by the operator.", false);


    protected SRMixedTree tree;

    @Override
    public void initAndValidate(){
        this.tree = treeInput.get();
    }


    /**
     * @param parent the parent
     * @param child  the child that you want the sister of
     * @return the other child of the given parent.
     */
    protected Node getOtherChild(final Node parent, final Node child) {
        if (parent.getLeft().getNr() == child.getNr()) {
            return parent.getRight();
        } else {
            return parent.getLeft();
        }
    }

    /**
     * replace child with another node
     *
     * @param node
     * @param child
     * @param replacement
     */
    public void replace(final Node node, final Node child, final Node replacement) {
        node.removeChild(child);
        node.addChild(replacement);
        node.makeDirty(Tree.IS_FILTHY);
        replacement.makeDirty(Tree.IS_FILTHY);
    }

}
