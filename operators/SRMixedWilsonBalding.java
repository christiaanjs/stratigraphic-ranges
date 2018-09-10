package operators;

import beast.evolution.tree.Node;
import beast.evolution.tree.SRMixedTree;
import beast.evolution.tree.SRTree;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import sranges.StratigraphicRange;

import java.util.ArrayList;
import java.util.Random;

public class SRMixedWilsonBalding extends SRMixedTreeOperator {
    @Override
    public void initAndValidate() {
    }

    /**
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {

        SRMixedTree tree = treeInput.get(this);

        //double x0 = 10;

        double oldMinAge, newMinAge, newRange, oldRange, newAge, fHastingsRatio, dimensionCoefficient;
        int newDimension, oldDimension;

        // choose a random node avoiding root and leaves that are direct ancestors
        int nodeCount = tree.getNodeCount();

        ArrayList<Integer> allowableNodeIndices = new ArrayList<Integer>();
        ArrayList<Integer> sRangeInternalNodeNrs = tree.getSRangesInternalNodeNrs();

        // select child
        for (int index=0; index<nodeCount; index++) {
            Node node = tree.getNode(index);
            // the node is not the root,
            //  it is not a sampled ancestor on a zero branch, (we select fake node)
            //  it is not an internal node of a stratigraphic range
            if (!node.isRoot() &&
                    !node.isDirectAncestor() &&
                    !sRangeInternalNodeNrs.contains(node.getNr()))
                allowableNodeIndices.add(index);
        }

        Node i;

        int allowableNodeCount = allowableNodeIndices.size();

        if (allowableNodeCount == 0) {
            return Double.NEGATIVE_INFINITY;
        }

        i = tree.getNode(allowableNodeIndices.get(Randomizer.nextInt(allowableNodeCount)));

        Node iP = i.getParent();
        Node CiP; // Sibling of i
        if (iP.getLeft().getNr() == i.getNr()) {
            CiP = iP.getRight();
        } else {
            CiP = iP.getLeft();
        }

        // make sure that there is at least one candidate edge to attach node iP to
        // If iP is the root then i must be attached as a leaf below CiP
        if (iP.getParent() == null && CiP.getHeight() <= i.getHeight()) {
            return Double.NEGATIVE_INFINITY; //
        }

        // choose another random node to insert i above or to attach i to this node if it is a leaf
        Node j;
        Node jP;

        final int leafNodeCount = tree.getLeafNodeCount();

        if (leafNodeCount != tree.getExternalNodes().size()) {
            System.out.println("node counts are incorrect. NodeCount = " + nodeCount + " leafNodeCount = " +
                    leafNodeCount + " external node count = " + tree.getExternalNodes().size());
        }

        // make sure that the target branch <jP, j> or target leaf j is above the subtree being moved

        int nodeNumber;
        double newParentHeight;
        boolean attachingToLeaf;
        boolean adjacentEdge;
        //boolean adjacentLeaf;
        do {
            adjacentEdge = false;
            nodeNumber = Randomizer.nextInt(nodeCount + leafNodeCount);
            if (nodeNumber < nodeCount) { // Attaching to branch
                j = tree.getNode(nodeNumber);
                jP = j.getParent();
                if (jP != null)
                    newParentHeight = jP.getHeight();
                else // Attaching to root branch
                    newParentHeight = Double.POSITIVE_INFINITY;
                if (!CiP.isDirectAncestor())
                    adjacentEdge = (CiP.getNr() == j.getNr() || iP.getNr() == j.getNr());
                attachingToLeaf = false;
            } else { // Attaching to leaf
                j = tree.getExternalNodes().get(nodeNumber - nodeCount);
                jP = j.getParent();
                newParentHeight = j.getHeight();
                attachingToLeaf = true;
            }
        } while (j.isDirectAncestor() || // Must not attach to zero length branch
                (newParentHeight <= i.getHeight()) || // New parent must be above selected branch
                (i.getNr() == j.getNr()) ||
                adjacentEdge );  // Adjacent edge must not be selected

        if (attachingToLeaf && iP.getNr() == j.getNr()) {
            System.out.println("Proposal failed because j = iP");
            return Double.NEGATIVE_INFINITY;
        }

        if (jP != null && jP.getNr() == i.getNr()) {
            System.out.println("Proposal failed because jP = i. Heights of i = " + i.getHeight() + " Height of jP = " + jP.getHeight());
            return Double.NEGATIVE_INFINITY;
        }

        //oldDimension = nodeCount - tree.getDirectAncestorNodeCount() - 1;
        oldDimension = allowableNodeCount;
        StratigraphicRange pruningRange = null;
        StratigraphicRange attachingRange = null;

        //classify the type of move being performed before changing the tree structure
        boolean pruningFromSA = CiP.isDirectAncestor();
        boolean pruningFromSRange = !CiP.isDirectAncestor() && sRangeInternalNodeNrs.contains(iP.getNr());
        if (pruningFromSRange || pruningFromSA) {
            pruningRange = tree.getRangeOfNode(iP);
        }
        boolean attachingToSRange = !attachingToLeaf && jP != null && tree.belongToSameSRange(jP.getNr(),j.getNr());
        if (attachingToLeaf || attachingToSRange) {
            attachingRange = tree.getRangeOfNode(j);
        }



        //Hastings numerator calculation + newAge of iP
        if (attachingToLeaf) {
            newRange = 1;
            newAge = j.getHeight();
        } else {
            if (jP != null) {
                newMinAge = Math.max(i.getHeight(), j.getHeight());
                newRange = jP.getHeight() - newMinAge;
                newAge = newMinAge + (Randomizer.nextDouble() * newRange);
            } else {
                double randomNumberFromExponential;
                randomNumberFromExponential = Randomizer.nextExponential(1);
                //newRange = x0 - j.getHeight();
                //randomNumberFromExponential = Randomizer.nextDouble() * newRange;
                newRange = Math.exp(randomNumberFromExponential);
                newAge = j.getHeight() + randomNumberFromExponential;
            }
        }

        Node PiP = iP.getParent();



        //Hastings denominator calculation
        if (CiP.isDirectAncestor()) {
            oldRange = 1;
        }
        else {
            oldMinAge = Math.max(i.getHeight(), CiP.getHeight());
            if (PiP != null) {
                oldRange = PiP.getHeight() - oldMinAge;
            } else {
                oldRange = Math.exp(iP.getHeight() - oldMinAge);
                //oldRange = x0 - oldMinAge;
            }
        }

        boolean oldPossiblySymmetric = !pruningFromSA && !pruningFromSRange;
        boolean oldSymmetric = oldPossiblySymmetric && tree.getNodeIsSymmetric(iP);

        if(oldPossiblySymmetric){ // Remove symmetry of old event before moving
            // Is this necessary?
            tree.setNodeIsSymmetric(iP, false);
        }

        boolean newPossiblySymmetric = !attachingToSRange && !attachingToLeaf;
        boolean newSymmetric = newPossiblySymmetric && Randomizer.nextBoolean();

        //update
        if (iP.getNr() != j.getNr() && CiP.getNr() != j.getNr()) { // Not special case

            // Remove from old location

            iP.removeChild(CiP); //remove <iP, CiP>

            if (PiP != null) { // iP is not root
                boolean left = PiP.getLeft().getNr() == iP.getNr();
                Node anotherChild;
                if (left) {
                    anotherChild = PiP.getRight();
                } else {
                    anotherChild = PiP.getLeft();
                }

                PiP.removeChild(iP);   // remove <PiP,iP>
                CiP.setParent(PiP);
                if (left) {
                    PiP.setLeft(CiP);
                    PiP.setRight(anotherChild);
                } else {
                    PiP.setLeft(anotherChild);
                    PiP.setRight(CiP);
                }  // add <PiP, CiP> keeping the orientation of the edge <PiP,iP>, that is, Or(CiP)=Or(iP)
                PiP.makeDirty(Tree.IS_FILTHY);
                CiP.makeDirty(Tree.IS_FILTHY);
            } else { // iP is root
                CiP.setParent(null); // completely remove <iP, CiP>
                tree.setRootOnly(CiP);
            }

            // Attach to new location
            boolean jLeft= (jP != null) && jP.getLeft().getNr() == j.getNr();

            if (jP != null) {
                Node CjP;
                if (jLeft) {
                    CjP = jP.getRight();
                } else {
                    CjP = jP.getLeft();
                }
                jP.removeChild(j);  // remove <jP, j>
                iP.setParent(jP);
                if (jLeft) {
                    jP.setLeft(iP);
                    jP.setRight(CjP);
                } else {
                    jP.setLeft(CjP);
                    jP.setRight(iP);
                } // add <jP, iP> keeping the orientation of the edge <jP, j>

                jP.makeDirty(Tree.IS_FILTHY);
            } else {
                iP.setParent(null); // completely remove <PiP, iP>
                tree.setRootOnly(iP);
            }
            j.setParent(iP);

            if(newPossiblySymmetric){
                if(newSymmetric){
                    iP.setLeft(i);
                    iP.setRight(j);
                    tree.setNodeIsSymmetric(iP, true);
                    iP.sort();
                } else { // New event is asymmetric
                    if(Randomizer.nextBoolean()){ // Randomly choose the orientation
                        iP.setLeft(j);
                        iP.setRight(i);
                    } else {
                        iP.setLeft(i);
                        iP.setRight(j);
                    }
                }
            } else {
                if (attachingToSRange || (attachingToLeaf && !jLeft)) {
                    iP.setLeft(j);
                    iP.setRight(i);
                } else {
                    iP.setLeft(i);
                    iP.setRight(j);
                }
            }

            iP.makeDirty(Tree.IS_FILTHY);
            j.makeDirty(Tree.IS_FILTHY);
        } else {
            // special case 1: iP = j when pruning from sampled ancestor iP and attaching to the branch above (PiP, iP)
            // special case 2: CiP = j when pruning from a branch (PiP, CiP) and attaching to a leaf CiP
            // In both cases the internal tree structure does not change, only the height of iP
            if (iP.getNr() == j.getNr()) {
                if (attachingToSRange) {
                    //in special case 1: when attaching to the range
                    //make i right
                    //otherwise choose randomly
                    iP.setLeft(CiP);
                    iP.setRight(i);
                } else {
                    if(newSymmetric){
                        iP.setLeft(i);
                        iP.setRight(j);
                        tree.setNodeIsSymmetric(iP, true);
                        iP.sort();
                    } else { // New event is asymmetric
                        if(Randomizer.nextBoolean()){ // Randomly choose the orientation
                            iP.setLeft(j);
                            iP.setRight(i);
                        } else {
                            iP.setLeft(i);
                            iP.setRight(j);
                        }
                    }
                    iP.setLeft(i);
                    iP.setRight(CiP);
                }
            }
            if (CiP.getNr() == j.getNr()) {
                if (PiP == null || PiP.getLeft().getNr() == iP.getNr()) {
                    //in special case 2: when iP is not the root
                    //make i the same orientation as iP
                    //otherwise make i left
                    // Is attachingToLeaf always true???
                    iP.setLeft(i);
                    iP.setRight(CiP);
                } else {
                    iP.setLeft(CiP);
                    iP.setRight(i);
                }
            }
        }
        iP.setHeight(newAge);

        // remove or add nodes to the ranges
        if (pruningFromSA) {
            pruningRange.removeNodeNr(iP.getNr());
            pruningRange.addNodeNr(CiP.getNr());
        }
        if (attachingToLeaf) {
            attachingRange.removeNodeNr(j.getNr());
            attachingRange.addNodeNr(iP.getNr());
        }

        //newDimension = nodeCount - tree.getDirectAncestorNodeCount() - 1;

        newDimension = 0;
        sRangeInternalNodeNrs = tree.getSRangesInternalNodeNrs();
        for (int index=0; index<nodeCount; index++) {
            Node node = tree.getNode(index);
            //the node is not the root, it is not a sampled ancestor on a zero branch, it is not an internal node of a stratigraphic range
            if (!node.isRoot() && !node.isDirectAncestor() && !sRangeInternalNodeNrs.contains(node.getNr()))
                newDimension++;
        }
        dimensionCoefficient = (double) oldDimension / newDimension;

        double symmetryCoefficient = (oldPossiblySymmetric ? 0.5 : 1.0)/(newPossiblySymmetric ? 0.5 : 1.0);
        double orientationCoefficient = (oldSymmetric ? -2.0 : 1.0)/(newSymmetric ? 2.0 : 1.0);

//        for (StratigraphicRange range:sRangeSet.getRanges()) {
//            ArrayList<Node> nodes = (ArrayList) range.getNodes();
//            for (int index=0; index< nodes.size(); index++) {
//                Node node = nodes.get(index);
//                int nodeNr = node.getNr();
//                if (!tree.getNode(nodeNr).equals(node)) {
//                    System.out.println("Node "+ node.toString() + " is not equal to " + tree.getNode(nodeNr));
//                    System.out.println("in range " + range.getFirstOccurrenceID() + ". Resulting tree: ");
//                    System.out.println(tree.getRoot().toString());
//                    System.out.println("node i is " + i.getNr() + ", node j is " + j.getNr());
//                }
//            }
//        }



        fHastingsRatio = Math.abs(symmetryCoefficient * orientationCoefficient * dimensionCoefficient * newRange / oldRange);

        return Math.log(fHastingsRatio);

    }
}
