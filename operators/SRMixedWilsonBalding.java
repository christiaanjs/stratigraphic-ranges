package operators;

import beast.evolution.tree.Node;
import beast.evolution.tree.SRMixedTree;
import beast.evolution.tree.SRTree;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import sranges.StratigraphicRange;

import java.util.*;
import java.util.stream.Collectors;

public class SRMixedWilsonBalding extends SRMixedTreeOperator {


    /**
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        Set<Integer> srangeInternalNodes = new HashSet<>(tree.getSRangesInternalNodeNrs());


        // Select pruned edge
        List<Node> eligibleNodes = getEligibleNodes(srangeInternalNodes);
        int eligibleNodeCount = eligibleNodes.size();
        Node i = eligibleNodes.get(Randomizer.nextInt(eligibleNodeCount));

        Node iP = i.getParent();
        Node CiP = iP.getLeft().getNr() == i.getNr() ? iP.getRight() : iP.getLeft();
        boolean pruningFromSA = iP.isFake();
        boolean pruningFromSRange = srangeInternalNodes.contains(iP.getNr());
        boolean oldPossiblySymmetric = !pruningFromSA && !pruningFromSRange;
        boolean oldSymmetric = oldPossiblySymmetric && tree.getNodeIsSymmetric(iP);

        List<Node> eligibleTargets = getEligibleTargets(i);
        int eligibleTargetCount = eligibleTargets.size();

        if (eligibleTargetCount == 0) {
            return Double.NEGATIVE_INFINITY;
        }

        // Select target edge
        Node j = eligibleTargets.get(Randomizer.nextInt(eligibleTargetCount));
        Node jP = j.getParent();

        boolean attachingToLeaf = j.isLeaf();
        boolean attachingToRange = !attachingToLeaf && !j.isRoot() && tree.belongToSameSRange(jP.getNr(),j.getNr());
        boolean newPossiblySymmetric = !attachingToLeaf && !attachingToRange;
        boolean newSymmetric = newPossiblySymmetric && Randomizer.nextBoolean();

        double oldHeightRange;
        if (iP.isDirectAncestor()) {
            oldHeightRange = 1;
        }
        else {
            double oldMinAge = Math.max(i.getHeight(), CiP.getHeight());
            if (iP.isRoot()) {
                oldHeightRange = Math.exp(iP.getHeight() - oldMinAge);
            } else {
                oldHeightRange = iP.getParent().getHeight() - oldMinAge;
            }
        }

        // Update

        //  // Disconnect pruned
        iP.removeChild(CiP);
        if(iP.isRoot()){
            CiP.setParent(null);
            tree.setRootOnly(CiP);
        } else {
            Node PiP = iP.getParent();
            if(PiP.getLeft().getNr() == iP.getNr()){ //Maintain orientation
                PiP.setLeft(CiP);
            } else {
                PiP.setRight(CiP);
            }
            iP.setParent(null);
        }
        if(pruningFromSRange){
            StratigraphicRange pruningRange = tree.getRangeOfNode(iP);
            pruningRange.removeNodeNr(iP.getNr());
            pruningRange.addNodeNr(CiP.getNr());
        }
        if(oldPossiblySymmetric){
            tree.setNodeIsSymmetric(iP, false);
        }
        CiP.makeDirty(Tree.IS_FILTHY);

        double newHeightRange;

        if(j.isRoot()){
            tree.setRootOnly(iP);
            double delta = Randomizer.nextExponential(1.0);
            iP.setHeight(j.getHeight() + delta);
            newHeightRange = Math.exp(delta);
        } else {
            if(jP.getLeft().getNr() == j.getNr()){ //Maintain orientation
                jP.setLeft(iP);
            } else {
                jP.setRight(iP);
            }
            if(attachingToLeaf){
                iP.setHeight(j.getHeight());
                newHeightRange = 1.0;
            } else {
                double minNewAge = Math.max(i.getHeight(), j.getHeight());
                newHeightRange = jP.getHeight() - minNewAge;
                iP.setHeight(minNewAge + Randomizer.nextDouble()*newHeightRange);
            }
        }

        if(Randomizer.nextBoolean()){
            iP.setLeft(j); // TODO: In possibly symmetric, but asymmetric, orientation is unidentifiable?
            iP.setRight(i);
        } else {
            iP.setLeft(i);
            iP.setRight(j);
        }
        if(newPossiblySymmetric){
            tree.setNodeIsSymmetric(iP, newSymmetric);
        }
        if(attachingToRange){
            StratigraphicRange attachingRange = tree.getRangeOfNode(j);
            attachingRange.addNodeNrAfter(jP.getNr(), iP.getNr());
        }
        if(attachingToLeaf){
            StratigraphicRange attachingRange = tree.getRangeOfNode(j);
            if(attachingRange != null){ // TODO: Why is this sometimes null?
                attachingRange.removeNodeNr(j.getNr());
                attachingRange.addNodeNr(iP.getNr());
            }
        }

        int newEligibleNodeCount = getEligibleNodes(new HashSet<>(tree.getSRangesInternalNodeNrs())).size();

        double dimensionCoefficient = ((double) eligibleNodeCount)/newEligibleNodeCount;
        double symmetryCoefficient = (oldPossiblySymmetric ? 0.5 : 1.0)/(newPossiblySymmetric ? 0.5 : 1.0);
        double orientationCoefficient = (oldSymmetric ? -2.0 : 1.0)/(newSymmetric ? 2.0 : 1.0);

        double fHastingsRatio = Math.abs(symmetryCoefficient * orientationCoefficient * dimensionCoefficient * newHeightRange / oldHeightRange);

        return Math.log(fHastingsRatio);
    }

    private boolean eligibleTargetExists(Node i){
        Node iP = i.getParent();
        Node CiP = iP.getLeft().getNr() == i.getNr() ? iP.getRight() : iP.getLeft();
        return !(iP.isRoot() && CiP.getHeight() <= i.getHeight()); // or only eligible is iP?
    }


    private List<Node> getEligibleNodes(Set<Integer> srangeInternalNodes){
        return Arrays.stream(tree.getNodesAsArray())
                .filter(n -> !n.isDirectAncestor()
                        && !n.isRoot()
                        && !srangeInternalNodes.contains(n.getNr()))
                .collect(Collectors.toList());
    }

    private List<Node> getEligibleTargets(Node i){
        return Arrays.stream(tree.getNodesAsArray())
                .filter(n -> n.getNr() != i.getNr()
                        && (n.isRoot() || n.getParent().getNr() != i.getParent().getNr())
                        && (i.isRoot() || i.getParent().getNr() != n.getNr()) // TODO: Can we remove this condition by adding special case?
                        && !n.isDirectAncestor()
                        && n.getHeight() >= i.getHeight())
                .collect(Collectors.toList());
    }
}
