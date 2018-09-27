package speciation;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.SRMixedTree;
import beast.evolution.tree.SRTree;
import beast.evolution.tree.TreeInterface;
import sranges.StratigraphicRange;

public class SRangesMixedBirthDeathModel extends SRangesBirthDeathModel {
    public Input<RealParameter> symmetricSpeciationProbability = new Input<RealParameter>("symmetricSpeciationProbability",
            "Probability that a speciation event is symmetric", Input.Validate.REQUIRED);

    public Input<RealParameter> anageneticSpeciationRate = new Input<RealParameter>("anageneticSpeciationRate",
            "Rate that anagenetic speciation occurs", Input.Validate.REQUIRED);

    private double lambda_a;
    private double beta;


    @Override
    public void updateParameters(){
        super.updateParameters();
        lambda_a = anageneticSpeciationRate.get().getValue();
        beta = symmetricSpeciationProbability.get().getValue();
    }

    @Override
    protected double q_tilde(double t, double c1, double c2) { // Reciprocal from paper
        return Math.exp((lambda_a + beta * (lambda + mu + psi))*t)*Math.pow(super.q_tilde(t, c1, c2), 1 - beta);
    }

    @Override
    protected double log_q_tilde(double t, double c1, double c2) {
        return (lambda_a + beta * (lambda + mu + psi))*t + (1 - beta)*super.log_q_tilde(t, c1, c2);
    }


    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree)
    {
        final SRMixedTree mixedTree = (SRMixedTree) tree;
        updateParameters();

        //region Constant factors
        int nodeCount = tree.getNodeCount();
        updateParameters();
        if (lambdaExceedsMu && lambda <= mu) {
            return Double.NEGATIVE_INFINITY;
        }

        if (lambda < 0 || mu < 0 || psi < 0) {
            return Double.NEGATIVE_INFINITY;
        }

        double x0 = origin;

        double x1 = tree.getRoot().getHeight();

        if (x0 < x1 ) {
            return Double.NEGATIVE_INFINITY;
        }

        double logPost;
        if (!conditionOnRootInput.get()){
            logPost = -log_q(x0, c1, c2);
        } else {
            if (tree.getRoot().isFake()){   //when conditioning on the root we assume the process
                //starts at the time of the first branching event and
                //that means that the root can not be a sampled ancestor
                return Double.NEGATIVE_INFINITY;
            } else {
                logPost = -log_q(x1, c1, c2);
            }
        }

        if (conditionOnSamplingInput.get()) {
            logPost -= log_oneMinusP0(x0, c1, c2);
        }

        if (conditionOnRhoSamplingInput.get()) {
            if (conditionOnRootInput.get()) {
                logPost -= Math.log(lambda) + log_oneMinusP0Hat(x1, c1, c2)+ log_oneMinusP0Hat(x1, c1, c2);
            }  else {
                logPost -= log_oneMinusP0Hat(x0, c1, c2);
            }
        }
        //endregion

        for (int i = 0; i < nodeCount; i++) {
            Node node = tree.getNode(i);
            if (node.isLeaf()) { //
                if  (!node.isDirectAncestor())  {
                    if (node.getHeight() > 0.000000000005 || rho == 0.) { // Fossil leaf
                        Node fossilParent = node.getParent();
                        if (mixedTree.belongToSameSRange(i, fossilParent.getNr())) { // Stratigraphic range branch above
                            logPost += Math.log(psi) + log_q_tilde(node.getHeight(), c1, c2) + log_p0s(node.getHeight(), c1, c2);
                        } else { // Non-stratigraphic range branch above
                            logPost += Math.log(psi) + log_q(node.getHeight(), c1, c2) + log_p0s(node.getHeight(), c1, c2);
                        }
                    } else {
                        logPost += Math.log(4*rho);
                    }
                }
            } else {
                if (node.isFake()) { // Sampled ancestor point on branch, k
                    logPost += Math.log(psi);
                    Node parent = node.getParent();
                    Node child = node.getNonDirectAncestorChild();
                    if (parent != null && mixedTree.belongToSameSRange(parent.getNr(),i)) { // Stratigraphic range branch above
                        logPost += log_q_tilde(node.getHeight(), c1, c2) - log_q(node.getHeight(), c1, c2);
                    }
                    if (child != null && mixedTree.belongToSameSRange(i,child.getNr())) { // Stratigraphic range branch below
                        logPost += log_q(node.getHeight(), c1, c2)-  log_q_tilde(node.getHeight(), c1, c2);
                    }
                } else { // Speciation event
                    logPost += Math.log(lambda) - log_q(node.getHeight(), c1, c2);
                    logPost += mixedTree.getNodeIsSymmetric(node) ? Math.log(2) + Math.log(beta) : Math.log(1 - beta); // Marginalise orientation
                }
            }
        }
        for (StratigraphicRange range: mixedTree.getSRanges()) {
            Node first =  tree.getNode(range.getNodeNrs().get(0));
            if (!range.isSingleFossilRange()) {
                double tFirst =first.getHeight();
                double tLast = tree.getNode(range.getNodeNrs().get(range.getNodeNrs().size()-1)).getHeight();
                logPost += psi*(tFirst - tLast); //Length of range
            }
            Node ancestralLast = findAncestralRangeLastNode(first);
            if (ancestralLast != null) { // i in I - has ancestral stratigraphic range
                double tOld = ancestralLast.getHeight();
                double tYoung = first.getHeight();
                logPost += log_lambda_times_int_limits_p(tOld, tYoung, c1, c2);
            }
        }

        return logPost;
    }
}
