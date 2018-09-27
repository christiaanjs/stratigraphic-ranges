package speciation;

import beast.evolution.tree.Node;
import beast.evolution.tree.SRMixedTree;
import beast.evolution.tree.SRTree;
import beast.evolution.tree.Tree;
import sranges.StratigraphicRange;

public class SRangesGradientBirthDeathModel extends AbstractSRangesGradientBirthDeathModel {

    @Override
    public double[] calculateTreeLogLikelihoodParamGradient(Tree tree) {
        int nodeCount = tree.getNodeCount();
        updateParameters();
        if (lambdaExceedsMu && lambda <= mu) {
            return new double[PARAM_DIM];
        }

        if (lambda < 0 || mu < 0 || psi < 0) {
            return new double[PARAM_DIM];
        }

        //double x0 = tree.getRoot().getHeight() + origToRootDistance;
        double x0 = origin;

//        if (taxonInput.get() != null) { //TODO rewrite this part for SA sRanges
//
//            if (taxonAge > origin) {
//                return Double.NEGATIVE_INFINITY;
//            }
//            double logPost = 0.0;
//
//            if (conditionOnSamplingInput.get()) {
//                logPost -= Math.log(oneMinusP0(x0, c1, c2));
//            }
//
//            if (conditionOnRhoSamplingInput.get()) {
//                logPost -= Math.log(oneMinusP0Hat(x0, c1, c2));
//            }
//
//            if (SATaxonInput.get().getValue() == 0) {
//                logPost += Math.log(1 - oneMinusP0(taxonAge, c1, c2));
//            } else {
//                logPost += Math.log(oneMinusP0(taxonAge, c1, c2));
//            }
//
//            return logPost;
//        }

        double x1=tree.getRoot().getHeight();

        if (x0 < x1 ) {
            return new double[PARAM_DIM];
        }

        double[] grad = new double[PARAM_DIM];
        if (!conditionOnRootInput.get()){
            increment(grad, -1.0, log_q_grad(x0));
            grad[ORIGIN] -= log_q_grad_t(x0);
        } else {
            if (tree.getRoot().isFake()){   //when conditioning on the root we assume the process
                //starts at the time of the first branching event and
                //that means that the root can not be a sampled ancestor
                return new double[PARAM_DIM];
            } else {
                increment(grad, -1.0, log_q_grad(x1));
            }
        }

        if (conditionOnSamplingInput.get()) {
            // logPost -= log_oneMinusP0(x0, c1, c2);
            increment(grad, -1.0, log_oneMinusP0_grad(x0));
            grad[ORIGIN] -= log_oneMinusP0_grad_t(x0);
        }

        if (conditionOnRhoSamplingInput.get()) {
            if (conditionOnRootInput.get()) {
                //logPost -= Math.log(lambda) + log_oneMinusP0Hat(x1, c1, c2)+ log_oneMinusP0Hat(x1, c1, c2);
                grad[LAMBDA] -= 1.0/lambda;
                increment(grad, -2.0, log_oneMinusP0Hat_grad(x1, c1, c2)); //TODO: Should this really be twice?
            }  else {
                //logPost -= log_oneMinusP0Hat(x0, c1, c2);
                increment(grad, -1.0, log_oneMinusP0Hat_grad(x0, c1, c2)); //TODO: Should this really be twice?
                grad[ORIGIN] -= log_oneMinusP0Hat_grad_t(x0);
            }
        }

        // TODO: Include root branch?
        //increment(grad, log_q_grad(x0));
        //grad[ORIGIN] += log_q_grad_t(x0);

        for (int i = 0; i < nodeCount; i++) {
            Node node = tree.getNode(i);
            if (node.isLeaf()) {
                if  (!node.isDirectAncestor())  {
                    if (node.getHeight() > 0.000000000005 || rho == 0.) {
                        Node fossilParent = node.getParent();
                        if (fossilParent != null && ((SRTree)tree).belongToSameSRange(i, fossilParent.getNr())) {
                            // logPost += Math.log(psi) + log_q_tilde(node.getHeight(), c1, c2) + log_p0s(node.getHeight(), c1, c2);
                            grad[PSI] += 1.0/psi;
                            increment(grad, log_q_tilde_grad(node.getHeight()));
                            increment(grad, log_p0s_grad(node.getHeight()));
                        } else {
                            //logPost += Math.log(psi) + log_q(node.getHeight(), c1, c2) + log_p0s(node.getHeight(), c1, c2);
                            grad[PSI] += 1.0/psi;
                            increment(grad, log_q_grad(node.getHeight()));
                            increment(grad, log_p0s_grad(node.getHeight()));
                        }
                    } else {
                        //logPost += Math.log(4*rho);
                        grad[RHO] += 1.0/rho;
                    }
                }

            } else {
                if (node.isFake()) {
                    //logPost += Math.log(psi);
                    grad[PSI] += 1.0/psi;
                    Node parent = node.getParent();
                    Node child = node.getNonDirectAncestorChild();
                    if (parent != null && ((SRTree)tree).belongToSameSRange(parent.getNr(),i)) {
                        //logPost += log_q_tilde(node.getHeight(), c1, c2) - log_q(node.getHeight(), c1, c2);
                        increment(grad, log_q_tilde_grad(node.getHeight()));
                        increment(grad, -1.0, log_q_grad(node.getHeight()));
                    }
                    if (child != null && ((SRTree)tree).belongToSameSRange(i,child.getNr())) {
                        //logPost += log_q(node.getHeight(), c1, c2)-  log_q_tilde(node.getHeight(), c1, c2);
                        increment(grad, log_q_grad(node.getHeight()));
                        increment(grad, -1.0, log_q_tilde_grad(node.getHeight()));
                    }
                } else {
                    //logPost += Math.log(lambda) - log_q(node.getHeight(), c1, c2);
                    grad[LAMBDA] += 1.0/lambda;
                    increment(grad, -1.0, log_q_grad(node.getHeight()));
                }
            }
        }

        double[] srangeContributions = new double[((SRTree)tree).getSRanges().size()];
        for (StratigraphicRange range:((SRTree)tree).getSRanges()) {
            Node first =  tree.getNode(range.getNodeNrs().get(0));
            if (!range.isSingleFossilRange()) {
                double tFirst =first.getHeight();
                double tLast = tree.getNode(range.getNodeNrs().get(range.getNodeNrs().size()-1)).getHeight();
                //logPost += psi*(tFirst - tLast);
                grad[PSI] += tFirst - tLast;
            }
            Node ancestralLast = findAncestralRangeLastNode(first);
            if (ancestralLast != null) {
                double tOld = ancestralLast.getHeight();
                double tYoung = first.getHeight();
                //logPost += log_lambda_times_int_limits_p(tOld, tYoung, c1, c2);
                increment(grad, log_lambda_times_int_limits_p_grad(tOld, tYoung));
            }
        }


        return grad;
    }

    private double log_oneMinusP0_grad_t(double t) {
        //return
        // Math.log(1 - (lambda + mu + psi - c1 * ((1 + c2) - Math.exp(-c1 * t) * (1 - c2)) / ((1 + c2) + Math.exp(-c1 * t) * (1 - c2))) / (2 * lambda));
        double oneMinusP0 = oneMinusP0(t, c1, c2);
        double denomTerm = 1 + c2 + Math.exp(-c1 * t)*(1 - c2);
        return (c1 * c1 / (lambda * oneMinusP0)) * (1 + c2) * (1 - c2) * Math.exp(-c1 * t)*(1 - c2) / (denomTerm * denomTerm);
    }

    private double log_q_grad_t(double t) {
        // return Math.log(Math.exp(c1 * t) * (1 + c2) * (1 + c2) + Math.exp(-c1 * t) * (1 - c2) * (1 - c2) + 2 * (1 - c2 * c2));
        double q_grad = c1 * Math.exp(c1 * t) * (1 + c2) * (1 + c2) - c1 * Math.exp(-c1 * t) * (1 - c2) * (1 - c2);
        return q_grad/q(t, c1, c2);
    }

    private double log_oneMinusP0Hat_grad_t(double t) {
        throw new RuntimeException("Not implemented");
    }

    private double[] log_q_tilde_grad(double t) {
        //0.5*(t*(lambda + mu + psi) + log_q(t,c1,c2));
        double[] grad = log_q_grad(t);
        grad[LAMBDA] += t;
        grad[MU] += t;
        grad[PSI] += t;
        scale(grad, 0.5);
        return grad;
    }

    private double[] log_lambda_times_int_limits_p_grad(double tOld, double tYoung) {
        double firstTerm = int_limits_term(tYoung);
        double secondTerm = int_limits_term(tOld);
        double lambda_times_int_limits_p = (lambda+mu+psi-c1)*(tOld - tYoung) + 2*Math.log(firstTerm)
           + 2*Math.log(secondTerm);

        double[] grad = new double[PARAM_DIM];
        System.arraycopy(c2_grad, 0, grad, 0, PARAM_DIM);
        grad[LAMBDA] -= 1;
        grad[MU] -= 1;
        grad[PSI] -= 1;
        scale(grad, tYoung - tOld);

        double[] firstTermGrad = int_limits_term_grad(tYoung);
        double[] secondTermGrad = int_limits_term_grad(tOld);

        for(int i = 0; i < PARAM_DIM; i++){
            grad[i] += 2.0/firstTerm*firstTermGrad[i];
            grad[i] -= 2.0/secondTerm*secondTermGrad[i];
        }

        scale(grad, 1.0/lambda_times_int_limits_p);
        return grad;
    }

    private double int_limits_term(double t){
        return Math.exp(-c1*t)*(1-c2)+(1+c2);
    }

    private double[] int_limits_term_grad(double t){
        double[] grad = new double[PARAM_DIM];
        double expMinusC1T = Math.exp(-c1 * t);
        for(int i = 0; i < PARAM_DIM; i++){
            grad[i] = -t * expMinusC1T * (1 - c2) * c1_grad[i] + (1 - expMinusC1T)*c2_grad[i];
        }
        return grad;
    }

    private double[] log_oneMinusP0Hat_grad(double x1, double c1, double c2) {
        throw new RuntimeException("Not implemented");
    }

    private double[] log_p0s_grad(double t) {
        double p0 = (lambda + mu + psi - c1 * ((1 + c2) - Math.exp(-c1 * t) * (1 - c2)) / ((1 + c2) + Math.exp(-c1 * t) * (1 - c2))) / (2 * lambda);
        double p0s = r + (1 - r) * p0;

        // c1 * ((1 + c2) - Math.exp(-c1 * t) * (1 - c2)) / ((1 + c2) + Math.exp(-c1 * t) * (1 - c2)))
        double[] grad = new double[PARAM_DIM];

        double expMinusC1T = Math.exp(-c1 * t);
        double f = (1 + c2);
        double g = expMinusC1T*(1 - c2);
        double fPlusG = f + g;
        double ratio = (f - g)*(f + g);
        double ratioGradDenom = fPlusG * fPlusG;
        for(int i = 0; i < PARAM_DIM; i++){
            double fGrad = c2_grad[i];
            double gGrad = -expMinusC1T*t*c1_grad[i]*(1 - c2)  - expMinusC1T*c2_grad[i];
            double ratioGrad = 2.0*(g*fGrad - 2*f*gGrad)/ratioGradDenom;
            grad[i] = ratio*c2_grad[i] + ratioGrad*c2;
        }
        double num = lambda + mu + psi - c1 * ((1 + c2) - Math.exp(-c1 * t) * (1 - c2)) / ((1 + c2) + Math.exp(-c1 * t) * (1 - c2));
        grad[LAMBDA] += 1;
        grad[LAMBDA] = grad[LAMBDA] - 2*num/lambda;
        grad[MU] += 1;
        grad[PSI] += 1;
        scale(grad, 1.0/(2.0 * lambda));

        grad[REMOVAL_PROBABILITY] = (1.0 - p0)/(1.0 - r);
        scale(grad,  (1.0 - r)/p0s);

        return grad;
    }

    private double[] log_oneMinusP0_grad(double t) {
        //return Math.log(1 - (lambda + mu + psi - c1 * ((1 + c2) - Math.exp(-c1 * t) * (1 - c2)) / ((1 + c2) + Math.exp(-c1 * t) * (1 - c2))) / (2 * lambda));
        throw new RuntimeException("Not implemented");
    }

    private double[] log_q_grad(double t) {
        //return Math.log(Math.exp(c1 * t) * (1 + c2) * (1 + c2) + Math.exp(-c1 * t) * (1 - c2) * (1 - c2) + 2 * (1 - c2 * c2));
        double oneOverQ = 1.0/q(t, c1, c2);
        double[] grad = new double[PARAM_DIM];

        // First term
        increment(grad, Math.exp(c1 * t) * t * (1 + c2) * (1 + c2), c1_grad);
        increment(grad, Math.exp(c1 * t) * 2.0 * (1 + c2), c2_grad);

        // Second term
        increment(grad, -Math.exp(-c1 * t) * t * (1 - c2) * (1 - c2), c1_grad);
        increment(grad, -Math.exp(-c1 * t) * 2.0 * (1 - c2), c2_grad);

        // Third term
        increment(grad, -4.0 * (1 - c2), c2_grad);

        scale(grad, oneOverQ);

        return grad;
    }

    private double[] c1_grad;
    private double[] c2_grad;

    @Override
    protected void updateParameters(){
        super.updateParameters();
        // c1 = Math.sqrt((lambda - mu - psi) * (lambda - mu - psi) + 4 * lambda * psi);
        // c2 = -(lambda - mu - 2*lambda*rho - psi) / c1;
        double netRate = lambda - mu - psi;
        c1_grad = new double[PARAM_DIM];
        c1_grad[LAMBDA] = (2.0*netRate + 4.0*psi)/(2.0 * c1);
        c1_grad[MU] = (-2.0*netRate)/(2.0 * c1);
        c1_grad[PSI] = (-2.0*netRate + 4.0*lambda)/(2.0 * c1);

       double c2_num = -(lambda - mu - 2*lambda*rho - psi);
       double c1_sq = c1 * c1;
       c2_grad = new double[PARAM_DIM];
       c2_grad[LAMBDA] = ((-1.0 - 2.0 * rho) * c1 - c2_num * c1_grad[LAMBDA])/c1_sq;
       c2_grad[MU] = (-c1 - c2_num * c1_grad[MU])/c1_sq;
       c2_grad[PSI] = (-c1 - c2_num*c1_grad[PSI])/c1_sq;
       c2_grad[RHO] = (-2.0 * lambda)/c1_sq;
    }

    private static final int ORIGIN = 0;
    private static final int LAMBDA = 1;
    private static final int MU = 2;
    private static final int PSI = 3;
    private static final int RHO = 4;
    private static final int REMOVAL_PROBABILITY = 5;

    private static void increment(double[] target, double scale, double[] x) {
        for (int i = 0; i < target.length; i++) {
            target[i] += scale * x[i];
        }
    }

    private static void increment(double[] target, double[] x) {
        for (int i = 0; i < target.length; i++) {
            target[i] += x[i];
        }
    }

    private static void scale(double[] x, double scale){
        for(int i = 0; i < x.length; i++){
            x[i] *= scale;
        }
    }

}
