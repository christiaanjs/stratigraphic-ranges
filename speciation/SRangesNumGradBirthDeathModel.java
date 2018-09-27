package speciation;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;

import java.util.Arrays;
import java.util.stream.Collectors;

public class SRangesNumGradBirthDeathModel extends AbstractSRangesGradientBirthDeathModel {

    @Override
    public double[] calculateTreeLogLikelihoodParamGradient(Tree tree) {
        String[] paramNames = getParamInputNames();
        double[] x = getParams();
        double[] h = Arrays.stream(x).map(y -> Math.sqrt(Math.ulp(y)) * y).toArray();

        double[] grad = new double[PARAM_DIM];
        for (int i = 0; i < PARAM_DIM; i++) {
            setInputValue(paramNames[i], Double.toString(x[i] - h[i]));
            double a = calculateTreeLogLikelihood(tree);
            setInputValue(paramNames[i], Double.toString(x[i] + h[i]));
            double b = calculateTreeLogLikelihood(tree);
            grad[i] = (b - a) / (2.0 * h[i]);

            setInputValue(paramNames[i], Double.toString(x[i]));

        }
        return grad;
    }

    private double[] getParams(){
        return Arrays.stream(getParamInputNames())
                .mapToDouble(p -> ((RealParameter) getInputValue(p)).getValue() )
                .toArray();
    }
}
