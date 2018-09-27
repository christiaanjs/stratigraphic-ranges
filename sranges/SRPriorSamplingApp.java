package sranges;

import beast.app.beastapp.BeastMain;
import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import speciation.SRangesMixedBirthDeathModel;

import java.io.IOException;

public class SRPriorSamplingApp {
    public static void main(String[] args) throws IOException {

//        double x0 = 3.0;
//        double lambda = 1.0;
//        double mu = 0.2;
//        double psi = 1.0;
//        double rho = 0.8;
//        double beta = 0.2;
//        double lambda_a = 0.3;
//
//        SRangesMixedBirthDeathModel l = new SRangesMixedBirthDeathModel();
//
//        //l.setInputValue("tree", tree);
//        setRealInput(l, l.originInput, x0);
//        setRealInput(l, l.birthRateInput, lambda);
//        setRealInput(l, l.deathRateInput, mu);
//        setRealInput(l, l.samplingRateInput, psi);
//        setRealInput(l, l.samplingProportionInput, rho);
//        setRealInput(l, l.removalProbability, 0.0);
//        setRealInput(l, l.symmetricSpeciationProbability, beta);
//        setRealInput(l, l.anageneticSpeciationRate, lambda_a);
//
//        l.initAndValidate();
//
//        System.out.println("Calculated likelihood");
        //System.out.println(l.calculateTreeLogLikelihood(tree));
        BeastMain.main(new String[]{ "-seed", "1234", "-overwrite", "stratigraphic-ranges/examples/sample-prior.xml" });
    }


    public static <T> void setRealInput(BEASTObject beastObject, Input<T> input, double value){
        beastObject.setInputValue(input.getName(), new RealParameter(Double.toString(value)));
    }
}
