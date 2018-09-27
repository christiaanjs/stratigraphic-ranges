package sranges;

import beast.evolution.tree.SRMixedTree;
import beast.util.Randomizer;
import speciation.AbstractSRangesGradientBirthDeathModel;
import speciation.SRangesGradientBirthDeathModel;
import speciation.SRangesNumGradBirthDeathModel;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static sranges.SRPriorSamplingApp.setRealInput;

public class SRangesGradientVerification {
    public static void main(String[] args) throws FileNotFoundException {
        double x0 = 3.0;
        double lambda = 1.0;
        double mu = 0.2;
        double psi = 1.0;
        double rho = 0.8;
        double beta = 0.0;
        double lambda_a = 0.0;
        double[] simArgs = new double[]{x0, lambda, mu, psi, rho, beta, lambda_a};

        Randomizer.setSeed(123);

        int nSims = 10000;
        int printEvery = 10;

        Stream<SRMixedTree> simStream = Stream.generate(() -> SRMixedTreeSimulator.doSimulation(simArgs));

        String outFile = "srange-num-gradient.log";
        PrintWriter out = new PrintWriter(outFile);
        out.println(String.join(",", new SRangesGradientBirthDeathModel().getParamInputNames()));

        final AtomicInteger simCount = new AtomicInteger(0);

        Stream.generate(() -> SRMixedTreeSimulator.doSimulation(simArgs))
                .filter(t -> t.getSRanges().size() > 0) // Condition on sampling
                .limit(nSims)
                .parallel()
                .peek(t -> {
                    int i = simCount.getAndIncrement();
                    if (i % printEvery == 0) System.out.printf("Performed %d simulations\n", i);
                })
                .map(t -> getGradient(x0, lambda, mu, psi, rho, t) )
                .map(g -> Arrays.stream(g).mapToObj(x -> Double.toString(x)).collect(Collectors.joining(",")))
                .forEach(out::println);

        out.close();
    }

    private static AbstractSRangesGradientBirthDeathModel newLikelihood(){
        return new SRangesNumGradBirthDeathModel();
    }

    private static double[] getGradient(double x0, double lambda, double mu, double psi, double rho, SRMixedTree tree) {
        AbstractSRangesGradientBirthDeathModel l = newLikelihood();
        setRealInput(l, l.originInput, x0);
        setRealInput(l, l.birthRateInput, lambda);
        setRealInput(l, l.deathRateInput, mu);
        setRealInput(l, l.samplingRateInput, psi);
        setRealInput(l, l.rhoProbability, rho);
        setRealInput(l, l.removalProbability, 0.0);
        l.setInputValue("tree", tree);
        l.initAndValidate();
        return l.calculateTreeLogLikelihoodParamGradient(tree);
    }
}
