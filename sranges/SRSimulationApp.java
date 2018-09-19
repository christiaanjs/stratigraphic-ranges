package sranges;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.SRMixedTree;
import beast.util.Randomizer;
import com.sun.org.apache.xpath.internal.operations.Bool;
import javafx.util.Pair;
import speciation.SRangesMixedBirthDeathModel;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static java.util.Comparator.comparing;

public class SRSimulationApp {
    public static void main(String[] args) throws IOException {
        double x0 = 3.0;
        double lambda = 1.0;
        double mu = 0.2;
        double psi = 1.0;
        double rho = 0.8;
        double beta = 0.2;
        double lambda_a = 0.3;
        double[] simArgs = new double[]{ x0, lambda, mu, psi, rho, beta, lambda_a };

        Randomizer.setSeed(123);

        int nSims = 10000;
        int printSimsEvery = 10000;
        int printAcceptedEvery = 100;
        int numberOfRanges = 3;

        Function<SRMixedTree, Function<StratigraphicRange, Pair<Boolean, Integer>>> predicateMapper = t ->
                (r -> new Pair<>(t.rangeHasRhoSample(r), r.getNodeNrs().size()));

        List<Pair<Boolean, Integer>> targetRanges = new ArrayList<>(numberOfRanges); // Has rho samples, number of samples
        targetRanges.add(new Pair<>(false, 1));
        targetRanges.add(new Pair<>(false, 2));
        targetRanges.add(new Pair<>(true, 2));

        Comparator<Pair<Boolean, Integer>> canonicalOrder = comparing((Pair<Boolean, Integer> p) -> p.getKey())
                .thenComparing(p -> p.getValue());
        targetRanges.sort(canonicalOrder);

        Predicate<SRMixedTree> condition = t -> (t.getNumberOfRanges() == numberOfRanges) &&(
            listEqual(
                    t.getSRanges().stream()
                            .map(predicateMapper.apply(t))
                            .sorted(canonicalOrder)
                            .collect(Collectors.toList())
                    , targetRanges));

        String outFile = "conditioned-sranges.trees";
        PrintWriter printWriter = new PrintWriter(new FileWriter(outFile));

        String header = String.join(",", IntStream.range(0, numberOfRanges)
                .mapToObj((int i) -> String.format("hasPsi_%d,nSamples_%d", i, i)).collect(Collectors.toList()));
        printWriter.println(header);

        final AtomicInteger simCount = new AtomicInteger(0);
        final AtomicInteger acceptedCount = new AtomicInteger(0);

        Stream.generate(() -> SRMixedTreeSimulator.doSimulation(simArgs))
                .parallel()
                .peek(t -> {
                    int i = simCount.getAndIncrement();
                    if(i % printSimsEvery == 0) System.out.printf("Performed %d simulations\n", i);
                })
                .filter(condition)
                .limit(nSims)
                .map(t -> new Pair<>(acceptedCount.getAndIncrement(), t))
                .peek(p -> {
                    int i = p.getKey();
                    if(i % printAcceptedEvery == 0) System.out.printf("%d/%d accepted\n", i, nSims);
                })
                .forEach(p -> printWriter.println(String.format(
                        "tree SIM_%d = %s;",
                        p.getKey(),
                        p.getValue().toString()
                )));

        printWriter.close();

    }

    private static <T> boolean listEqual(List<T> a, List<T> b){
        if(a.size() != b.size()) return false;
        for(int i = 0; i < a.size(); i++){
            if(!a.get(i).equals(b.get(i))){
                return false;
            }
        }
        return true;
    }

}
