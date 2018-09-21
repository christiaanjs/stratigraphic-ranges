package sranges;

import beast.core.Loggable;
import beast.evolution.tree.SRMixedTree;
import beast.util.Randomizer;
import javafx.util.Pair;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class SRSimulationApp {

    public static void main(String[] args) throws IOException {
        double x0 = 3.0;
        double lambda = 1.0;
        double mu = 0.2;
        double psi = 1.0;
        double rho = 0.8;
        double beta = 0.2;
        double lambda_a = 0.3;
        double[] simArgs = new double[]{x0, lambda, mu, psi, rho, beta, lambda_a};

        Randomizer.setSeed(123);

        int nSims = 10;
        int printSimsEvery = 10000;
        int printAcceptedEvery = 1;
        int numberOfRanges = 3;

        List<Pair<Boolean, Integer>> targetRanges = new ArrayList<>(numberOfRanges); // Has rho samples, number of samples
        targetRanges.add(new Pair<>(false, 1));
        targetRanges.add(new Pair<>(false, 2));
        targetRanges.add(new Pair<>(true, 2));

        targetRanges.sort(SRangesUtil.CANONICAL_ORDER);

        Predicate<SRMixedTree> condition = t -> (t.getNumberOfRanges() == numberOfRanges) && (
                SRangesUtil.listEqual(
                        t.getSRanges().stream()
                                .map(SRangesUtil.getPredicateMapper(t))
                                .sorted(SRangesUtil.CANONICAL_ORDER)
                                .collect(Collectors.toList())
                        , targetRanges));

        SRMixedTree initTree = Stream.generate(() -> SRMixedTreeSimulator.doSimulation(simArgs))
                .filter(condition)
                .iterator().next();

        sortAndAddRangeIds(initTree);

        SRMixedStatsLogger logger = new SRMixedStatsLogger(initTree);
        String outFile = "conditioned-sranges.log";
        PrintStream out = new PrintStream(outFile);

        logger.init(out);
        out.println();

        final AtomicInteger simCount = new AtomicInteger(0);
        final AtomicInteger acceptedCount = new AtomicInteger(0);

        Stream.generate(() -> SRMixedTreeSimulator.doSimulation(simArgs))
                .parallel()
                .peek(t -> {
                    int i = simCount.getAndIncrement();
                    if (i % printSimsEvery == 0) System.out.printf("Performed %d simulations\n", i);
                })
                .filter(condition)
                .limit(nSims)
                .peek(SRSimulationApp::sortAndAddRangeIds)
                .map(t -> new Pair<>(acceptedCount.getAndIncrement(), t))
                .peek(p -> {
                    int i = p.getKey();
                    if (i % printAcceptedEvery == 0) System.out.printf("%d/%d accepted\n", i, nSims);
                })
                .forEach(p -> {
                    logger.setTree(p.getValue());
                    logger.log(p.getKey(), out);
                    out.println();
                });

        out.close();

    }

    private static void sortAndAddRangeIds(SRMixedTree tree){
        List<StratigraphicRange> sranges = tree.getSRanges();
        sranges.sort(SRangesUtil.getCanonicalOrder(tree));
        for(int i = 0; i < sranges.size(); i++){
            sranges.get(i).setID("range" + i);
        }
    }



}
