package sranges;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.SRMixedTree;
import javafx.util.Pair;

import java.io.PrintStream;
import java.util.Comparator;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Stream;

public class SRMixedStatsLogger extends CalculationNode implements Loggable {
    public Input<SRMixedTree> treeInput = new Input<>("tree", "Mixed mode stratigraphic range tree to calculate statistics from",
            Input.Validate.REQUIRED);

    private SRMixedTree tree;
    private List<StratigraphicRange> sranges;

    public SRMixedStatsLogger(){}

    public SRMixedStatsLogger(SRMixedTree tree){
        this.tree = tree;
        initSranges();
    }

    public void setTree(SRMixedTree tree){
        this.tree = tree;
        tree.getLeafNodeCount(); // Ensure leaf node count has been cached
        initSranges();
    }

    private void treeStatsInit(PrintStream out){
        String id = tree.getID();
        if (id == null || id.matches("\\s*"))
            id = "tree";
        out.print(
            id + ".rootHeight\t" +
            id + ".sampledAncestorCount\t"
        );
    }

    private void srangeInit(PrintStream out, StratigraphicRange range){
        String id = range.getID();
        out.print(id + ".oldest\t");
        if(!range.isSingleFossilRange() && !tree.rangeHasRhoSample(range)){
            out.print(id + ".youngest\t");
        }
    }

    // TODO: MRCA of pairs of ranges
    // TODO: Speciation events

    @Override
    public void init(PrintStream out) {
        treeStatsInit(out);
        sranges.forEach(r -> srangeInit(out, r));
    }

    private void logTreeStats(PrintStream out){
        out.print(
            tree.getRoot().getHeight() + "\t" +
            tree.getDirectAncestorNodeCount() + "\t"
        );
    }

    private void logRangeStats(PrintStream out, StratigraphicRange range){
        out.print(tree.getNode(range.getFirstOccurrenceNr()).getHeight() + "\t");
        if(!range.isSingleFossilRange() && !tree.rangeHasRhoSample(range)){
            out.print(tree.getNode(range.getLastOccurrenceNr()).getHeight() + "\t");
        }
    }

    @Override
    public void log(long sample, PrintStream out) {
        logTreeStats(out);
        sranges.forEach(r -> logRangeStats(out, r));
    }

    @Override
    public void close(PrintStream out) { }

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        initSranges();
    }

    private void initSranges() {
        sranges = tree.getSRanges();
        sranges.sort(SRangesUtil.getCanonicalOrder(tree));
    }
}
