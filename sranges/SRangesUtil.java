package sranges;

import beast.evolution.tree.Node;
import beast.evolution.tree.SRMixedTree;
import beast.evolution.tree.SRTree;
import javafx.util.Pair;

import java.util.Comparator;
import java.util.List;
import java.util.function.Function;

import static java.util.Comparator.comparing;

public class SRangesUtil {
    public static final Comparator<Pair<Boolean, Integer>> CANONICAL_ORDER = comparing((Pair<Boolean, Integer> p) -> p.getKey())
            .thenComparing(p -> p.getValue());

    public static final Function<StratigraphicRange, Pair<Boolean, Integer>> getPredicateMapper(SRTree tree){
        return r -> new Pair<>(tree.rangeHasRhoSample(r), r.getNodeNrs().size());
    }

    public static final Comparator<StratigraphicRange> getCanonicalOrder(SRTree tree){
        Function<StratigraphicRange, Pair<Boolean, Integer>> mapper = getPredicateMapper(tree);
        return (r1, r2) -> CANONICAL_ORDER.compare(mapper.apply(r1), mapper.apply(r2));
    }

    static <T> boolean listEqual(List<T> a, List<T> b){
        if(a.size() != b.size()) return false;
        for(int i = 0; i < a.size(); i++){
            if(!a.get(i).equals(b.get(i))){
                return false;
            }
        }
        return true;
    }

    public static Pair<Integer, Integer> numberNodes(Node root){
        return numberNodes(root, new Pair(0, root.getLeafNodeCount()));
    }

    public static Pair<Integer, Integer> numberNodes(Node node, Pair<Integer, Integer> counter){ // Leaf count, internal node count
        if(node.isLeaf()){
            node.setNr(counter.getKey());
            return new Pair<>(counter.getKey() + 1, counter.getValue());
        } else {
            for(Node child: node.getChildren())
                counter = numberNodes(child, counter);
            node.setNr(counter.getValue());
            return new Pair<>(counter.getKey(), counter.getValue() + 1);
        }
    }
}
