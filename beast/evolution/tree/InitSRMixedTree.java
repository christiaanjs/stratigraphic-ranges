package beast.evolution.tree;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import sranges.SRangesUtil;
import sranges.StratigraphicRange;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class InitSRMixedTree extends SRMixedTree implements StateNodeInitialiser {

    public Input<RealParameter> originInput = new Input<>("origin", "height of origin of process", Input.Validate.REQUIRED);

    @Override
    public void initAndValidate(){
        super.initAndValidate();
        initStateNodes();
    }

    @Override
    public void initStateNodes() {
        List<String> taxa = m_initial.get().m_taxonset.get().asStringList();
        List<TraitSet> traits = m_initial.get().m_traitList.get();
        Map<String, Node> taxaNodes = new HashMap<String, Node>();

        final double maxHeight = originInput.get().getValue();
        double maxTaxonHeight = Double.NEGATIVE_INFINITY;

        for(String t: taxa){
            Node newNode = newNode();
            newNode.setID(t);
            for(TraitSet ts: traits)
                newNode.setMetaData(ts.getTraitName(), ts.getValue(t));
            maxTaxonHeight = Math.max(maxTaxonHeight, newNode.getHeight());
            taxaNodes.put(t, newNode);
        }

        List<StratigraphicRange> sranges = stratigraphicRangeInput.get();
        nodeCount = taxa.size() * 2 - 1;
        leafNodeCount = taxa.size();
        internalNodeCount = taxa.size() - 1;

        List<Node> candidates = new LinkedList<>();

        for(StratigraphicRange r: sranges){
            if(r.isSingleFossilRange()) {
                candidates.add(taxaNodes.get(r.getFirstOccurrenceID()));
            } else {
                Node first = taxaNodes.get(r.getFirstOccurrenceID());
                Node last = taxaNodes.get(r.getLastOccurrenceID());
                Node fake = newNode();
                fake.setHeight(first.getHeight());
                fake.addChild(first);
                fake.addChild(last);
                candidates.add(fake);
            }
        }

        int remainingNodeCount = candidates.size() - 1;
        final double minHeight = maxTaxonHeight;
        List<Double> nodeTimes = Stream.generate(() -> Randomizer.nextDouble()*(maxHeight - minHeight) + minHeight)
                .limit(remainingNodeCount)
                .sorted()
                .collect(Collectors.toList());

        for(double height: nodeTimes){
            Node child1 = candidates.remove(Randomizer.nextInt(candidates.size()));
            Node child2 = candidates.remove(Randomizer.nextInt(candidates.size()));
            Node newNode = newNode();
            newNode.setHeight(height);
            newNode.addChild(child1);
            newNode.addChild(child2);
            candidates.add(newNode);
        }

        root = candidates.get(0);
        SRangesUtil.numberNodes(root);
        initArrays();
        initSymmetric();


        if (m_initial.get() != null) {
            m_initial.get().assignFromWithoutID(this);
        }

    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) { }
}
