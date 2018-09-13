package beast.evolution.tree;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import sranges.StratigraphicRange;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Alexandra Gavryushkina
 *
 * A labeled oriented tree on stratigraphic ranges under budding speciation.
 * At every internal node, the branch leading to the left child represents the ancestral species
 * and the branch leading to the right child represents the descendant species.
 * Branches starting with sampled ancestors inheret the orientation from the parental branches:
 * if a node is fake then the non-direct ancestor child gets the same left/right orientation as the fake node.
 * If the root is a fake node then the non-direct ancestor child is always left.
 */
public class SRTree extends Tree {

    public Input<List<StratigraphicRange>> stratigraphicRangeInput = new Input<>("stratigraphicRange", "all stratigraphic ranges", new ArrayList<>());

    public Input<Tree> treeInput = new Input<>("tree", "tree to start with");

    protected ArrayList<StratigraphicRange> sRanges;
    protected ArrayList<StratigraphicRange> storedSRanges;

    public SRTree(){}

    public SRTree(Node root, List<StratigraphicRange> sRanges){
        super(root);
        this.sRanges = new ArrayList<>(sRanges);
    }

    @Override
    public void initAndValidate() {
        if (treeInput.get() != null) {
            assignFromWithoutID(treeInput.get());
        }
        super.initAndValidate();
    }

    protected void initSRanges() {
        if (stratigraphicRangeInput.get() != null) {
            sRanges = (ArrayList) stratigraphicRangeInput.get();
            List<Node> externalNodes = getExternalNodes();
            for (StratigraphicRange range:sRanges) {
                range.removeAllNodeNrs();
                for (Node node:externalNodes) {
                    if(node.getID().equals(range.getFirstOccurrenceID()) && !range.isSingleFossilRange()) {
                        if (!node.isDirectAncestor()) {
                            throw new RuntimeException("The first occurrence always has to be a sampled ancestor but " +
                                    range.getFirstOccurrenceID() + " is not a sampled ancestor. Something went wrong in " +
                                    "initializing the stratigraphic range tree."  );
                        }
                        range.setFirstOccurrenceNodeNr(node.getParent().getNr());
                    }
                    if (node.getID().equals(range.getLastOccurrenceID())) {
                        if (node.isDirectAncestor()) {
                            range.setLastOccurrenceNodeNr(node.getParent().getNr());
                        } else {
                            range.setLastOccurrenceNodeNr(node.getNr());
                        }
                    }
                }
                range.initAndValidate();
            }
        } else {
            sRanges = new ArrayList<>();
            ArrayList<StratigraphicRange> firstRanges = new ArrayList<>();
            ArrayList<StratigraphicRange> lastRanges = new ArrayList<>();
            for (Node node:getExternalNodes()) {
                String ID = node.getID();
                String IDwithoutPrefix = ID;
                String prefix = "";
                int i=ID.lastIndexOf('_');
                if (i>0) {
                    IDwithoutPrefix = ID.substring(0,i);
                    prefix = ID.substring(i+1);
                }
                if (prefix.equals("first")) {
                    if (!node.isDirectAncestor()) {
                        throw new RuntimeException("The first occurrence always has to be a sampled ancestor but " +
                                node.getID() + " is not a sampled ancestor.");
                    }
                    boolean found = false;
                    for (StratigraphicRange candidateRange:lastRanges) {
                        if (candidateRange.getID().equals(IDwithoutPrefix)) {
                            candidateRange.setFirstOccurrenceID(ID);
                            candidateRange.setFirstOccurrenceNodeNr(node.getParent().getNr());
                            sRanges.add(candidateRange);
                            lastRanges.remove(candidateRange);
                            found=true;
                            break;
                        }
                    }
                    if (!found) {
                        StratigraphicRange range = new StratigraphicRange();
                        range.setID(IDwithoutPrefix);
                        range.setFirstOccurrenceID(ID);
                        if (node.isDirectAncestor()) {
                            range.setFirstOccurrenceNodeNr(node.getParent().getNr());
                        } else {
                            range.setFirstOccurrenceNodeNr(node.getNr());
                        }
                        firstRanges.add(range);
                    }
                } else {
                    boolean found = false;
                    for (StratigraphicRange candidateRange:firstRanges) {
                        if (candidateRange.getID().equals(IDwithoutPrefix)) {
                            if (!prefix.equals("last")) {
                                throw new RuntimeException("Taxa " + candidateRange.getFirstOccurrenceID() + " and " +
                                        ID  + " are found in the tree. If " + ID + " is the last occurrence then add " +
                                        "_last at the end.");
                            }
                            candidateRange.setLastOccurrenceID(ID);
                            candidateRange.setLastOccurrenceNodeNr(node.getNr());
                            sRanges.add(candidateRange);
                            firstRanges.remove(candidateRange);
                            found=true;
                            break;
                        }
                    }
                    if (!found) {
                        StratigraphicRange range = new StratigraphicRange();
                        range.setID(IDwithoutPrefix);
                        range.setLastOccurrenceID(ID);
                        range.setLastOccurrenceNodeNr(node.getNr());
                        lastRanges.add(range);
                    }
                }
            }
            if (!firstRanges.isEmpty()) {
                throw new RuntimeException("There are taxa with first occurrence only " + firstRanges.toString() +". " +
                        "Single fossil ranges can not have _first at the end." );
            }
            for (StratigraphicRange range:lastRanges) {
                range.makeSingleFossilRange();
            }
            sRanges.addAll(lastRanges);
        }

        initStoredRanges();
    }

    public void initStoredRanges() {
        storedSRanges = new ArrayList<>();
        for (int i=0; i<sRanges.size(); i++) {
            StratigraphicRange range_src = sRanges.get(i);
            StratigraphicRange range_sink = new StratigraphicRange();
            ArrayList<Integer> nodeNrs_src = (ArrayList) range_src.getNodeNrs();
            for (int j=0; j<nodeNrs_src.size(); j++) {
                range_sink.addNodeNr(nodeNrs_src.get(j));
            }
            range_sink.setFirstOccurrenceID(range_src.getFirstOccurrenceID());
            range_sink.setLastOccurrenceID(range_src.getLastOccurrenceID());
            storedSRanges.add(range_sink);
        }

    }

    /**
     * copy of all values from existing tree *
     */
    @Override
    public void assignFrom(final StateNode other) {
        super.assignFrom(other);
        initSRanges();
    }

    /**
     * as assignFrom, but only copy tree structure *
     */
    @Override
    public void assignFromFragile(final StateNode other) {
        super.assignFromFragile(other);
        if(stratigraphicRangeInput.get()!= null) {
            initSRanges();
        }
    }


    /**
     * StateNode implementation *
     */
    @Override
    protected void store() {

        super.store();

        for (StratigraphicRange range_src:sRanges) {
            int index = sRanges.indexOf(range_src);
            StratigraphicRange range_sink = storedSRanges.get(index);
            range_sink.removeAllNodeNrs();
            for (int i=0; i< range_src.getNodeNrs().size(); i++) {
                int nodeNr = range_src.getNodeNrs().get(i);
                range_sink.addNodeNr(nodeNr);
            }
        }
    }

    @Override
    public void restore() {

        super.restore();

        ArrayList<StratigraphicRange> tmp_ranges = storedSRanges;
        storedSRanges = sRanges;
        sRanges=tmp_ranges;
    }


    // SRange methods:

    public ArrayList<StratigraphicRange> getSRanges() {
        return  sRanges;
    }

    public ArrayList<Integer> getSRangesInternalNodeNrs() { // nodes that do not represent the first occurrences, does not include nodes
        // from single fossil range
        ArrayList<Integer> internalNodeNrs = new ArrayList<>();
        for (StratigraphicRange range: sRanges) {
            internalNodeNrs.addAll(range.getInternalNodeNrs());
        }
        return internalNodeNrs;
    }

    public StratigraphicRange sRangesContainsID(String taxonName) {
        for (StratigraphicRange range:sRanges) {
            if (range.getFirstOccurrenceID().equals(taxonName)) {
                return range;
            }
            if (range.getLastOccurrenceID().equals(taxonName)) {
                return range;
            }
        }
        return null;
    }

    public StratigraphicRange getRangeOfNode(Node node) {
        int nodeNr = node.getNr();
        for (StratigraphicRange candidate_range:sRanges) {
            if (candidate_range.containsNodeNr(nodeNr)) {
                return candidate_range;

            }
        }
        return null;
    }



    public boolean belongToSameSRange(int node1Nr, int node2Nr) {
        for (StratigraphicRange range:sRanges) {
            if (range.containsNodeNr(node1Nr) && range.containsNodeNr(node2Nr)) {
                return true;
            }
        }
        return false;
    }

}
