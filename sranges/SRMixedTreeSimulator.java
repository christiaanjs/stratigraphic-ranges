package sranges;


import beast.evolution.tree.Node;
import beast.evolution.tree.SRMixedNode;
import beast.evolution.tree.SRMixedTree;
import beast.util.Randomizer;

import java.util.*;
import java.util.stream.Collectors;

public class SRMixedTreeSimulator {

    private double x0, lambda, mu, psi, rho, beta, lambda_a, totalRate;

    private double t;
    private int nSpecies;
    private Node root;
    private List<Node> activeNodes;
    private Set<Node> symmetricNodes;
    private Map<Node, Integer> speciesMap;

    private int numberOfSimulations;

    public static SRMixedTree doSimulation(double... args){
        SRMixedTreeSimulator simulator = new SRMixedTreeSimulator(args);
        return simulator.doSimulation();

    }

    public SRMixedTreeSimulator(double[] args){
        this(args[0], args[1], args[2], args[3], args[4], args[5], args[6]);
    }

    public SRMixedTreeSimulator(double x0, double lambda, double mu, double psi, double rho, double beta, double lambda_a){
        this.x0 = x0;
        this.lambda = lambda;
        this.mu = mu;
        this.psi = psi;
        this.rho = rho;
        this.beta = beta;
        this.lambda_a = lambda_a;

        totalRate = lambda + mu + psi + lambda_a;

        numberOfSimulations = 0;

        initSimulation();
    }

    public int getNumberOfSimulations(){
        return numberOfSimulations;
    }

    public SRMixedTree doSimulation(){
        numberOfSimulations++;
        initSimulation();
        simulateFullTree();
        sampleLeaves();
        numberNodes();
        return buildTree();
    }

    private void initSimulation(){
        t = x0;
        nSpecies = 0;

        activeNodes = new ArrayList<>();
        symmetricNodes = new HashSet<>();
        speciesMap = new HashMap<>();

        root = createNode();
        assignNewSpecies(root);
        activeNodes.add(root);

    }

    public Node getRoot(){
        return root;
    }

    public int getNActiveLineages(){
        return activeNodes.size();
    }

    public double timeAtNextEvent(){
        return (t -= Randomizer.nextExponential(totalRate * getNActiveLineages()));
    }

    public Node simulateFullTree(){
        while(timeAtNextEvent() > 0 && getNActiveLineages() > 0){
            double u = Randomizer.nextDouble() * totalRate;
            double F = 0;
            if(u < (F += lambda)){
                createSpeciationEvent();
            } else if(u < (F += mu)){
                createExtinctionEvent();
            } else if(u < (F += psi)){
                createSamplingEvent();
            } else {
                createAnageneticSpeciationEvent();
            }
        }

        for(Node leaf: activeNodes){
            leaf.setHeight(0.0);
        }

        return getRoot();
    }

    public Node sampleLeaves(){
        return sampleLeaves(root);
    }

    public Node sampleLeaves(Node node){
        if(node.isDirectAncestor()){
            return node;
        } else if(node.isLeaf()){
            if(node.getHeight() == 0.0 && Randomizer.nextDouble() < rho){
                return node;
            } else {
                return null; // Extinct or unsampled lineage
            }
        } else { // Node is fake or speciation
            Node leftSampled = sampleLeaves(node.getLeft());
            Node rightSampled = sampleLeaves(node.getRight());

            if(leftSampled != null && rightSampled != null){
                node.removeAllChildren(false);
                node.setLeft(leftSampled);
                node.setRight(rightSampled);
                return node;
            } else if(leftSampled != null || rightSampled != null){
                return leftSampled != null ? leftSampled : rightSampled;
            } else {
                if(node.isFake())
                    throw new RuntimeException("Something went wrong pruning sampled ancestor");
                else
                    return null;
            }
        }
    }

    public int numberNodes(){
        return numberNodes(root, 0);
    }

    public int numberNodes(Node node, int counter){
        if(!node.isLeaf()){
            for(Node child: node.getChildren())
                counter = numberNodes(child, counter);
        }
        node.setNr(counter++);
        return counter;
    }

    private List<Node> createNewStratigraphicRange(Node node){
        LinkedList<Node> newRange = new LinkedList<>();
        newRange.add(node);
        return newRange;
    }

    public Map<Integer, List<Node>> collectStratigraphicRanges(){
        return collectStratigraphicRanges(root);
    }

    public Map<Integer, List<Node>> collectStratigraphicRanges(Node node){
        if(node.isLeaf()){ // Leaf
            int speciesId = speciesMap.get(node);
            List<Node> newRange = createNewStratigraphicRange(node);
            Map<Integer, List<Node>> map = new HashMap<>();
            map.put(speciesId, newRange);
            return map;
        } else if(node.isFake()){ // Sampled ancestor
            Node sampledAncestor = node.getRight();
            int speciesId = speciesMap.get(sampledAncestor);
            Map<Integer, List<Node>> map = collectStratigraphicRanges(node.getLeft());
            if(map.containsKey(speciesId)){
                map.get(speciesId).add(0, node); // Expects fake node
            } else {
                map.put(speciesId, createNewStratigraphicRange(node));
            }
            return map;
        } else { // Speciation event
            Map<Integer, List<Node>> leftMap = collectStratigraphicRanges(node.getLeft());
            Map<Integer, List<Node>> rightMap = collectStratigraphicRanges(node.getRight());
            leftMap.putAll(rightMap);
            return leftMap;
        }
    }

    private StratigraphicRange buildStratigraphicRange(List<Node> sortedRangeNodes){
        return new StratigraphicRange(sortedRangeNodes.stream()
                .map((Node n) -> n.getNr()).collect(Collectors.toList()));
    }

    public List<StratigraphicRange> buildStratigraphicRanges(){
        Map<Integer, List<Node>> srangeMap = collectStratigraphicRanges();
        return srangeMap.values().stream()
                .map(this::buildStratigraphicRange)
                .collect(Collectors.toList());
    }

    public SRMixedTree buildTree(){
        return new SRMixedTree(root, buildStratigraphicRanges(), new ArrayList<>(symmetricNodes));
    }

    private void createSpeciationEvent(){
        Node parent = popActiveNode();
        boolean symmetric = Randomizer.nextDouble() < beta;

        Node child1 = createNode();
        Node child2 = createNode();

        parent.setHeight(t);
        parent.addChild(child1);
        parent.addChild(child2);

        activeNodes.add(child1);
        activeNodes.add(child2);

        if(symmetric){
            symmetricNodes.add(parent);
            assignNewSpecies(parent.getLeft());
        } else {
            propagateSpecies(parent.getLeft());
        }

        assignNewSpecies(parent.getRight());
    }

    private void createExtinctionEvent(){
        Node node = popActiveNode();
        node.setHeight(t);
    }

    private void createSamplingEvent(){
        Node node = popActiveNode();
        node.setHeight(t);

        Node sample = createNode();
        Node child = createNode();

        node.setRight(sample);
        node.setLeft(child);

        propagateSpecies(sample);
        propagateSpecies(child);

        sample.setHeight(t);

        activeNodes.add(child);
    }

    private void createAnageneticSpeciationEvent(){
        Node node = popActiveNode(); // TODO: Do we need to create an event for this?
        assignNewSpecies(node);
        activeNodes.add(node);
    }

    private void propagateSpecies(Node child){
        speciesMap.put(child, speciesMap.get(child.getParent()));
    }

    private void assignNewSpecies(Node node){
        int newSpecies = nSpecies++;
        speciesMap.put(node, newSpecies);
    }

    private Node popActiveNode(){
        return activeNodes.remove(Randomizer.nextInt(activeNodes.size()));
    }

    private Node createNode(){
        Node node = new SRMixedNode();
        return node;
    }

    public Set<Node> getSymmetricNodes(){
        return symmetricNodes;
    }

}
