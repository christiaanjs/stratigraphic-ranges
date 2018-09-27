package speciation;

import beast.evolution.speciation.SpeciesTreeDistribution;
import beast.evolution.tree.Tree;

public interface SpeciesTreeGradientDistribution {
    public double[] calculateTreeLogLikelihoodParamGradient(Tree tree);
    public String[] getParamInputNames();
}
