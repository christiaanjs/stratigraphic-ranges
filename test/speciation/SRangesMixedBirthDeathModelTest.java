package test.speciation;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.tree.SRMixedTree;
import beast.evolution.tree.SRTree;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import junit.framework.TestCase;
import org.junit.Test;
import speciation.SRangesBirthDeathModel;
import speciation.SRangesMixedBirthDeathModel;
import sranges.StratigraphicRange;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by gavryusa on 24/07/17.
 */
public class SRangesMixedBirthDeathModelTest extends TestCase {

    @Test
    public void testLikelihoodNonMixedSpecialCase() throws Exception {

        ArrayList<String> taxa = new ArrayList<String>(Arrays.asList("1", "2", "3"));
        String newick = "(((((A:3.4,2_last:0.0):1.0,2_first:0.0):0.7,(B:3.5,(3_last:1.7,3_first:0.0):0.8):1.6):0.55,1_last:0.0):0.85,1_first:0.0):0.5";
        Tree tree_initial = new TreeParser(newick, false);

        StratigraphicRange sr1 = new StratigraphicRange();
        sr1.setID("range1");
        Taxon taxon1_first = new Taxon("1_first");
        Taxon taxon1_last = new Taxon("1_last");
        sr1.setInputValue("firstOccurrence", taxon1_first);
        sr1.setInputValue("lastOccurrence", taxon1_last);
        StratigraphicRange sr2 = new StratigraphicRange();
        Taxon taxon2_first = new Taxon("2_first");
        Taxon taxon2_last = new Taxon("2_last");
        sr2.setInputValue("firstOccurrence", taxon2_first);
        sr2.setInputValue("lastOccurrence", taxon2_last);
        sr2.setID("range2");
        StratigraphicRange sr3 = new StratigraphicRange();
        Taxon taxon3_first = new Taxon("3_first");
        Taxon taxon3_last = new Taxon("3_last");
        sr3.setInputValue("firstOccurrence", taxon3_first);
        sr3.setInputValue("lastOccurrence", taxon3_last);
        sr3.setID("range3");
        ArrayList<StratigraphicRange> sranges = new ArrayList<>();
        sranges.add(sr1);
        sranges.add(sr2);
        sranges.add(sr3);
        SRMixedTree tree = new SRMixedTree();
        tree.setInputValue("stratigraphicRange", sranges);
        tree.setInputValue("nodetype", "beast.evolution.tree.SRMixedNode");
        tree.assignFrom(tree_initial);

        SRangesBirthDeathModel model = new SRangesMixedBirthDeathModel();
        model.setInputValue("tree", tree);
        model.setInputValue("origin", new RealParameter("7.0"));
        model.setInputValue("birthRate", new RealParameter("1.5"));
        model.setInputValue("deathRate", new RealParameter("0.5"));
        model.setInputValue("samplingRate", new RealParameter("0.1"));
        model.setInputValue("removalProbability", new RealParameter("0.0"));
        model.setInputValue("rho", new RealParameter("0.5"));
        model.setInputValue("symmetricSpeciationProbability", new RealParameter("0.0"));
        model.setInputValue("anageneticSpeciationRate", new RealParameter("0.0"));

        model.initAndValidate();

        assertEquals(-33.29951335631795, model.calculateTreeLogLikelihood(tree), 1e-14);

    }
}
