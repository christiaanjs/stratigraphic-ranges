package speciation;

public abstract class AbstractSRangesGradientBirthDeathModel extends SRangesBirthDeathModel implements SpeciesTreeGradientDistribution {
    protected static final int PARAM_DIM = 6;

    @Override
    public String[] getParamInputNames() {
        return new String[]{
                originInput.getName(),
                birthRateInput.getName(),
                deathRateInput.getName(),
                samplingRateInput.getName(),
                rhoProbability.getName(),
                removalProbability.getName()
        };
    }
}
