# Python-code used to generate simulate data with msprime
```
import msprime
import numpy as np

#########################################################################


def simulation(length, ne, n1, n2, n3, nb, debug=False):

    mutation_rate=1e-8
    recombination_rate=1e-8
    sample_per_pop=40

    N_e    = ne
    N_0    = n1
    N_1    = n3
    N_1_bg = n2
    N_2    = n1
    N_3    = n1
    N_3_bn = nb
    N_4    = n1
    N_5    = n2
    N_6    = n1

    population_configurations = [
        msprime.PopulationConfiguration(sample_size=sample_per_pop, initial_size = N_0),
        msprime.PopulationConfiguration(sample_size=sample_per_pop, initial_size = N_1),
        msprime.PopulationConfiguration(sample_size=sample_per_pop, initial_size = N_2),
        msprime.PopulationConfiguration(sample_size=sample_per_pop, initial_size = N_3),
        msprime.PopulationConfiguration(sample_size=sample_per_pop, initial_size = N_4),
        msprime.PopulationConfiguration(sample_size=sample_per_pop, initial_size = N_5),
        msprime.PopulationConfiguration(sample_size=sample_per_pop, initial_size = N_6)]

    
    T_m23_starts     = 4000

    T_p3_bn_ends     = 4000    
    T_p3_p4_split    = 5000

    T_m02_m20_starts = 5000

    T_p1_p2_split    = 6000

    T_p4_p5_split    = 6000

    T_m01_m10_starts = 7000

    T_p2_p4_split    = 8000

    T_m04_m40_starts = 20000

    T_p0_p4_split    = 40000

    T_p4_p6_split    = 60000

    
    demographic_events = [
        msprime.MassMigration(
            time=T_m23_starts, source = 3, destination = 2, proportion = 0.25),

        msprime.PopulationParametersChange(
            time=T_m23_starts, initial_size=N_1_bg,population=1),

        msprime.PopulationParametersChange(
            time=T_p3_bn_ends, initial_size=N_3_bn,population=3),

        msprime.MassMigration(
            time = T_p3_p4_split, source = 3, destination = 4, proportion = 1.0),

        msprime.MassMigration(
            time=T_m02_m20_starts, source = 0, destination = 2, proportion = 0.15),
        msprime.MassMigration(
            time=T_m02_m20_starts, source = 2, destination = 0, proportion = 0.15),

        msprime.MassMigration(
            time = T_p1_p2_split, source = 1, destination = 2, proportion = 1.0),

        msprime.MassMigration(
            time = T_p4_p5_split, source = 5, destination = 4, proportion = 1.0),

        msprime.MassMigration(
            time=T_m01_m10_starts, source = 0, destination = 2, proportion = 0.25),
        msprime.MassMigration(
            time=T_m01_m10_starts, source = 2, destination = 0, proportion = 0.25),

        msprime.MassMigration(
            time = T_p2_p4_split, source = 2, destination = 4, proportion = 1.0),

        msprime.MassMigration(
            time=T_m04_m40_starts, source = 0, destination = 4, proportion = 0.25),
        msprime.MassMigration(
            time=T_m04_m40_starts, source = 4, destination = 0, proportion = 0.25),

        msprime.MassMigration(
            time = T_p0_p4_split, source = 0, destination = 4, proportion = 1.0),

        msprime.MassMigration(
            time = T_p4_p6_split, source = 4, destination = 6, proportion = 1.0)
    ]
    
    dd = msprime.DemographyDebugger(
            Ne=N_e,
            population_configurations=population_configurations,
            demographic_events=demographic_events)

    if debug:
        dd.print_history() #can comment out when we are happy demography is correct

    ts = msprime.simulate(
        population_configurations = population_configurations,
        demographic_events = demographic_events,
        mutation_rate = mutation_rate,
        length = length,
        recombination_rate = recombination_rate,
        record_migrations = False)
    
    ts = msprime.mutate(
          ts, rate=mutation_rate,
         model=msprime.InfiniteSites(alphabet=msprime.NUCLEOTIDES))

    return ts

#########################################################################


length=int(10e6)
ne = 100000
n1 = 100000
n2 = 50000
n3 = 10000
nb = 25000

for i in range(1,26):
    
    ts = simulation(length, ne, n1, n2, n3, nb, False)

    ofile = "simulations/replicate"+str(i)+".vcf"

    with open(ofile, "w") as vcf_file:
        ts.write_vcf(vcf_file, 2)
```
