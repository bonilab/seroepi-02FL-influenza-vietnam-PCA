## README FILE FOR C++ CODE

The C++ code in this folder simulates a population of 500,000 individuals experiencing influenza infections for 36 years.  The last 9 years of this simulation are taken as the years we use to output population titer data.  In main_sim.cpp, lines 126-165 show which weeks samples are taken from the simulation.  And, lines 200-250 show the main infection loop for the 36-year simulation which involves simulating infection patterns as they appeared in the Vietname influenza sentinel surveillance data from 2006 to 2014.

The person object file (person.cpp) contains all of the functions computing the probability that an individual gets infected, the antibody titer and titer waning after infection, and any boosting that occurs via the flag "person::boosting_is_on".

Lines 260-340 in main_sim.cpp show the structure of the main week-to-week infection loop.

To run the simulation, first download all the files into a single folder in a linux/unix environment.  Type 'make' to compile the C++ code into the executable file 'runsim'.  Then type './runsim' to run a single simulation.  Each simulation run will generate one "B Matrix" which is a file of 6700 individuals (one individual in each row) that shows 11 antibody titers for each individual.

The batch scripts in this folder run multiple simulations at once to carry out the validation procedures that are summarized in Figures S2 and S3 in the supplementary materials. 