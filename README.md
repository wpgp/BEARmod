# BEARmod
Basic Epidemic, Activity, and Response COVID-19 model

This model implements a basic SEIR simulation model, accounting for variable daily movement patterns, recovery rates, and contact rates. Demonstration of this model can be seen in a recent Nature paper [1]

For a placeholder dummy dataset and example simulation run, please see "run_model_small.R", which uses a dummy movement dataset "testmove.csv"

## Overall model
This model is a metapopulation model of COVID-19 transmission, based on an SEIR modeling framework. Within each patch, this model follows a fairly simple SEIR framework. The primary complexities this model is designed to describe are daily movement patterns, and spatially and temporally heterogeneous reductions in movement and contact rates. Specifically, this model is particularly suited for data that generally come from mobile phone companies.

### Baseline patch-level processes
Within each patch, this model first calculates the number of infected people who recovered or were otherwise removed from the infectious population (ie. through self-isolatuion) at an average rate r, where r is equal to the inverse of the average infectious period. This is explicitly incorporated as a Bernoulli trial for each infected person with a probability of recovering 1-exp⁡(-r). 
Then, the model converts exposed people to infectious by similarly incorporating a Bernoulli trial for each exposed individual, where the daily probability of becoming infectious 1-exp⁡(-ε), where ε was the inverse of the average time spent exposed but not infectious. 
Finally, to end the exposure, infection, and recovery step of the model, newly exposed people are calculated for each city based on the number of infectious people in the city I_i, and the average number of daily contacts that lead to transmission that each infectious person has c. This model then simulates the number of newly exposed people through a random draw from a Poisson distribution for each infectious person where the mean number of new infections per person was c, which was then multiplied by the fraction of people in the patch who are susceptible.
The infection processes within each patch therefore approximate the following deterministic, continuous-time model, where c and r varied through time:
dS/dt=S-c SI/N
dE/dt=c SI/N-εE
dI/dt=εE-rI
dR/dt=rI

### Movement between patches
After completing the infection-related processes, the model moves infectious people between cities, using the proportion of people who went from each patch to each other patch measured in the input OD matrix. Infectious people are moved from their current location to each possible destination (including remaining in the same place) using Bernoulli trials for each infectied person, and each possible destination city. 
Through this model, stochasticity in the numbers and places where COVID-19 appears between simulation runs in this model through variance in numbers of people becoming exposed, infectious, and removed/recovered, as well as variance in numbers of people moving from one city to another.

## Input options and formats
Note: These parameter specifications are relevant for v 0.92, denoted at the top of the bearmod_fx.R file.

First, you will create an empty population list HPop, using InitiatePop(). This function takes as inputs:
- pat_locator: A data frame with variables "patNames", "patIDs" (numeric; sequential from 1:number of patches), and "pop" (population per patch)
- initialInf: A vector of initially infected people per patch, length equal to the number of patches
- initialExp: A vector of initially exposed people per patch, length equal to the number of patches

The initial HPop is then fed into the runSim function, which has the following inputs:
- HPop
- pat_info: This is the same as pat_locator
- movement_reduction_df: a data frame with 3 variables, "date", "name", and "relative_movement". "name" corresponds to the patNames ID for the patch, and "relative_movement" indicates the relative proportion of movement for that day--.3 means all movement for that patch in that day (both incoming and outgoing) will be 30% of the baseline value (specified in mobmat later). This is specified on a per-day basis, and does not have to be complete--any missing day/patch pairs will have 100% of the baseline movement patterns
- contact_reduction_df: a data frame with 3 variables, "date", "name", and "relative_contact". Same as movement_reduction_df except this refers to the relative contact rate within a patch for a given day--ie. .5 means half as many contacts per person
- mobmat: A data frame with variables "date", "fr_pat", "to_pat", "move_prop". fr_pat and to_pat refer to the patch IDs of the origin and destination patches (see patIDs from pat_locator), and move_prop is the proportion of people who move from each origin to each destination on the given day in "date". If stayers are not denoted (origin = destination), then the model will designate this as 1 - sum(movement elsewhere) for a given patch. 

--more parameter definitions coming soon--

Contact:
Nick W Ruktanonchai; 
nrukt00 at gmail.com

[1] Lai, S., Ruktanonchai, N.W., Zhou, L. et al. Effect of non-pharmaceutical interventions to contain COVID-19 in China. Nature (2020). https://doi.org/10.1038/s41586-020-2293-x
