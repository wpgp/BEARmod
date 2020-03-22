# BEARmod
Basic Epidemic, Activity, and Response COVID-19 model

This model implements a basic SEIR simulation model, accounting for variable daily movement patterns, recovery rates, and contact rates. The timestep defaults to daily, but can be changed to better approximate continuous time. See "bearmod_fx.R" and "updates and notes.txt" for more details for now.

For a placeholder dummy dataset and example simulation run, please see "run_model_small.R", which uses a dummy movement dataset "testmove.csv"

# Preliminary model text
(Directly pasted from working manuscript, please forgive tone changes)

During each timestep, infected people first recovered or were removed at an average rate r, where r was equal to the inverse of the average infectious period, and removal represents self-isolation and effective removal from the population as a potential transmitter of disease. Explicitly, this was incorporated as a Bernoulli trial for each infected person with a probability of recovering 1-exp⁡(-r). Then, the model converted exposed people to infectious by similarly incorporating a Bernoulli trial for each exposed individual, where the daily probability of becoming infectious 1-exp⁡(-ε), where ε was the inverse of the average time spent exposed but not infectious. Finally, to end the exposure, infection, and recovery step of the model, newly exposed people were calculated for each city based on the number of infectious people in the city I_i, and the average number of daily contacts that lead to transmission that each infectious person has c. We simulated the number of exposed in a patch on a given day through a random draw from a Poisson distribution for each infectious person where the mean number of new infections per person was c, which was then multiplied by the fraction of people in the city that were susceptible.

The infection processes within each patch therefore approximate the following deterministic, continuous-time model, where c and r varied through time:

dS/dt=S-c SI/N

dE/dt=c SI/N-εE

dI/dt=εE-rI

dR/dt=rI


Movement
After the model completed the infection-related processes, we moved infectious people between cities. To do this, we moved infected people from their current location to each possible destination (including remaining in the same place) using Bernoulli trials for each infection person, and each possible destination city. Because movement was calculated independently for each possible destination, it was possible for the total number of movers and stayers to be greater or fewer than the number of infectious people in the patch. To account for this, we rescaled the total numbers of people calculated to move and stay in their resident patch using the total number of infectious people in that patch.

Through this model, stochasticity in the numbers and places where COVID-19 appears between simulation runs in this model through variance in numbers of people becoming exposed, infectious, and removed/recovered, as well as variance in numbers of people moving from one city to another.

<a rel="license" href="http://creativecommons.org/licenses/by/3.0/us/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/3.0/us/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/3.0/us/">Creative Commons Attribution 3.0 United States License</a>.

Contact:
Nick W Ruktanonchai; 
nrukt00 at gmail.com
