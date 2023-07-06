# bethedging_stochastic
Simulation code for research on how stochastic models of bet hedger evolution defy deterministic predictions. [Link to paper].

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About This Research</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
        <li><a href="#julia-sims">Julia Stochastic Simulations</a></li>
        <li><a href="#matlab-sims">Matlab Markov Numerics</a></li>
      </ul>
    </li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## About This Research <a name="about-the-project"></a>

Bet hedging is a ubiquitous strategy for risk reduction in the face of unpredictable environmental change where a lineage lowers its variance in fitness across environments at the expense of also lowering its arithmetic mean fitness. Previous, deterministic research has quantified this trade-off using geometric mean fitness (GMF), and has found that bet hedging is expected to evolve only if it has a higher GMF than the wild-type. We introduce a novel stochastic framework that leverages both individual-based simulations and Markov chain numerics to capture the effects of stochasticity in the phenotypic distribution of diversified bet hedger offspring, environmental regime, and reproductive output. We find that modeling stochasticity can alter the sign of selection for the bet hedger when the GMF of the bet hedger is greater than, equal to, and less than the wild-type. We show that stochasticity in phenotypic distribution and in environment drive the sign of selection to differ from the deterministic prediction in different ways: phenotypic stochasticity causes bet hedging to be less beneficial than predicted, while environmental stochasticity causes bet hedging to be more beneficial than predicted. We conclude that existing, deterministic methods may not be sufficient to predict when bet hedging traits are adaptive.

### Built With <a name="built-with"></a>

Simulation code is written in the following languages and packages:

* [Julia Programming Language](https://julialang.org/)
  * [Distributions.jl Package](https://juliastats.org/Distributions.jl/stable/)
  * [Ticktock.jl Package](https://github.com/cormullion/TickTock.jl)
  * [Interpolations.jl Package](http://juliamath.github.io/Interpolations.jl/latest/)
  * [Dates.jl Package](https://docs.julialang.org/en/v1/stdlib/Dates/)
  * [Base.Threads.jl Package](https://docs.julialang.org/en/v1/base/multi-threading/)
* [MatLab Programming Language](https://www.mathworks.com/products/matlab/programming-with-matlab.html)

### Julia Stochastic Simulations <a name="julia-sims"></a>

Stochastic, individual based simulations written in Julia. We model the evolution of bet hedging in asexual populations of constant size *N* evolving in discrete non-overlapping generations under the Wright-Fisher model.

Simulation code:
* [bhpfix_sim.jl](https://github.com/mweissman97/bethedging_stochastic/blob/2427366dea02216dcdd6033350f6a0abaec06516/julia_sims/bhpfix_sim.jl) constructs a normalized probability of fixation curve as a function of population size.
  * The code takes in 4 arguments: 1) environmental stochasticity type (on 1, or off 0); 2) phenotypic stochasticity type (on 1, or off 0); 3) pSpec of the invading bet hedger; and 4) Delta GMF.
  * Used to construct Fig. 2, Fig. 3A
* [paramspace.jl](https://github.com/mweissman97/bethedging_stochastic/blob/2427366dea02216dcdd6033350f6a0abaec06516/julia_sims/paramspace.jl) extends the previous NPfix simulation across parameter space, for all values of pSpec and Delta GMF. The sign of selection for the NPFix curve at every point in parameter space is automatically classified as deleterious, neutral, sign inversion, beneficial, or unknown.
  * Used to construct Fig. 4, Supplemental Fig. 2

### Matlab Markov Numerics <a name="matlab-sims"></a>

Markov transition matrices written in the MatLab Programming Language that numerically estimate the probability of fixation.

Functions:
* [construct_markov.m](https://github.com/mweissman97/bethedging_stochastic/blob/2427366dea02216dcdd6033350f6a0abaec06516/matlab_sims/construct_markov.m) constructs a generic single generation Markovian transition matrix.
  * Function inputs: N = population size, wbh = fitness of invading population, wwt = fitness of wild-type
* [e0_meanmat.m](https://github.com/mweissman97/bethedging_stochastic/blob/2427366dea02216dcdd6033350f6a0abaec06516/matlab_sims/e0_meanmat.m) constructs an effective generation matrix for the P0E0 model
  * Function inputs: wa = fitness of specialist phenotype in A, wb = fitness of specialist phenotype in B, s = Delta GMF, pS = pSpec of BH, N = population size
* [e1_meanmat.m](https://github.com/mweissman97/bethedging_stochastic/blob/2427366dea02216dcdd6033350f6a0abaec06516/matlab_sims/e1_meanmat.m) constructs an effective generation matrix for the P0E1 model
  * Function inputs: wa = fitness of specialist phenotype in A, wb = fitness of specialist phenotype in B, s = Delta GMF, pS = pSpec of BH, N = population size
* [p1_meanmat.m](https://github.com/mweissman97/bethedging_stochastic/blob/2427366dea02216dcdd6033350f6a0abaec06516/matlab_sims/p1_meanmat.m) constructs an effective generation matrix for either P1 model (P1E0 or P1E1)
  * Function Inputs: wa = fitness of specialist phenotype in A, wb = fitness of specialist phenotype in B, s = Delta GMF, pS = pSpec of BH, N = population size, emodel = environmental stochasticity on or off

Simulation code:
* [markovnpf_multistoch.m](https://github.com/mweissman97/bethedging_stochastic/blob/2427366dea02216dcdd6033350f6a0abaec06516/matlab_sims/markovnpf_multistoch.m) numerically estimates the probability of fixation for a single bet hedging mutant by finding the stationary distribution of the transition matrix.
  * Used to construct Fig. 2, Fig. 3A
* [markovmodel_genstats.m](https://github.com/mweissman97/bethedging_stochastic/blob/2427366dea02216dcdd6033350f6a0abaec06516/matlab_sims/markovmodel_genstats.m) analyzes Markov matrices at the end of the first generation and first effective generation
  * Used to construct Fig. 3B-C, Supplemental Fig. 1

<!-- CONTACT -->
## Contact <a name="contact"></a>

Maya Weissman - [@maya_weissman](https://twitter.com/maya_weissman) - maya_weissman@brown.edu
