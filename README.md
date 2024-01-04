# Beyond the (geometric) mean: stochastic models undermine deterministic predictions of bet hedger evolution
[Link to preprint]((https://www.biorxiv.org/content/10.1101/2023.07.11.548608v1)).

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About This Research</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
        <li><a href="#julia-sims">Julia Stochastic Simulations</a></li>
        <li><a href="#MATLAB-sims">R Code</a></li>
        <li><a href="#R-code">MATLAB Markov Numerics</a></li>
      </ul>
    </li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## About This Research <a name="about-the-project"></a>

Bet hedging is a ubiquitous strategy for risk reduction in the face of unpredictable environmental change where a lineage lowers its variance in fitness across environments at the expense of also lowering its arithmetic mean fitness. Classically, the benefit of bet hedging has been quantified using geometric mean fitness (GMF); bet hedging is expected to evolve if and only if it has a higher GMF than the wild-type. We build upon previous research on the effect of incorporating stochasticity in phenotypic distribution, environment, and reproduction to investigate the extent to which these sources of stochasticity will impact the evolution of real-world bet hedging traits. We utilize both individual-based simulations and Markov chain numerics to demonstrate that modeling stochasticity can alter the sign of selection for the bet hedger compared to deterministic predictions. We find that bet hedging can be deleterious at small population sizes and beneficial at larger population sizes. This non-monotonic dependence of the sign of selection on population size, known as sign inversion, exists across parameter space for both conservative and diversified bet hedgers. We apply our model to published data of bet hedging strategies to show that sign inversion exists for biologically relevant parameters in two study systems: *Papaver dubium*, an annual poppy with variable germination phenology that grows in Central England, and *Salmonella typhimurium*, a pathogenic bacteria that exhibits antibiotic persistence. Taken together, our results suggest that GMF is not enough to predict when bet hedging is adaptive.

### Built With <a name="built-with"></a>

Simulation code is written in the following languages and packages:

* [Julia Programming Language](https://julialang.org/)
  * [Distributions.jl Package](https://juliastats.org/Distributions.jl/stable/)
  * [Ticktock.jl Package](https://github.com/cormullion/TickTock.jl)
  * [Interpolations.jl Package](http://juliamath.github.io/Interpolations.jl/latest/)
  * [Dates.jl Package](https://docs.julialang.org/en/v1/stdlib/Dates/)
  * [Base.Threads.jl Package](https://docs.julialang.org/en/v1/base/multi-threading/)
* [MATLAB Programming Language](https://www.mathworks.com/products/MATLAB/programming-with-MATLAB.html)
* [R Programming Language](https://www.r-project.org/)
  * [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
  * [tidyr](https://tidyr.tidyverse.org/)
  * [ggplot2](https://ggplot2.tidyverse.org/)

### Julia Stochastic Simulations <a name="julia-sims"></a>

Stochastic, individual based simulations written in Julia. We model the evolution of bet hedging in asexual populations of constant size *N* evolving in discrete non-overlapping generations under the Wright-Fisher model. Julia code written for Release 1.6.2.

Simulation code:
* [bhpfix_sim.jl](https://github.com/mweissman97/bethedging_stochastic/blob/2427366dea02216dcdd6033350f6a0abaec06516/julia_sims/bhpfix_sim.jl) constructs a normalized probability of fixation curve as a function of population size.
  * The code takes in 4 arguments: 1) environmental stochasticity type (on 1, or off 0); 2) phenotypic stochasticity type (on 1, or off 0); 3) pSpec of the invading bet hedger; and 4) Delta GMF.
  * Used to construct Fig. 2
* [paramspace.jl](https://github.com/mweissman97/bethedging_stochastic/blob/2427366dea02216dcdd6033350f6a0abaec06516/julia_sims/paramspace.jl) extends the previous NPfix simulation on the "high risk / low risk" bet hedging strategy across parameter space, for all values of pSpec and Delta GMF. The sign of selection for the NPFix curve at every point in parameter space is classified as deleterious, neutral, sign inversion, beneficial, or unknown.
  * Used to construct Fig. 3A, 3B, S2, and S3.
* [hihi_paramspace.jl]([julia_sims/hihi_paramspace.jl](https://github.com/mweissman97/bethedging_stochastic/blob/e8101cf72d5199c734076d28a7e80926c9a36d1a/julia_sims/hihi_paramspace.jl)) applies the previously developed parameter space survey to a "high risk / high risk" model of bet hedging.
* [papaver_npfix.jl]([julia_sims/papaver_npfix.jl](https://github.com/mweissman97/bethedging_stochastic/blob/e8101cf72d5199c734076d28a7e80926c9a36d1a/julia_sims/papaver_npfix.jl)) and [salmonella_npfix.jl](https://github.com/mweissman97/bethedging_stochastic/blob/e8101cf72d5199c734076d28a7e80926c9a36d1a/julia_sims/salmonella_npfix.jl) adapt the NPfix simulations to the *Papaver dubium* and *Salmonella typhimurium* examples respectively.
  * Parameters for *Papaver dubium* germination phenology estimated from [Arthur 1973](https://www.nature.com/articles/hdy197321)
  * Parameters for *Salmonella typhimurium* antibiotic persistence estimated from [Arnoldini 2014](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001928)

### MATLAB Markov Numerics <a name="MATLAB-sims"></a>

Markov transition matrices written in the MATLAB Programming Language that numerically estimate the probability of fixation. MATLAB Code written for Release 2019B.

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
* [markovnpf_multistoch.m](https://github.com/mweissman97/bethedging_stochastic/blob/2427366dea02216dcdd6033350f6a0abaec06516/MATLAB_sims/markovnpf_multistoch.m) numerically estimates the probability of fixation for a single bet hedging mutant by finding the stationary distribution of the transition matrix.
  * Used to construct Fig. 2, Fig. 3A
* [markovmodel_genstats.m](https://github.com/mweissman97/bethedging_stochastic/blob/2427366dea02216dcdd6033350f6a0abaec06516/MATLAB_sims/markovmodel_genstats.m) analyzes Markov matrices at the end of the first generation and first effective generation
  * Used to construct Fig. 3B-C, Supplemental Fig. 1

### R *Papaver* Climate Estimates <a name="R-code"></a>

* [papaver_climate.R](https://github.com/mweissman97/bethedging_stochastic/blob/f1d5db755cb8260c03c13f89899943a4cd397101/R-code/papaver_climate.R) constructs distributions for the realized mean winter temperature, number of days below 0C, and number of days below -1.6C.
  * Central England historical climate data taken from [Parker 1992](https://www.metoffice.gov.uk/hadobs/hadcet/).
  * Threshold between "mild" and "harsh" winters were calculated for each climatic variable as the mean of the 1968 and 1966 measurmenets, which were the observed "harsh" and "mild" years respectively in [Arthur 1973](https://www.nature.com/articles/hdy197321)
  * Used to construct Fig. S1

<!-- CONTACT -->
## Contact <a name="contact"></a>

[Personal Website](https://sciencemaya.com) - [Twitter](https://twitter.com/maya_weissman) - maya_weissman@brown.edu
