# clarkia-bet-hedging

### Bet hedging is not sufficient to explain germination patterns of a winter annual plant.

### Authors

  - Gregor-Fausto Siegmund, Cornell University, <gs589@cornell.edu>, <gfsiegmund@gmail.com>
  - David A. Moeller, University of Minnesota
  - Vincent M. Eckhart, Grinnell College
  - Monica A. Geber, Cornell University

This repository contains results associated with the project on bet hedging in <i>Clarkia xantiana</i> ssp. <i>xantiana</i>.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7734839.svg)](https://doi.org/10.5281/zenodo.7734839) 2nd release, associated with resubmission of manuscript after one round of peer review.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7082021.svg)](https://doi.org/10.5281/zenodo.7082021) 1st release, associated with first submission of manuscript.
[https://github.com/gregor-fausto/clarkia-bet-hedging](https://github.com/gregor-fausto/clarkia-bet-hedging) Link to Github repository associated with the project. 
-----

### Description

The repository is a response to the following request from Maya Weissman, by email on 05/31/24:  I [...] wanted to include some of your results in an analysis I am doing. I would like to know the arithmetic mean population growth rate, variance in population growth rate, and long-term stochastic population growth rate with and without a seed bank (aka Figure 3)" and "If you could just send me the csv used to generate Figure 3, that would be hugely helpful."

To generate the requested CSV, I used the following script `scripts/006_testHypotheses/01_demographicBetHedgingTest.R'. This script was used to estimate the values used in Figure 3, save an RDS file with the estimates, and plot Figure 3. I loaded the RDS file, renamed a handful of variables to make them more legible to users (renamed renamed `lambda` to `lambda.s`, and renamed `lamba.nosb` to `lambda.s.nosb`). I then saved the CSV.

---

### File description

#### siegmundAmNat2023-demographicBetHedgingTest-estimatesFromFigure3.csv

The file contains values used in the test of the demographic patterns expected with bet hedging (Figure 3 in manuscript). Methods for calculating these values are described in the main and supplementary texts of the manuscript. 

Specifically (quoting from Methods: Analysis: Demographic test of bet hedging) "To calculate the arithmetic mean population growth rate, we used the average of annual population growth rates, lambda_a (Evans et al. 2007). We obtained values for the average population growth rate with the field estimates of germination as well as with the seed bank eliminated (g_1 = 1). We used the posterior modes for annual estimates of per capita reproductive success, Y(t), and for population-level estimates of seed mortality and germination. We assumed that each estimate for per capita reproductive success was equally likely, calculated annual population growth rates with (eq. [1]) and without (eq. [2]) a seed bank, and computed la as the average. 

"To calculate temporal variability in population growth rate, we drew 1,000 samples from the 15 years of per capita reproductive success estimates with replacement. We paired these resampled years of estimates with the population-level values for germination and seed survival rates to calculate annual population growth rates. For both the case with and the case without a seed bank, we calculated the variance of the sequence of population growth rates.

"To calculate the long-term stochastic population growth rate, we used the same sequence of population growth rates that we used to calculate temporal variability in fitness. We calculated the long-term stochastic population growth rate with the field estimates of germination as well as with the seed bank eliminated (g_1 = 1)." 

See equation 3 in the manuscript for the equation used to calculate the stochastic population growth rate. The table below describes the variables in the CSV file. The table is a wide-form presentation of the values used for Figure 3. Panel A plots lambda.a vs. lambda.a.nosb. Panel B plots var.lambda vs. var.lambda.nosb. Panel C plots lambda.s vs. lambda.s.nosb.

| variableName                      | description   |
| --------------------------------- | ------------- |
| site                              | Site acronym              |
| lambda.a          | Arithmetic mean population growth rate with a seed bank, calculated as the average of annual population growth rates and assuming that each estimate for per capita reproductive success was equally likely |
| lambda.a.nosb          | Arithmetic mean population growth rate without a seed bank, calculated as the average of annual population growth rates and assuming that each estimate for per capita reproductive success was equally likely |
| var.lambda                            | Variance in population growth rate with a seed bank, calculated by drawing 1,000 samples from the 15 years of per capita reproductive success estimates with replacement and pairing these resampled years of estimates with the population-level values for germination and seed survival rates to calculate annual population growth rates. We then calculated the variance of the sequence of population growth rates.                               |
| var.lambda.nosb                            | Variance in population growth rate without a seed bank, calculated by drawing 1,000 samples from the 15 years of per capita reproductive success estimates with replacement and pairing these resampled years of estimates with the population-level values for germination and seed survival rates to calculate annual population growth rates. We then calculated the variance of the sequence of population growth rates.                               |
| lambda.s                            | Long-term stochastic population growth rate with a seed bank, calculated by drawing 1,000 samples from the 15 years of per capita reproductive success estimates with replacement and pairing these resampled years of estimates with the population-level values for germination and seed survival rates to calculate annual population growth rates (same sequence used for calculated var.lambda). We then calculated the long-term stochastic population growth rate using equation 3 in the text.                               |
| lambda.s.nosb                            | Long-term stochastic population growth rate without a seed bank, calculated by drawing 1,000 samples from the 15 years of per capita reproductive success estimates with replacement and pairing these resampled years of estimates with the population-level values for germination and seed survival rates to calculate annual population growth rates (same sequence used for calculated var.lambda). We then calculated the long-term stochastic population growth rate using equation 3 in the text.      