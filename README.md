# SACE_PS_LMM

The R code implements the methods in the article 
"A mixed model approach to estimate the survivor average causal effect in cluster-randomized trials"

The authors develop a mixed model approach along with an expectation-maximization algorithm to estimate the survivor average causal effect (SACE) in cluster-randomized trials. We model the continuous outcome measure with a random intercept to account for intracluster correlations due to cluster-level randomization, and model the principal strata membership both with and without a random intercept. In simulations, we compare the performance of our approaches with an existing fixed-effects approach to illustrate the importance of accounting for clustering in cluster-randomized trials.

me2code: SACE estimation with random intercepts in both the membership models and the outcome models.

mecode: SACE estimation with random intercepts in the outcome models.

fecode: SACE estimation without random intercepts.
