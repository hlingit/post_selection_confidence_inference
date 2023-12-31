# Post selection confidence inference

The goal of this research is to build corrected confidence interval for linear estimators if the model used in linear regression is determined by information-based model selection criteria, such as AIC, BIC, etc.

Reason why this post-selection correction is needed: classical confidence inference of OLS assumes the model is pre-determined. If the same data is used for both model selection and fitting, then the inference is made conditioning on the model selection event, which itself is random and cannot be ignored.

In the `example.R`, we showcase how to build a corrected confidence interval for a predicted mean at new data point $x_{new}$ for an OLS whose model is determined by AIC. Here we consider a linear model $Y=X\beta+\epsilon$. The sample size is $n=50$, the number of candidate predictors $p=10$, and the ground truth $\beta^*=(1,2,3,0,0,\cdots,0)$ has only three nonzero elements.

The result shows that the corrected confidence interval is wider than the uncorrected, problematic one.

In the `motivation.R`, we provide the R code to replicate the result in the motivating example. It can be seen that the uncorrected, classical confidence intervals have an empirical coverage significantly smaller than the target level, regardless of whether $\sigma$ is known or not.

In the `simulation-knownSigma.R`, we provide the R code to replicate the simulation results in sec 4 of the paper, for the case of known sigma2. It can be seen that after post-selection correction, the proposed confidence intervals have coverage close to the target level.

In the `simulation-unknownSigma.R`, we provide the R code to replicate the simulation results in sec 4 of the paper, for the case of unknown sigma2.

In the `realdata.R`, we provide the R code to replicate the US consumption example in sec 5 of the paper. The US consumption dataset we used is the `us_change.csv` file.

