# post_selection_confidence_inference

The goal of this research is to build corrected confidence interval for linear estimators if the model used in linear regression is determined by information-based model selection criteria, such as AIC, BIC, etc.

Reason why this post-selection correction is needed: classical confidence inference of OLS assumes the model is pre-determined. If the same data is used for both model selection and fitting, then the inference is made conditioning on the model selection event, which itself is random and cannot be ignored.

In the `example.R`, we showcase how to build a corrected confidence interval for a predicted mean at new data point $x_{new}$ for an OLS whose model is determined by AIC. Here we consider a linear model $Y=X\beta+\epsilon$. The sample size is $n=50$, the number of candidate predictors $p=10$, and the ground truth $\beta^*=(1,2,3,0,0,\cdots,0)$ has only three nonzero elements.

The result shows that the corrected confidence interval is wider than the uncorrected, problematic one.
