# From averaging to bootstrapping: Rethinking statistical approaches for testing theories of xylem conduit widening

The scaling of xylem conduits with path length has fascinated plant physiologists for decades. This scaling relationship has been observed in hundreds of different species of terrestrial plants spanning a tremendous range of form and function. There are several mechanistic models which predict the optimal scaling exponent for this relationship, ranging from b=0 to b=0.5, and several models predict exponents which are quantitatively similar or identical. The conventional approach to evaluating empirical support for predictions from theory involves estimation of b and 95% confidence intervals (CI) from linear regression of the log-transformed equation log(D) = log(D0) + b log(L).  If the 95% CI of b include a particular value predicted by theory, then the empirical data cannot reject the predicted value, and this is interpreted as support for the theory. If the 95% CI exclude a particular value predicted by theory, this is interpreted as lack of support for the theory. However, this conventional approach could lack the sensitivity to distinguish between two similar predictions for b due to data preprocessing decisions such as data averaging or combining averaging with weighting. Further, statistical choices, such as which regression model (OLS, MA, or SMA) to use to fit the data could also substantially affect the calculated slope and 95% CIs. Here we simulate the effects on data preprocessing and statistical decisions on the calculated model coefficients in order to demonstrate that these choices are important to consider when evaluating mechanistic model predictions. 

## Getting Started

To run the code provided, you need the latest version of R and R Studio (both can be downloaded from here: https://posit.co/download/rstudio-desktop/). The code uses several packages which are listed at the very top of the code. The code is set up to automatically download any missing packages and to load them prior to running the simulations. 
