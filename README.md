This is the repo for the R package boundHTS, developed to accompany my ongoing research on constructing coherent hierarchical predictive distributions for bounded time series.

The package addresses a common limitation in probabilistic hierarchical forecasting: many existing approaches rely on Gaussian assumptions that are poorly suited to data constrained to bounded supports, such as proportions, rates, and counts. This can lead to forecasts that violate natural constraints (e.g., producing values outside [0,1]) or fail to respect distributional characteristics inherent to the data.

boundHTS provides a framework for generating hierarchical forecasts that are both coherent across aggregation levels and consistent with the support of the underlying distributions. In particular, it enables probabilistic forecasting for bounded variables while preserving interpretability and respecting structural constraints imposed by hierarchical aggregation.
