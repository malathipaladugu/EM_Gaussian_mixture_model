EM_Gaussian_mixture_model
--------------------------
EM loop updates means, covariances and mixing co-efficients until change in mean is within the tolerance of 1e-4
After convergence, the data points classified into 2 clusters from a mixture of 2 Gaussian distributions
The algorithm will iterate, updating the parameters (mus, sigmas, pi) until convergence. During this, it will print the iteration number to the console
A scatter plot showing the data points classified into two clusters by the EM algorithm, where different colors represent different clusters
