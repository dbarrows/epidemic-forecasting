---
title: "Stochastic SIR Forecasting Showdown Part 1: Parameter Fitting"
author: "Dexter Barrows"
date: "18:32 04 January 2016"
#bibliography: bibliography.bib
output: 
  html_document: 
    fig_width: 9.5
#output: pdf_document 
#fig_width: 9.5
#fontsize: 12pt
---





***

# Setup

Now that we have established which methods we wish to evaluate the efficacy of for epidemic forecasting, it is prudent to see how they perform when fitting parameters for a known epidemic model. We have already seen how they perform when fitting parameters for a model with a deterministic evolution process and observation noise, but a more realistic model will have both process and observation noise.

To form such a model, we will take a deterministic SIR ODE model given by

$$
\begin{aligned}
	\frac{dS}{dt} & = - \beta S I  \\
	\frac{dI}{dt} & = \beta S I - \gamma \\
	\frac{dR}{dt} & = \gamma I,
\end{aligned}
$$

and add process noise by allowing $\beta$ to embark on a geometric random walk given by

$$
\beta_{t+1} = \exp \left( \log(\beta_{t}) + \eta (\log(\bar{\beta}) - \log(\beta_{t})) + \epsilon_{t} \right).
$$

We will take $\epsilon_{t}$ to be normally distributed with variance $\rho$ such that $\epsilon_{t} \sim \mathcal{N}(0,\rho)$. The geometric attraction term constrains the random walk, the force of which is $\eta \in [0,1]$. If we take $\eta = 0$ then the walk will be unconstrained; if we let $\eta = 1$ then all values of $\beta_t$ will be independent from the previous value (and consequently all other values in the sequence).

We can observe the effect of $\eta$ visually by setting $\eta$ to each extreme value and then to an intermediate value and plotting the results. The anchoring value $\beta_0$ was set to $\beta_0 = R_0 \cdot r / N$ where $R_0 = 3.0$, $r = 0.1$ and $N = 500$, giving $B_0 = 6 \times 10^{-4}$. The noise parameter $\rho$ was set to $0.5$.

For $\eta = 0$ we have

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png) 
Here we can see how the unconstrained walk allows the sequence to climb several orders of magnitude higher than the starting value, which is undesirable for the model formulation and unrealistic from a biological perspective.

Now if we set $\eta = 1$ we have

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png) 

and corresponding density plot

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) 

Notice that we have frequent oscillations between smaller and larger values within only a few time steps, but the range of values is desirable.

Finally choosing an intermediate value of $\eta = 0.5$ gives us

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 

and corresponding density plot

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png) 

Now we see a density plot similar in shape to the desired density, but now the geometric random walk displays dependence on previous values.

If we take the full stochastic SIR system and evolve it using an Euler stepping scheme with a step size of $h = 1/7$, for 1 step per day, we obtain the following plot

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png) 

where the solid line is the true trajectory and the dots show the data points obtained by adding in observation error defined as $\epsilon_{obvs} = \mathcal{N}(0,10)$. Note that this is only one possible realization of the system, and that different random numbers drawn from the system's random number generator will produce varying results. The plot below shows 9 additional example trajectories in addition to the one shown in the previous plot

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png) 

We can see the average behaviour of the system by taking the average over 100 trajectories. This is shown below with the 10 sample trajectories from the previous plot overlaid.

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png) 

Here the solid line is the mean behaviour, the dotted lines are the sample trajectories from the previous plot, and the grey ribbon is centre 95th quantile.


***

# Calibrating Samples

In order to compare HMCMC and IF2 we need to set up a fair and theoretically justified way to select the number of samples to draw form the HMCMC iterations and the number of particles to use for IF2. We assume that we are working with a problem that has an unknown real solution, so we use the Monte Carlo Standard Error (MCSE).

Suppose we are using a Monte-Carlo based method to obtain an estimate $\hat{\mu}_{n}$ for a quantity $\mu$, where $n$ is the number of samples. Then the Law of Large Numbers says that $\hat{\mu}_{n} \rightarrow \mu$ as $n \rightarrow \infty$. Further, the Central Limit Theorem says that the error $\hat{\mu}_{n} - \mu$ should shrink with number of samples such that $\sqrt{n} (\hat{\mu}_{n} - \mu) \rightarrow \mathcal{N}(0,\sigma^2)$ as $n \rightarrow \infty$, where $\sigma^2$ is the variance of the samples drawn.

We of course do not know $\mu$, but the above allows us to obtain an estimate $\hat{\sigma}_n$ for $\sigma$ given a number of samples $n$ as

$$
\hat{\sigma}_n = \sqrt{\frac{1}{n} \sum_{i=1}^{n} (X_i - \hat{\mu}) }
$$

which is known as the Monte Carlo Standard Error.

We can modify this formula to account for multiple variables by replacing the single variance measure sum by

$$
\Theta^* V (\Theta^*)^T
$$

where $\Theta^*$ is a row vector containing the reciprocals of the means of the parameters of interest, and $V$ is the variance-covariance matrix with respect to the same parameters. This in effect scales the variances with respect to their magnitudes and accounts for covariation between parameters in one fell swoop. We also divide by the number of parameters, yielding

$$
\hat{\sigma}_n = \sqrt{\frac{1}{n} \frac{1}{P} \Theta^* V (\Theta^*)^T }
$$

where $P$ is the number of particles.

The goal here is to then pick the number of HMCMC samples and IF2 particles to yield similar MCSE values. To do this we picked a combination of parameters for RStan that yielded decent results when applied to the stochastic SIR model specified above, calculated the resulting mean MCSE across several model fits, and isolated the expected number of IF2 particles needed to obtain the same value. This was used as a starting value to "titrate" the IF2 iterations to the same point.

The resulting values were 1000 HMCMC warm-up iterations with 2000 samples drawn post-warm-up, and 2000 IF2 particles sent through 15 passes, each method giving an approximate MCSE of 0.006.


***

# IF2 Fitting

Now we will use an implementation of the IF2 algorithm to attempt to fit the stochastic SIR model to the data from Figure [?]. The goal here is just parameter inference, but since IF2 works by applying a series on particle filters we essentially get the average system state estimates for a very small additional computational cost. Hence, we will will also look at that estimated behaviour in addition the the parameter estimates.

The code used here is a mix of R and C++ implemented using RCpp. The fitting was undertaken using $10^4$ particles with 10 IF2 passes and a cooling schedule given by a reduction in particle spread determined by $0.8^{p}$, where p is the pass number starting with 0.



The total runtime for the fitting was determined using the R `system.time()` function, and yielded


```
##    user  system elapsed 
##  22.833   0.203  23.216
```

The MLE parameter estimates, taken to be the mean of the particle swarm values after the final pass, were


```
##        R0         r        I0     sigma       eta      berr 
## 3.4212456 0.1034095 7.0442374 8.8072686 0.7757032 0.1538346
```

giving a relative error of


```
##          R0           r          I0       sigma         eta        berr 
##  0.14041519  0.03409544  0.40884749 -0.11927314  0.55140636 -0.69233078
```

From last IF2 particle filtering iteration, the mean state values from the particle swarm at each time step are shown with the true underlying state and data in the plot below

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png) 

where the solid line shows the true states, the points are the data, and the dashed line shows the state estimates from the last IF2 filtering pass.

# IF2 convergence plots

Since IF2 is an iterative algorithm where each pass through he data is expected to push the parameter estimates towards the MLE, we can see the evolution of these estimates as a function of the pass number. Plots showing this evolution are shown below for the six most critical parameters

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png) 

The horizontal axis shows the IF2 pass number. The solid black lines show the evolution of the ML estimates, the solid grey lines show the true value, and the dashed grey lines show the mean parameter estimates from the particle swarm after the final pass. 

# IF2 Densities

Of diagnostic importance are the densities of the parameter estimates given by the final parameter swarm. These are shown below.

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png) 

As before, the solid grey lines show the true parameter values and the dashed grey lines show the density means. It is worth noting that the IF2 parameters chosen were in part chosen so as to not artificially narrow these densities; a more aggressive cooling schedule and/or an increased number of passes would have resulted in much narrower densities, and indeed have the potential to collapse them to point estimates.


# HMC Fitting

We can use the Hamiltonian Monte Carlo algorithm implemented in the `Rstan` package to fit the stochastic SIR model as above. This was done with a single HMC chain of 3000 iterations with 1000 of those being warm-up iterations and a tinning value of 5.



The runtime retrieved again using R's `system.time()` shows


```
##    user  system elapsed 
## 103.947   2.009 107.199
```

which is significantly slower than either the custom IF2 algorithm or the POMP implementation.

The MLE parameter estimates, taken to be the means of the samples in the chain were


```
##         R0          r         I0      sigma        eta       berr 
## 3.36560451 0.09977558 6.26791501 8.65157364 0.46259628 0.11848018
```

giving a relative error of


```
##           R0            r           I0        sigma          eta 
##  0.121868168 -0.002244153  0.253583002 -0.134842636 -0.074807443 
##         berr 
## -0.763039636
```


# HMC Densities

The densities produced by Stan's HMC implementation are shown below:

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18-1.png) 

As before the solid grey lines show the true values and the dashed grey lines show the density means.

It is worth noting that the densities shown here represent a "true" MLE density estimate in that they represent HMC's attempt to directly sample from the parameter space according to the likelihood surface, unlike IF2 which is in theory only trying to get a ML point estimate. Further, these densities are much smoother and potentially more robust than those produced by the custom IF2 implementation.


# HMC Bootstrapping

Unlike particle particle-filtering-based approaches, HMC does not produce state estimates as a by-product of parameter fitting, but we can use information about the stochastic nodes related to the noise in the $\beta$ geometric random walk to reconstruct state estimates. The results of 100 bootstrap trajectories is shown below:

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-1.png) 

As before the solid line shows the true states, the dots show the data, the dotted line shows the average system behaviour, the dashed line shows the bootstrap mean, and the grey ribbon shows the centre 95th quantile of the bootstrap trajectories.


# Multi-trajectory parameter estimation

Here we fit the stochastic SIR model to 100 _(enough?)_ random independent trajectories using each method and examine the density of the estimates.




