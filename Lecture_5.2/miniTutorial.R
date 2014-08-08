## # Tutorial for the R package *bayou*
## The purpose of *bayou* is to fit Bayesian models of adaptive evolution to phylogenetic comparative data. Specifically, *bayou* provides a flexible framework for fitting multi-optima Ornstein-Uhlenbeck models to phylogenetic comparative data. This tutorial demonstrates some of the options for running *bayou*.

## ## Reversible-jump MCMC over regime placement
## In this example, we will fit a reversible-jump MCMC model to an included dataset (the Chelonia dataset of Jaffe et al. 2012). To start, we will specify a model in which no parameters are fixed. This will estimate the posterior of shift number, location and magnitude as well as all other parameters.

## We begin by loading the package and simulating the data. 
#require(devtools)
#install_github("bayou", username="uyedaj")

## Now load the package and an example dataset
require(bayou)

## Simulate a phylogeny under a pure birth model
tree <- sim.bdtree(1, 0, n = 100)

## Rescale the age of the tree to 1
tree$edge.length <- tree$edge.length/max(branching.times(tree))
plot(tree)
true.pars <- identifyBranches(tree, 3)
true.pars$alpha <- 2
true.pars$sig2 <- 1
true.pars$k <- 3
true.pars$ntheta <- 4
true.pars$theta <- c(0, -2, 2, 5)
true.pars$t2 <- 2:4
dat <- dataSim(true.pars, model="OU", tree)$dat

## ### Defining a prior function
## We now need to define prior function to set up our model. This can be the trickiest part. We will set half-cauchy priors for alpha and sigma^2, a normal prior for theta, a conditional Poisson for the number of shifts, and "dsb" controls how many shifts can be per branch (either 0, 1 or Inf) and the probability of a shift being on that branch. Since we have set bmax = 1 and prob = 1; we are specifying a model with a maximum of 1 shift per branch, and equal probability across all branches.
prior <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy",dsb="dsb", dk="cdpois", dtheta="dnorm"), param=list(dalpha=list(scale=1), dsig2=list(scale=1), dk=list(lambda=10, kmax=200), dsb=list(bmax=1,prob=1), dtheta=list(mean=0, sd=2)))
## The figure produced gives a rough visual of the chosen prior distributions. 

## ### Running the MCMC
## Now we are going to run the mcmc. We are going to output files to the R temporary directory. If you want to specify another director, replace "getwd()" with the path of the directory. By default, *bayou outputs to the R temporary directory. We will run a relatively short chain of only 10,000 generations.
ngen <- 20000

par(mfrow=c(2,3))
fit1 <- bayou.mcmc(tree, dat, SE=SE, model="OU", prior, ngen=ngen, new.dir=TRUE, plot.freq=2000, ticker.freq=1000)

## Most of the output is saved in a file, only an overview of the run is saved in *fit1*. 
fit1

## We can load the actual chains by running the following code:
chain <- load.bayou(fit1, save.Rdata=FALSE, cleanup=TRUE)
chain <- set.burnin(chain, 0.3)

## We can return a summary of our MCMC results by summarizing the chain. Notice the very small effective sample sizes for parameters for only 10,000 generations, need to get more like 100 for each by running the MCMC for more generations.
## How close did we get to the true values?
out <- summary(chain)

## We can visualize the chains by viewing traces for each parameter using the plotting utilities of the R package coda:
plot(chain)

## We can view where there are shifts of high probability.
par(mfrow=c(1,1))
plotSimmap.mcmc(pars2simmap(true.pars, tree)$tree, chain, burnin=0.3, circle=TRUE,fsize=0.4)

## And we can view the density of phenotypic optima and location of highly supported shifts on the phenogram. Here we show all shifts with posterior probabilities greater than *pp.cutoff = 0.3*. 
phenogram.density(tree, dat, chain=chain, burnin=0.3, pp.cutoff=0.3)

## ## Fitting models with fixed parameters
## *bayou* also allows the user to fit models with fixed parameters, including non-reversible jump models. Let's fit the true regime placements and see what happens.
startpar <- list(alpha=1, sig2=1, k=3, ntheta=4, theta=rep(0,4), sb=true.pars$sb, loc=rep(0, 3), t2=true.pars$t2)
fixed.hypothesis <- list(k=startpar$k, sb=startpar$sb)

## Now that we know how to specify fixed models, lets set up the prior function:
prior.fixed <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy",dsb="fixed", dk="fixed", dtheta="dnorm", dloc="dloc"), 
                          param=list(dalpha=list(scale=1), dsig2=list(scale=1), dtheta=list(mean=0, sd=2)),
                          fixed=fixed.hypothesis)

## Now we can run our MCMC chains using our priors. Right now, you have to specify the starting parameters. This will soon be unnecessary...
fit.fixed <- bayou.mcmc(tree, dat, SE=SE, model="OU", prior=prior.fixed, startpar=startpar, ngen=ngen, new.dir=TRUE)

chain.fixed <- load.bayou(fit.fixed, save.Rdata=FALSE, cleanup=FALSE)
chain.fixed <- set.burnin(chain.fixed, 0.3)
plot(chain.fixed)
out.fixed <- summary(chain.fixed)
phenogram.density(tree, dat, chain=chain.fixed, burnin=0.3, pp.cutoff=0.5)
plotSimmap.mcmc(tree, chain.fixed, burnin=0.3, circle=TRUE, fsize=0.5)

## ## Model comparison
## We can now compare the two model using Bayes Factors estimated using stepping stone estimation of the marginal likelihood.

## ### Estimating marginal likelihoods
## To compare models, we want to estimate the marginal likelihood, which is performed using the stepping stone algorithm in *bayou*. The method first estimates a reference function by fitting a series of curves to the posterior of an MCMC chain. The function *steppingstone* will output a graphic showing the best-fitting density functions overlaying the posterior distribution. Then the stepping stone MCMC's are run for every value in the vector *Bk*, corresponding to each value of the power posterior function ranging from the reference function (*Bk = 0*) to the posterior (*Bk = 1*). 
## To speed things up, you can specify multiple cores. Default is set to 2 cores, but if you have more you can use them, I'm going to make use
## of 5 cores.

ss <- steppingstone(Bk=seq(0,1,length.out=5), chain, tree, dat, SE=SE, prior=prior, new.dir=TRUE, ngen=ngen, cores=5)
ss
ss <- set.burnin(ss, 0.3)
plot(ss)

ss.fixed <- steppingstone(Bk=seq(0,1,length.out=5), chain.fixed, tree, dat, SE=0, startpar=startpar, prior=prior.fixed, ngen=ngen, cores=5)
ss.fixed <- set.burnin(ss.fixed, 0.3)
ss.fixed$lnr
plot(ss.fixed)

BF <- 2*(ss.fixed$lnr-ss$lnr)
BF

## ### We can also fit QG models, or halflife/Vy parameterizations with associated priors.

## Rescale tree to generation time
gen.time <- 5
years <- 20 * 10^6
QGtree <- tree
QGtree$edge.length <- QGtree$edge.length * years/gen.time

prior.QG <- make.prior(QGtree, model="QG", dists=list(dh2="dbeta", dP="dlnorm", dw2="dlnorm", dNe="dlnorm",dk="cdpois", dtheta="dnorm", dloc="dloc"), 
                          param=list(dh2=list(shape1=20, shape2=25), dP=list(meanlog=-0.25, sdlog=0.5), dw2=list(meanlog=2, sdlog=1), dNe=list(meanlog=10, sdlog=1),
                                     dk=list(lambda=10, kmax=32), dtheta=list(mean=0, sd=2)))


fit.QG <- bayou.mcmc(QGtree, dat, SE=SE, model="QG", prior=prior.QG, ngen=ngen, new.dir=TRUE)
chain.QG <- load.bayou(fit.QG, save.Rdata=FALSE, cleanup=FALSE)
chain.QG <- set.burnin(chain.QG, 0.3)
plot(chain.QG)
out.QG <- summary(chain.QG)
phenogram.density(QGtree, dat, chain=chain.QG, burnin=0.3, pp.cutoff=0.5)
plotSimmap.mcmc(QGtree, chain.QG, burnin=0.3, circle=TRUE, fsize=0.5)

