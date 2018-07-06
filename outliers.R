# List of unique recorded concentration values
stuff = c(9.249, 31.656, 60.353, 130.997, 171.600, 
          215.184, 310.251, 709.624, 1339.552)

# Inverse 4-parameter logistic function
f = function(y,a,b,c,d) c*((a-d)/(y-d)-1)^(1/b)
# and the original function for reference
g = function(x,a,b,c,d) d + (a-d)/(1 + (x/c)^b)

# Goodness of fit function for optimising
ff = function(pars) 
{
  # picks the first bunch of modelled integers with the given parameters
  # (transformed so the input parameters can be negative for optimisation purposes
  # as the final parameters MUST be positive)
  data = f(1:30, exp(pars[1]), exp(pars[2]), exp(pars[3]), exp(pars[4])); 
  # initialise the result variable
  res=0; 
  # goes through each recorded concentration value
  # and finds the closest match in the modelled data
  # (using least squares)
  for (i in 1:9) res = res + min((stuff[i] - data)^2); 
  # output result
  res
}

# (cheated and used https://www.mycurvefit.com/ to get initial parameter values)
pars = log(c(0.2877245, 0.7128381, 93362130, 70055.13))
# initialise goodness of fit value
fit = 0
O$value = Inf
# optimise until it equilbriates properly
while (fit != O$value) 
{
  fit = O$value
  O = optim(pars, ff, control = list(maxit = 100000))
  pars = O$par
}

# should get these numbers when back-transformed
# (which are different to the ones I published in the comment, oops... same result though!)
# 2.875710e-01 7.127577e-01 3.090249e+11 2.258046e+07

# reinterpreted table with the discrete values
# assuming a "detection value" of 0 where concentration is listed as 0
table = 
  rbind(
    c(NA,1 ,0 ,0 ,1 ),
    c(NA,NA,9 ,0 ,NA),
    c(NA,NA,1 ,0 ,NA),
    c(16,NA,6 ,0 ,5 ),
    c(0 ,NA,0 ,NA,0 ),
    c(NA,NA,5 ,0 ,NA),
    c(NA,NA,NA,2 ,NA),
    c(NA,0 ,7 ,3 ,0 ),
    c(0 ,NA,1 ,0 ,0 ),
    c(NA,5 ,0 ,0 ,25),
    c(NA,0 ,3 ,1 ,1 ),
    c(NA,5 ,5 ,3 ,NA)
  )

# the smaller results look roughly like a Poisson distribution
hist(table, seq(-0.5,26.5), col='red', ylim = c(0,30))
# here's a Poisson distribution to illustrate the similarity (roughly chosen values)
points(0:10, dpois(0:10, 0.4)*25, col='blue', type='p', pch=19)

# reconstruct original concentrations
conc = f(table, exp(pars[1]), exp(pars[2]), exp(pars[3]), exp(pars[4]))
# place to put on the plots
coords = col(table)+0.05*(row(table)-6.5)
# convert the NaN values (not the NA values which are genuine NAs from the data)
# to some "detection limit" as plotted in the blog post - roughly 7?
conc[is.nan(conc)] = 7
# plot values
par(lwd=2, cex=1.5)
plot(t(coords), t(conc), ylim = c(5, 2500),
        pch=19, col='black', log='y', xlab = 'Visit', ylab = 'TNF concentration')
# show detection limit
abline(h=7, lty=3)

# assume Poisson distribution of "detection values", calculate Poisson confidence intervals
# and transform to obtain values for concentrations
lower = f(qgamma(0.025, shape = table, rate = 1), exp(pars[1]), exp(pars[2]), exp(pars[3]), exp(pars[4]))
# pick an arbitrary low concentration (off graph) where detection values include 0
# causing a NaN calculation in the transform
lower[is.nan(lower)] = 1 
upper = f(qgamma(0.975, shape = 1+table, rate = 1), exp(pars[1]), exp(pars[2]), exp(pars[3]), exp(pars[4]))

# plot confidence intervals
arrows(t(coords), t(lower), t(coords), t(upper), length=0.1, angle = 90, code = 3)
