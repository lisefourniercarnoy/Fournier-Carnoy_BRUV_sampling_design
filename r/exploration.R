# "Using R for Modelling and Quantitative Methods in Fisheries"
# Malcolm Haddon

# 2 A Non-Introduction to R ----

## 2.2 Programming in R ----

### 2.2.3 Getting Started with MQMF ----

install.packages("MQMF")
library("MQMF")
?MQMF

### 2.2.4 Examining Code within Functions ----

plot1
lm # Lists the contents of function "lm"
methods("mean") # List all available methods to call for a S3 generic function or S3 class ("mean" in this case)

# S3 generic functions are functions that work differently depending on the S3 class of object that you give it.
getS3method("print","table") # See the code of the function "print" for an object of class "table"

ls("package:MQMF") # Lists all objects that can be exported from package "MQMF"
MQMF::plot1 # Use the function "plot1" from package "MQMF" without loading package "MQMF"

### 2.2.5 Using Functions ----

# We can  write functions to make it easier to call the code.
# Below are two functions that, when given a vector with numbers, will count how many 1's there are.
# Two different syntaxes, but same aim.
countones2 <- function(x) return(length(which(x == 1)))
countones3 <- function(x) return(length(x[x == 1]))

vect <- c(1, 2, 3, 1, 2, 3, 1, 2, 3)

countones2(vect)
countones3(vect) # Both functions work, giving us the number of 1's in vector "vect", which is 3.
# Yay!


set.seed(7100809) # Telling R to create the same set of random numbers, so it's repeatable each time we run the code
# The number "7100809" is like an ID for the set of random numbers

matdat <- matrix(trunc(runif(40)*10), nrow=5, ncol=8)
matdat # A 5x8 matrix of random numbers between 0-9
apply(matdat, 2, countones3) # Apply function "countones3" to columns ("1" for rows, "2" for columns in a matrix) of object "matdat"

plotprep2 <- function(plots=c(1,1), # Making a function that's a variation of the function "plotprep", for plots.
                      width=6,
                      height=3.75,
                      usefont=7,
                      newdev=TRUE) {
  if((names(dev.cur()) %in% c("null device", "RStudioGD")) & # Code to open plot in a new window
     (newdev)) {
    dev.new(width=width, height=height, noRStudioGD = TRUE)
  }
  par(mfrow=plots, mai=c(0.45, 0.45, 0.1, 0.05), oma=c(0, 0, 0, 0))
  par(cex=0.75, mgp=c(1.35, 0.35, 0), font.axis=usefont, font=usefont, font.lab=usefont)
}
plotprep # IDK why we just made that but yea moving on

### 2.2.6 Random Number Generation ----

# Random numbers are hard to simulate, so in R, they're pseudo-random. It's a "set" series of numbers in a loop, each number dependent on the one previous
RNGkind() # The RNGkind() function gives 3 random number generators (or 3 loops of numbers) that are very long
# You have to start the loop somewhere, called the seed. If you use the same seed, you will get the same set of pseudo-random numbers after it.
# If you don't set a seed, each new session will use a different one. But if you do, every time you open the session, the previously set seed will be used
# That's important to remember to know what your analysis is doing

# When choosing a seed, it's best not to rely on your own random-number-generation abilities, because we're naturally biased
MQMF::getseed() # This function creates a pseudo-random seed (pseudo-random-ception lol) based on your system time on the computer

# By changing the seed (below, using the two lines of code below, or set.seed() with another number), you will get different sets of numbers when using the function "rnorm()"
seed <- getseed()
set.seed(seed)
round(rnorm(5), 5) # Set A of numbers,

set.seed(123456)
round(rnorm(5), 5) # Set B of numbers,

set.seed(seed)
round(rnorm(5), 5) # Set A of numbers again, using the system time's pseudo-random seed

### 2.2.7 Printing in R ----

# You can print an output to the console, or you can export it to a named text file
sink(file="filename", split=TRUE) # This will send the outputs to a file called "filename", and also show the output onscreen (with "split=TRUE"). The file will be saved in the working directory unless otherwise specified

### 2.2.8 Plotting in R ----

# We'll use base R for plotting, although other packages can be used if ya want
library("MQMF")
data("LatA") # We're using Length-at-age data to play around with plots

setpalette("R4")
par(mfrow=c(2,2), mai=c(0.45, 0.45, 0.1, 0.05))
par(cex=0.75, mgp=c(1.35, 0.35, 0), font.axis=7, font=7, font.lab=7)

hist(LatA$age) # basic plot
hist(LatA$age, breaks=20, col=3, main="") # playing around with colours etc.

# Okay I know all this

### 2.2.9 Dealing with Factors ----

# Factors can suck sometimes
DepCat <- as.factor(rep(seq(300,600,50), 2)) # Let's create a sequence of numbers from 300 to 600, with a 50 step (seq(300,600,50)), repeated twice (rep(__), 2), and make them factors
DepCat # This gives us 7 levels

as.numeric(DepCat) # Doesn't give us the numbers within the levels, but the levels (1 to 7)
as.numeric(levels(DepCat)) # Converts the 7 levels into numerics, not the replicates

?facttonum # This will convert all DepCat values into numerics
class(facttonum(DepCat)) # Success

# Okay we knew all that


### 2.2.10 Inputting Data ----

# tidyverse and dplr to the rescue, read.csv(), read.table() etc.



## 2.3 Writing Functions  ----

# The structure of functions is always the same
functionname <- function(argument1, fun, ...) {
  # body of the function
}


#### 2.3.1 Simple Functions ####

# Let's make a Von Bertlanffy function in R in different ways
ages <- 1:20
nages <- length(ages)
Linf <- 50
K <- 0.2
t0 <- -0.75 # Add a hypothetical species' parameters

# Making the Von Bertlanffy equation for the hypothetical species using a for loop
loopLt <- numeric(nages) # (idk how this works)
for (ag in ages) loopLt[ag] <- Linf * (1-exp(-K * (ag - t0)))

# Making the Von Bertlanffy equation for the hypothetical species by vectorising every age's length
vecLt <- Linf * (1 - exp(-K * (ages - t0))) # By putting all the lengths at age in an vector "vecLt"
vecLt

# Making the Von Bertlanffy equation for the hypothetical species by making a function
pars=c(Linf, K, t0) # I'm guessing "pars" is parameters?, and you're telling it these are the parameters we're using
vB <- function(pars, inages) { # this line tells us when calling for the function, we will need to give argument "pars" and argument "inages" for it to work
  return(pars[1] * (1 - exp(-pars[2] * (inages - pars[3]))))
}
funLt = vB(c(Linf, K, t0), ages)

# These are 3 ways to arrive at the same growth curve for the hypothetical species
ans <- cbind(ages, funLt, vecLt, loopLt); ans


#### 2.3.2 Function Input Values ####

# Functions work in an insulated environment. For example, you could change a variable called "popnum" in a function to whatever you like, and it wouldn't affect a variable also called "popnum" outside the function
# When using a function, take the time to write things out for the sake of clarity
# for example, in the function above (2.3.1), we wrote:
funLt = vB(c(Linf, K, t0), ages) # to be quick,
# but we should really write:
funLt = vB(pars=c(Linf, K, t0), inages=ages) # for clarity, so you know which argument you're giving the function
# Otherwise you can get confused

#### 2.3.3 R Objects ####

# Everything that exists in R is an object
# Everything that happens in R is a function call

#### 2.3.4 Scoping of Objects ####

# As mentioned in 2.3.2, functions have their own working environments.
# This means that if you have 2 variables (one within the function and one outside) with the same name, the function will first look within itself for a variable to work with.
# In 2.3.1, one of the parameters for the function is "inages", because it's better to give a different variable name, given we already have a variable "ages" outside the function. Just for clarity's sake.
# Conversely, a variable within a function doesn't really exist outside the function.

#### 2.3.5 Function Inputs and Outputs ####

# It's good practice to format any data you use into standard format, so that functions can more easily deal with them.
data(schaef)

# Vectors, data frames and matrices behave and are referenced differently, so you have to know what is what in order to deal with things correctly.
class(schaef) # schaef is a data frame, so it's best to transform into a matrix ?
a <- schaef[1:5, 2]
b <- schaef[1:5, "catch"]
c <- schaef$catch[1:5] # 3 ways of asking the first 5 rows of the "catch" column (the 2nd column in the dataframe)
cbind(a, b, c)

mschaef <- as.matrix(schaef) # Boom, matrix
mschaef[1:5, "catch"] # Another way, from a matrix

d <- try(mschaef$catch[1:5]) # This formulation doesn't work for matrices (mschaef), using the function "try" is like a test of a function before committing to it
d

# Okay I'm not sure what the point of that section was but moving on.

# 3 Simple Population Models ----

## 3.1 Introduction ----

# Fisheries modelling and ecological population modelling is pretty much a useless differenciation
# Here we try to understand fundamental dynamics in popuation biology

### 3.1.1 The Discrete Logistic Model ----

# The Discrete Logistic Model is a classical population model in ecology, that is surprisingly effective
# It's made of 4 terms (or components) that are time-step dependent

# Nt (the number in population at time t, N0 is the initial population size)
# r (growth rate of the population)
# K (the carrying capacity of the system for the species = max number of individuals the habitat can support basically)
# Ct (catches at time t)

# The equation goes:
# N(t+1) = Nt + rNt(1-(Nt/K))-Ct

# The (1-Nt/K) element limits the population growth by K, the carrying capacity. Without it, the population would grow endlessly from the Nt + rNt part of the equation, which is not realistic
# The more the population increases, the more the population growth is limited by K. This is a density-dependent effect.

surprod <- function(Nt, r, K) return ((r*Nt)*(1-(Nt/K))) # This creates a function for the Discrete Logistic Model
densdep <- function(Nt, K) return((1-Nt/K)) # This creates a function for the density-dependence
r <- 1.2; K <- 1000.0; Nt <- seq(10, 1000, 10) # These are parameters for a hypothetical species

par(mfrow = c(2,1), mai = c(0.4, 0.4, 0.05, 0.05), oma = c(0.0,0, 0.0, 0.0))
par(cex = 0.75, mgp = c(1.35, 0.35, 0), font.axis = 7, font = 7, font.lab = 7) # Setting up plot elements

plot1(Nt, surprod(Nt, r, K), xlab = "Population Nt", defpar = FALSE,
      ylab = "Production") # This plots the production of the hypothetical population based on how big the population is
# It's super productive until K is like "okay that's too many there's no more resources for you here" so growth declines

plot1(Nt, densdep(Nt, K), xlab = "Population Nt", defpar = FALSE,
      ylab = "Density-Dependence") # This plot basically explains how The density-dependent term deals with population growth
# I kind of read it like "growth brakes". At low population sizes, all good no brakes (Density-dependence at 0.8 or 1.0), then at big populations, hit the brakes too many people here (Density-dependence at 0.2 or whatever)

### 3.1.2 Dynamic Behaviour ----

discretelogistic # This function from MQMF uses a for loop.
# For loops I absolutely do not understand but we're gonna have to into it, because:
# Population dynamics are fundamentally sequential. The number at time t depends on the number at time t-1.
# Therefore vectorisation is not efficient.


yrs <- 100; rv = 2.8; Kv <- 1000.0; Nz = 100; catch = 0.0; p = 1.0
ans <- discretelogistic(r = rv, K = Kv, N0 = Nz, Ct = catch, Yrs = yrs, p = p)
avcatch <- mean(ans[(yrs-50):yrs, "nt"], na.rm = TRUE)
label <- paste0("r=", rv, " K=", Kv, " Ct=", catch, " N0=", Nz, " p=", p)
plot(ans, main = label, cex = 0.9, font = 7) # black dots are the last 20% of points

# Okay so. Above is a population dynamics plots for parameters we establish. The right plot is a phase plot, like a pendulum - the more points, the more chaotic.
# This basically tells us what would an unfished (catch = 0) population look like over 100 years if the starting population was 100 and the carrying capacity was 1000.
# Now the "rv" is the growth rate. When we change it (try anywhere between 0.5 and 2.8), the dynamics of the population change.
# Low growth rates give us a stable population (equilibrium), with very little oscillatory behaviour in the phase plot. Bigger growth rates are a lot more chaotic. You can get cycles with growth rates in between.
# If you change parameters like Nz (initial population), you'll get different points on the phase plot, but always on the parabola.

# Sometimes the model can have unusual behaviour that isn't really reflective of nature, like negative population numbers, or numbers way above K.
#You can force the minimum to be zero and the maximum to be K, or 1.1K, depending on how variable your population is

yrs = 600 # Let's try to run discretelogistic() for 600 years
ans <- discretelogistic(r=2.55, K = 1000, N0 = 100, Ct = 0.0, Yrs = yrs) # put it in an object we can have a look at
discretelogistic # The function calls for 601 years (look for "yr1 <- Yrs + 1" in the code) because in the matrix created, each year is given an Nt+1, but the last year (600) wouldn't be attributed one if the code didn't do that. IDK why tho
ans # We can have a look at the cycles created with our parameters by looking at the population numbers, rather than the plot.

### 3.1.3 Finding Boundaries between Behaviours

# Let's try to find what the cycle limit is for this scenario
yrs <- 600
ans <- discretelogistic(r=2.55, K = 1000, N0 = 100, Ct = 0.0, Yrs = yrs) # Same as above

avt <- round(apply(ans[(yrs-100):(yrs-1), 2:3], 1, mean), 2) # "For all years between 500 and 599, give me the mean of Nt - and round them to 2 decimal places". (the ",[2:3]), 1," section says "second and third columns of "ans", but only the first of the two", so the nt column)
count <- table(avt)
count[count > 1] # Boom we have 8 numbers that repeat themselves in the last 100 years, so an 8 cycle limit.

# Let's try to look for what r value gives stable limit cycles of different periods

testseq <- seq(1.9, 2.59, 0.01) # A test sequence of r values between 1.9 and 2.59, with a 0.01 step
nseq <- length(testseq)
result <- matrix(0, nrow = nseq, ncol = 2,
                 dimnames = list(testseq, c("r", "Unique")))

yrs <- 600
for (i in 1:nseq) { # "For every i between 1 and nseq..."
  rval <- testseq[i] # "the r values are given the test sequence" (see above)
  ans <- discretelogistic(r = rval, K = 1000, N0 = 100, Ct = 0.0, Yrs = yrs) # We make a discrete logistic model
  ans <- ans[-yrs,] # We remove the last year, not sure why
  ans[, "nt1"] <- round(ans[, "nt1"], 3) # "We round Nt to 3 decimal places"
  result[i, ] <- c(rval, length(unique(tail(ans[, "nt1"], 100)))) # "And we want the number of unique population numbers of the last 100 years for each r value as a matrix
  }
class(result)

result # Okay so what we have is the number of cycles in the last 100 years of a 600-year model, for growth rates between 1.9 to 2.59
# For lower growth rates, there's a 2-cycle stable limit, then 4, then 8 as r increases. 100 implies non-equilibrium.

### 3.1.4 Classical Bifurcation Diagram of Chaos ----

# From a very simple model you can get chaotic behaviour, which is fun and exciting
# May (1976) made a diagram (bifurcation diagram) of how equilibrium changes with just one parameter changed.
# The below code recreates the diagram (scary code idk how it works)

bifurcation <- function(testseq, tail = 100, yrs = 1000, limy = 0, incx = 0.001) {
  nseq <- length(testseq)
  result <- matrix(0, nrow = nseq, ncol = 2,
                   dimnames = list(testseq, c("r", "Unique Values")))
  result2 <- matrix(NA, nrow = nseq, ncol = tail)

  for (i in 1:nseq) {
    rval <- testseq[i]
    ans <- discretelogistic(r = rval, K = 1000, N0 = 100, Ct = 0.0, Yrs = yrs)
    ans[, "nt1"] <- round(ans[, "nt1"], 4)
    result[i, ] <- c(rval, length(unique(tail(ans[, "nt1"], tail))))
    result2[i,] <- tail(ans[, "nt1"], tail)
  }

if (limy[1] == 0) limy <- c(0, getmax(result2, mult = 1.02))
parset()
plot(rep(testseq[1], tail), result2[1,], type = "p", pch = 16, cex = 0.1,
     ylim = limy, xlim = c(min(testseq)*(1-incx), max(testseq)*(1+incx)),
     xlab = "r value", yaxs = "i", xaxs = "i", ylab = "Equilibrium Numbers",
     panel.first = grid())
for (i in 2:nseq)
  points(rep(testseq[i], tail), result2[i,], pch = 16, cex = 0.1)
return(invisible(list(result = result, result2 = result2)))
}
testseq <- seq(1.9, 2.975, 0.0005)
bifurcation(testseq, limy = 0)

# ANYWAY What is shows is the number of cycles (the more cycles the more chaotic) a population would go through with varying r values.
# The more branches, the more chaotic and unpredictable the population numbers get.
# In some areas it's a little less chaotic (areas with more white)

### 3.1.5 The Effect of Fishing on Dynamics ----

# So far we've only really looked at unfished populations. Let's catch some fish, make some money and see how the population evolves
yrs = 50; Kval = 1000
nocatch <- discretelogistic(r = 2.56, K = Kval, N0 = 500, Ct = 0, Yrs = yrs) # No catch
catch50 <- discretelogistic(r = 2.56, K = Kval, N0 = 500, Ct = 50, Yrs = yrs) # Catching 50 individuals per year
catch200 <- discretelogistic(r = 2.56, K = Kval, N0 = 500, Ct = 200, Yrs = yrs) # Catching 200 individuals per year
catch300 <- discretelogistic(r = 2.56, K = Kval, N0 = 500, Ct = 300, Yrs = yrs) # Catching 300 individuals per year

plottime <- function(x, ylab) {
  yrs <- nrow(x)
  plot1(x[, "year"], x[, "nt"], ylab = ylab, defpar = FALSE)
  avB <- round(mean(x[(yrs-40):yrs, "nt"], na.rn = TRUE), 3)
  mtext(avB, side = 1, outer = F, line = -1.1, font = 7, cex = 1.0)
}

par(mfrow = c(2,2), mai = c(0.25, 0.4, 0.05, 0.05), oma = c(1.0, 0, 0.25, 0))
par(cex = 0.75, mgp = c(1.35, 0.35, 0), font.axis = 7, font = 7, font.lab = 7) # Setting up plot elements

plottime(nocatch, "Catch = 0")
plottime(catch50, "Catch = 50")
plottime(catch200, "Catch = 200")
plottime(catch300, "Catch = 300")
mtext("years", side = 1, outer = TRUE, line = -0.2, font = 7, cex = 1.0)

# Holy shit this is really cool - This shows how the population would evolve and equilibriate (? lol) with different catches per year
# With increasing yearly catches, the population stabilises (at a lower number than unfished --864 v 907-- but still really cool to see)

plotphase <- function(x, label, ymax = 0) {
  yrs <- nrow(x)
  colnames(x) <- tolower(colnames(x))
  if (ymax[1] == 0) ymax <- getmax(x[, c(2:3)])
  plot(x[, "nt"], x[, "nt1"], type = "p", pch = 16, cex = 1.0, ylim = c(0,ymax),
       yaxs = "i", xlim = c(0, ymax), xaxs = "i", ylab = "nt1", xlab = "",
       panel.first = grid(), col = "darkgray")
  begin <- trunc(yrs * 0.6) # Last 40% of years
  points(x[begin:yrs, "nt"], x[begin:yrs, "nt1"], pch = 18, col = 1, cex = 1.2)
  mtext(label, side = 1, outer = F, line = -1.1, font = 7, cex = 1.2)
}

plotphase(nocatch, "Catch = 0", ymax = 1300)
plotphase(catch50, "Catch = 50", ymax = 1300)
plotphase(catch200, "Catch = 200", ymax = 1300)
plotphase(catch300, "Catch = 300", ymax = 1300)

# The Phase plots show how increasing catches stabilise the population (fewer clusters of values)

### 3.1.6 Determinism ----

# So far the model we've playing with is deterministic. If given the same parameters, the exact same patterns will emerge.
# In nature that's unlikely to happen because of environmental variation etc.
# We can add this random variation in models as a "process error" specific to each parameter, that will vary through time t.
# So in the Discrete Logistic Model we would replace r with (r + εt) and K with (K + ξt) where ε and ξ are random errors specific to each parameter.
# Including these would lead to population trajectories that aren't smooth, and that would require multiple iterations to capture the range of possible dynamics.


## 3.2 Age-Structured Modelling Concepts ----

# So far we've used simple models that disregard age- or size-structured differences in biology and behaviour
# But in a lot of dynamics, size or age affects productivity, fishing mortality etc. so it's good to get on it

### 3.2.1 Survivorship in a Cohort ----

# Alrighty so survivorship (S) is basically just the proportion of individuals from time t that make it to time t+1
# St = N(t+1)/Nt
# or : N(t+1) = St*Nt

# If mortality rate (natural and fishing combined) is Z, S is e^-Z

yrs <- 50; yrs1 <- yrs + 1
years <- seq(0, yrs, 1)
B0 <- 1000 # Biomass, or population size, at time 0
Z <- c(0.05, 0.1, 0.2, 0.4, 0.55)
nZ <- length(Z)
Bt <- matrix(0, nrow = yrs1, ncol = nZ, dimnames = list(years, Z))
Bt[1, ] <- B0

for (j in 1:nZ) for (i in 2:yrs1) Bt[i, j] <- Bt[(i-1), j]*exp(-Z[j])
plot1(years, Bt[, 1], xlab = "Years", ylab = "Population Size", lwd = 2)
if (nZ > 1) for (j in 2: nZ) lines(years, Bt[, j], lwd = 2, col = j, lty = j)
legend("topright", legend = paste0("Z =", Z), col = 1:nZ, lwd = 3,
       bty = "n", cex = 1, lty = 1:5)

# So this is just a plot of how the population would evolve under different mortality values, if there was no growth rate
# COol beans

### 3.2.2 Instantaneous vs Annual Mortality Rates ----

# Just above is instantaneous mortality rate. As in those that die, die with time steps super small, as opposed to dying all on Dec. 31st

# Annual proportion surviving is S = e^-Z
# So the proportion that die is A = 1 - e^-Z.

# If we ignore natural mortality for a sec, the above is basically harvest rate (proportion of the population that you fish out)
# Harvest rate: H = 1- e^-F (F being fishing mortality rate)
# Move things around to get instantaneous fishing mortality F = -log(1 - H)

Fi <- seq(0.001, 2, 0.001)
H <- 1 - exp(-Fi)
parset()
plot(Fi, H, type = "l", lwd = 2, panel.first = grid(),
     xlab = "Instantaneous Fishing Mortality F",
     ylab = "Annual Proportion Mortality H")
lines(c(0,1), c(0,1), lwd = 2, lty = 2, col = 2)

# There is a difference between Instantaneous and Annual mortality rates.
# Instantaneous mortality is just a subdivision of Annual mortality into tiny time steps, which with some math, ends up not just being divided by how many time steps you do in a year
# Math is getting to me I'm tired and I don't get it.

## 3.3 Simple Yield per Recruit ----

# Once recruited to a stock, a cohort can only decrease in size.
# You can view recruitment as post-larval arrival or as becoming available to a fishery. We'll view it as the first for now.
# Yield per Recruit is a notion in fisheries

age <- 1:11; nage <- length(age); N0 <- 1000
WaA <- c(NA, 0.082, 0.175, 0.283, 0.4, 0.523, 0.7, 0.85, 0.925, 0.99, 1.0)
H <- c(0.01, 0.06, 0.11, 0.16, 0.21, 0.26, 0.31, 0.36, 0.55, 0.8)
nH <- length(H)

NaA <- matrix(0, nrow = nage, ncol = nH, dimnames = list(age, H)) # This creates a blank matrix with row and column names taken from vectors "age" and "H" respectively.
CatchN <- NaA; CatchW <- NaA # This creates duplicate blank matrices of NaA

# If you want to for-loop a matrix, you need two for-loops in one another. One to code for rows, one to code for columns ?
# You can update multiple matrices with one for-loop call. As long as they are the same dimensions, you're all good.
# Here, we're updating NaA and CatchN
# Syntax is throwing me off, but I think: "NaA[age, i]" is like "i of value age in matrix NaA"

for (i in 1:nH) { # "For every i between 1 and nH..."
  NaA[1, i] <- N0 # "...Numbers at age is N0..." (set to 1000, see above)
  for (age in 2:nage) { # "...And for every age (row) between 2 and the max age of the species" (11, see above) - Age 1 is left alone
    NaA[age, i] <- NaA[(age - 1), i] * (1- H[i]) # "Every age i in NaA is given the value of __" or "The numbers at age (2 and up) is the previous age's numbers minus the harvest rate..." IN THE MATRIX NaA
    CatchN[age, i] <- NaA[(age - 1), i] - NaA[age, i] # "The Catch numbers at age (2 and up) is the difference between last year's and this year's population" IN THE MATRIX CatchN
  }
  CatchW[, i] <- CatchN[,i] * WaA # And then for Catch weight at age, do another matrix (CatchW), multiplying catch numbers at that age with the weight at that age.

}

totC <- t(colSums(CatchW, na.rm = TRUE))

NaA # This is a matrix that tells us: for different harvest rates (columns), here is the size of the cohort of 1, 2, 3 ... 11 year old individuals in a population.
CatchN # This one tells us how many fish of each age we would catch with different catch rates
CatchW # Same but for weight

# Essentially the more intense the harvest, the fewer older fish there are.

plot1(H, totC, xlab = "Harvest Rate", ylab = "Total Yield", lwd = 2) # This lumps all age classes together
# This is a simple yield per recruit plot that ignore natural mortality or any other influences.


### 3.3.1 Selectivity in Yield-per-Recruit ----

#  Usually in age-structured models in fisheries, you don't treat every age class as being as likely to be fished as each other
# In a real life situation we'd use fishery data (length- or age-composition data )to get selectivity curves, but here we'll just give values
# The reason for including selectivity is to determine the optimum age at which to begin applying fishing mortality. (setting mesh size etc.)

ages <- seq(0, 50, 1)
sel1 <- mature(-2.650425, 0.146017, sizeage = ages)
sel2 <- mature(-6, 0.2, ages)
sel3 <- mature(-6, 0.24, ages) # The mature() function isn't explained :/

plot1(ages, sel1, xlab = "Age Yrs", ylab = "Selectivity", cex = 0.75, lwd = 2)
lines(ages, sel2, col = 2, lxd = 2, lty = 3)
lines(ages, sel3, col = 3, lxd = 2, lty = 3)

abline(v = 25, col = "gray", lty = 2)
abline(h = c(0.25, 0.5, 0.75), col = "gray", lty = 2)


### 3.3.2 The Baranov Catch Equation ----

# In 3.2.1 we looked at survivorship, but that was assuming that all ages suffered that same mortality.
# We need to consider fishing mortality at age: S(a, t) = e^-(M + s(a)Ft)
# Survivorship of age a at time t = e^-(Natural Mortality + (selectivity of age a * Fully selected fishing mortality in year t))

# There's a whole bunch of math to arrive to the equation to estimate the numbers taken in the catch.


age <- 1:12; nage <- length(age)
sa <- mature(-4, 2, age)
H <- 0.2; M <- 0.3
FF <- -log(1 - H)
Ft <- sa * FF
N0 <- 1000
out <- cbind(bce(M, Ft, N0, age), "Select" = sa)
out
# Okay so object "out" is the application of the Baranov catch equation to a population with an annual harvest rate of 0.2 and instantaneous natural mortality rate of 0.3


### 3.3.3 Growth and Weight-at-Age ----

# To go from numvers at age in the catch to weight of the catch, we multiply by weight at age
# Usually this weight-at-age is derived from the Von Bertalanffy (vB) equation
# There's an assumption that there is a power relationship between length at age and weight at age so you need to add some parameters
# I get it but I don't get the math

## 3.4 Full Yield-per-Recruit ----

# Let's pull all of this together to have a more complete Yield-per-Recruit (YPR) analysis
# YPR analysis assumes constant recruitment, so we can follow a single cohort but also capture enough detail
# We're ignoring some variation and stochasticity in the recruitment here, so interpretation beware.

# Let's do a better YPR
age <- 0:20; nage <- length(age)
laa <- vB(c(50, 0.25, -1.5), age) # Length at age
WaA <- (0.015 * laa ^ 3.0) / 1000 # Weight at age
H <- seq(0.01, 0.65, 0.05); nH <- length(H) # Harvest rates
FF <- round(-log(1 - H), 5) # Fully selected fishing mortality
N0 <- 1000
M <- 0.1
numt <- matrix(0, nrow = nage, ncol = nH, dimnames = list(age, FF))
catchN <- matrix(0, nrow = nage, ncol = nH, dimnames = list(age, FF))
as50 <- c(1, 2, 3)
yield <- matrix(0, nrow = nH, ncol = length(as50), dimnames = list(H, as50))

for (sel in 1:length(as50)) {
  sa <- logist(as50[sel], 1.0, age) # That's just
  for (harv in 1:nH) {
    Ft <- sa * FF[harv]
    out <- bce(M, Ft, N0, age)
    numt[, harv] <- out[, "Nt"]
    catchN[, harv] <- out[, "Catch"]
    yield[harv, sel] <- sum(out[, "Catch"] * WaA, na.rm = TRUE)
  }
}

plot1(H, yield[,3], xlab ="Harvest Rate", ylab = "Yield", cex = 0.75, lwd = 2)
lines(H, yield[, 2], lwd = 2, col = 2, lty = 2)
lines(H, yield[, 1], lwd = 2, col = 3, lty = 3)
# The plot is the culmination of having :
# 1. Different ages (20)
# 2. Different harvest rates (between 0.01 anf 0.65)
# 3. Different selectivities (as in fishes become fishable at ages, 1, 2, and 3)

## 3.5 Concluding Remarks ----

# Okay we're barely alive here.
# These are just simulations, in reality you'd be struggling with getting the right parameter values

# 4 Model Parameter Estimation ----

# Good god here we go

## 4.1 Introduction ----

# Fitting a model to data is very important.
# Model fitting requires:
# data from a process of interest in nature (samples, observations)
# selection of a model structure good for the objective we have (example, a linear relationship)
# selection of function to represent the expected distribution of how model predictions will differ from reality (example, assuming normal errors and constant variance)
# searching for model parameters that increase model fit to observed datapoints (example, minimising sum-of-squares)

# There are different ways of improving fit, but that depends on what criteria you use when you describe "a good fit"

# You can assume anything in a model as long as you can defend these assumptions

### 4.1.1 Optimisation ----

# Okay. In Excel, you can find the optimum parameters for a model (the parameter values that yield the best fit to the observed data).
# It does that by applying every parameter value and testing the fit (kind of like the loops we made previously, kinda kinda)
# It's pretty easy for a non-dynamic process, but gets harder when you have recruitment, different mortalities, age-structure etc. (which is what I'm gonna have to do yay)
# Taking a guess at which parameters give the best fit is not advisable ofc
# We can do the excel thing in R too. (ofc)

## 4.2 Criteria of Best Fit ----

# Model fitting can involve the reduction of the sum-of-squared residuals (reducing the difference between observed and predicted observation - why you sum the squares of the difference idk)
# Or it can be the reduction of the negative log-likelihood negLL (log likelihood is the log of how close observed and predicted values are, as a probability - reducing the negative log likelihood is just the same as increasing the likelihood, but making it more readable)
# Or we can use Bayesian methods. Bayesian methods use prior knowledge of the most likely parameters and update them with the likelihoods of any new data. It's rescaling things as they go along


## 4.3 Model Fitting in R ----

# Grid searching for optimum parameters in Excel is unworkable when you have more than 2 parameters to fit
# R has some optimiser functions to test parameter value to fit the model to observed data
?nlm # You can give the nlm() function a guess at parameter values and it'll test around until a "best fit" combination is found
# Other functions will do essentially the same thing, the difference being how they vary the values

### 4.3.1 Model Requirements ----

# Let's do some practice of the theory
# To model-fit we need:
# observations from the system of study (CPUE, age- or length-composition of catch etc.)
# an R function of the model we need to fit, with a vector of parameters
# an R function calculating the criterion of best fit (minimum least squares or minimum negLL) returning a single value which reflects input parameters and data, which can be minimised by:
# an R function for automatically optimising the value of the selected criterion of best fit (nlm() for now)

### 4.3.2 A Length-at-Age Example ----

# Fitting a model = estimate the parameters so that predictions match the observations (according to chosen best-fit criteria)
# We'll use the vB length-at-age model as an example

### 4.3.3 Alternative Models of Growth ----

# vB is a good and popular curve for fish growth, but it doesn't mean it fits all species
# Two other models of growth are the Michaelis-Menten equation and the Gompertz growth curve
# But always be careful with biological interpretation

## 4.4 Sum of Squared Residual Deviation ----

ssq() # Sum of Squared is the most common criterion of model fit.
# It's a sum of the difference between observed and predicted values (the smaller the diff, the better the fit). Squared to avoid negative values

### 4.4.1 Assumption of Least-Squares ----

# A big assumption is that the residual error has Normal distribution and equal variance
# If data is transformed, this may violate this assumption (de-normalise residuals and mess up variance) or it can help data fit the assumption (make residuals normal and variance even)

### 4.4.2 Numerical Solutions ----

# The three growth curves, see 4.3.3
vB <- function (p, ages) return (p[1]*(1 - exp(-p[2]*(ages-p[3]))))
Gz <- function (p, ages) return (p[1]*exp(-p[2]*exp(p[3]*ages)))
mm <- function (p, ages) return ((p[1]*ages)/(p[2] + ages^p[3]))

# Let's give example parameters to test ssq
pars <- c("Linf" = 27.0, "K" = 0.15, "t0" = -2.0)

# And the the ssq gives us a value of the fit of the vB  curve to the hypothetical species in object "LatA"
ssq(p = pars, funk = vB, ages = LatA$age, observed = LatA$length)
ssq

### 4.4.3 Passing Functions as Arguments to Other Functions ----

# Above we give the function "vB" to the function "ssq", which can be done by putting "funk(...)" in the code the the main function
# I'm not sure I understand the "funk" thing here.

# Beware of misspelling object names, as it won't bring up an error, but result in zero.


### 4.4.4 Fitting the Models ----

# Let's try to fit the model to the hypothetical species data
parset()
ymax <- getmax(LatA$length)
plot(LatA$age, LatA$length, type = "p", pch = 16, cex = 1.2, xlab = "Age Years", ylab = "Length cm", col = rgb(1, 0, 0.1/5), ylim = c(0, ymax), yaxs = "i", xlim = c(0,44), panel.first = grid())


# Alright so- below is the model fitting for 3 growth curve models. The outputs give you, the minimum SSQ (top) and the parameters that yield that best fit model
# nlm will use the parameters you give it, and change them slightly to test whether it improves fit or not, so it's good to give it different sets of parameters to make sure the answers it gives you are stable

age <- 1:max(LatA$age)
pars <- c(27, 0.15, -2) # vB parameters
bestvB <- nlm(f = ssq, funk = vB, observed = LatA$length, p = pars, ages = LatA$age, typsize = magnitude(pars))
# We're using the nlm function to find the best parameters of vB  growth curve to match the observed LatA datapoints, using ssq (least squares) criteria for fit

outfit(bestvB, backtran = FALSE, title, "vB"); cat("\n") # idk what this does

pars <- c(26, 0.7, -0.5)
bestGz <- nlm(f = ssq, funk = Gz, observed = LatA$length, p = pars, ages = LatA$age, typsize = magnitude(pars))
# We're using the nlm function to find the best parameters of Gz growth curve to match the observed LatA datapoints, using ssq (least squares) criteria for fit

outfit(bestGz, backtran = FALSE, title, "Gz"); cat("\n") # idk what this does

pars <- c(26, 1, 1)
bestMM1 <- nlm(f = ssq, funk = mm, observed = LatA$length, p = pars, ages = LatA$age, typsize = magnitude(pars))
# We're using the nlm function to find the best parameters of mm growth curve to match the observed LatA datapoints, using ssq (least squares) criteria for fit

outfit(bestMM1, backtran = FALSE, title, "MM1"); cat("\n") # idk what this does

pars <- c(23, 1, 1)
bestMM2 <- nlm(f = ssq, funk = mm, observed = LatA$length, p = pars, ages = LatA$age, typsize = magnitude(pars))
# We're using the nlm function to find the best parameters of mm growth curve to match the observed LatA datapoints, using ssq (least squares) criteria for fit

outfit(bestMM2, backtran = FALSE, title, "MM2"); cat("\n") # idk what this does

# Never assume that the first optimum you're given is the best one. Play around with it first

# Let's plot the parameters to see how well it fitsssss

predvB <- vB(bestvB$estimate, ages)
predGz <- Gz(bestGz$estimate, ages)
predmm <- mm(bestMM2$estimate, ages)

ymax <- getmax(LatA$length)
xmax <- getmax(LatA$age)

plot(LatA$age, LatA$length, type = "p", pch = 16, col = rgb(1, 0, 0, 1/5), cex = 1.2, xlim = c(0, xmax), ylim = c(0, ymax), yaxs = "i", xlab = "Age", ylab = "Length (cm)", panel.first = grid())
lines(ages, predvB, lwd = 2, col = 4)
lines(ages, predGz, lwd = 2, col = 1)
lines(ages, predmm, lwd = 2, col = 3)
# The above gives you the observed points of length at age for the species, and the lines for the optimised ("best" parameters from nlm) growth curves

# Depending on the data you have (how well does the data cover very young and very old fish?) the curves may change

### 4.4.5 Objective Model Selection ----

# This section talks about AIC but I don't get the point

### 4.4.6 The Influence of Residual Error Choice on Model Fit ----

# In this example we used Normal residuals as an error but what if we used log-normal errors?
# To do this we just need to log-transform predicted and observed values before calculating SSQ

pars <- c(27.25, 0.15, -3)
bestvBN <- nlm(f = ssq, funk = vB, observed = LatA$length, p = pars, ages = LatA$age, typsize = magnitude(pars), iterlim = 1000) # Normal errors
outfit(bestvBN, backtran = FALSE, title = "Normal errors"); cat("\n")

ssqL <- function(funk, observed, ...) {
  predval <- funk(...)
  return(sum((log(observed) - log(predval))^2, na.rm = TRUE))
}

bestvBLN <- nlm(f = ssqL, funk = vB, observed = LatA$length, p = pars, ages = LatA$age, typsize = magnitude(pars), iterlim = 1000) # Log-Normal errors
outfit(bestvBLN, backtran = FALSE, title = "Log-Normal errors"); cat("\n")

# Essentially the best fit parameters change slightly, so the curve (see below) changes slightly too
predvBN <- vB(bestvBN$estimate, ages)
predvBLN <- vB(bestvBLN$estimate, ages)
ymax <- getmax(LatA$length)
xmax <- getmax(LatA$age)
parset()
plot(LatA$age, LatA$length, type = "p", pch = 16, col = rgb(1, 0, 0, 1/5), cex = 1.2, xlim = c(0, xmax), ylim = c(0, ymax), yaxs = "i", xlab = "Age", ylab = "Length (cm)", panel.first = grid())
lines(ages, predvBN, lwd = 2, col = 4)
lines(ages, predvBLN, lwd = 2, col = 1)

### 4.4.7 Remarks on Initial Model Fitting ----

# You can do a lot using Sum of Squares, but it requires normally distributed residuals and constant variance, so it can be a little constraining
# For non-normal residuals and unequal variance we can use maximum likelihood


## 4.5 Maximum Likelihood ----

# The aim with maximum likelihood is to search for model parameters that maximise the total likelihood of observation
# To use this criterion for model fit, you need the model to be defined to specify probabilities of each observation
# Max likelihood does not require normal residuals

### 4.5.1 Introductory Examples ----

# Okay so. Below we're simulating data with a normal distribution.

set.seed(12345)
x <- rnorm(10, mean = 5, sd = 1) # 10 observations, mean = 5, sd = 1
avx <- mean(x) # Getting the mean of those 10 observations (gonna be close to 5 but not quite)
sdx <- sd(x) # Getting the SD of the observations (close to 1 but not quite)
L1 <- dnorm(x, mean = 5, sd = 1) # Getting the likelihood of the observed and a "predicted" set (1)
L2 <- dnorm(x, mean = avx, sd = sdx)# Getting the likelihood of the observed and a "predicted" set (2)
# What we just did is try to compare two sets of predicted values. We want to know which is closer to the observed values (which has the highest likelihood)

result <- cbind(x, L1, L2, "L2gtL1" = (L2>L1))
result <- rbind(result, c(NA, prod(L1), prod(L2), 1)) # Making a matrix comparing the two likelihoods L1 and L2
rownames(result) <- c(1:10, "product")
colnames(result) <- c("x", "original", "estimated", "est > orig")
result # Look at the last row, saying that L2 (the second set of predicted values) is closer to the observed values
# Makes sense, given that the means and sd given to L2 are closer to the actual mean and sd of the simulated observations

?dnorm # Gives likelihood of an observation given its mean and sd
?pnorm # Cumulative probability density function.
?qnorm # Gives quantiles

# Individual likelihoods can be quite large, but when multiplied together, they can get quite small
# When you get quite small things, rounding errors can happen, so to avoid that we natural-log transform likelihoods and add them together (rather than doing the product
# Instead of maximising it, we try to minimise the sum of negative log-likelihoods
# Okay I kind of get it


## 4.6 Likelihoods from the Normal Distribution ----

# Likelihoods and probabilities are kind of different
# I don't get it

### 4.6.1 Equivalence with Sum-of-Squares ----

# If you estimate optimal model parameters of data with normally distibuted residuals, you'll get the same parameters using max likelihood and least-square methods.

### 4.6.2 Fitting a Model to Data Using Normal Likelihoods ----

# Previously we had fit models using SSQ. Let's use negLL (-log likelihood)
#vB
pars <- c(Linf = 27, K = 0.15, t0 = -3, sigma = 2.5)
ansvB <- nlm(f = negNLL, funk = vB, observed = LatA$length, p = pars, ages = LatA$age, typsize = magnitude(pars))
outfit(ansvB, backtran = FALSE, title = "vB by minimum negLL")
# MM
pars <- c(a = 23, b = 1, c = 1, sigma = 3)
ansMM <- nlm(f = negNLL, funk = mm, observed = LatA$length, p = pars, ages = LatA$age, typsize = magnitude(pars))
outfit(ansMM, backtran = FALSE, title = "vB by minimum negLL")
# OKay whatever

# There's then a residual plot too but I'm too lazy


## 4.7 Log-Normal Likelihoods ----

# Normal distribution is well known, but things in exploited populations (CPUE, catch ,effort) show highly skewed distributions
# A common distribution to describe these is the Log-Normal distribution
# In the equation to find the likelihood for this distribution, we're not substracting the values of observed and expected, but dividing them by each other. So a residual of 1 would mean a perfect mach between predicted and observed

### 4.7.1 Simplification of Log-Normal Likelihoods ----

# If you want to fit models using log-normal residuals, you can use functions (like negNLL()) that are made for normal likelihoods, as long as you transform the data and predictions before analysis

### 4.7.2 Log-Normal Properties ----

# In Normal distribution, predicted mean, median and mode are the same number. Not in log-normal distribution

x <- seq(0.05, 5, 0.01)
y <- dlnorm(x, meanlog = 0, sdlog = 1.2, log = FALSE)
y2 <- dlnorm(x, meanlog = 0, sdlog = 1, log = FALSE)
y3 <- dlnorm(x, meanlog = 0, sdlog = 0.6, log = FALSE)
y4 <- dlnorm(x,0.75, 0.6, log = FALSE)

parset(plots = c(1,2))
plot(x, y3, type = "l", lwd = 2, panel.first = grid(), ylab = "Log-Normal Likelihood")
lines(x, y, lwd = 2, col =2, lty = 2)
lines(x, y2, lwd = 2, col =3, lty = 3)
lines(x, y4, lwd = 2, col =4, lty = 4)
# This is a plot of likelihood distributions for different parameter sets

plot(log(x), y3, type = "l", lwd = 2, panel.first = grid(), ylab = "Log-Normal Likelihood")
lines(log(x), y, lwd = 2, col =2, lty = 2)
lines(log(x), y2, lwd = 2, col =3, lty = 3)
lines(log(x), y4, lwd = 2, col =4, lty = 4)
# This is the likelihood distributions for the same parameters, but log-transformed

### 4.7.3 Fitting a Curve Using Log-Normal Likelihoods ----

# Let's try to fit a model of recruitment using log-normal likelihoods as a parameter fit criterion;

data(tigers) # Tiger prawn recruitment data
lbh <- function(p, biom) return(log((p[1]*biom)/(p[2] + biom)))

# Testing the best parameters for this recruitment equation and this data
pars <- c("a" = 25, "b" = 4.5, "sigma" = 0.4)
best <- nlm(negNLL, pars, funk = lbh, observed = log(tigers$Recruit), biom = tigers$Spawn, typsize = magnitude(pars))
outfit(best, backtran = FALSE, title = "Beverton-Holt Recruitment")
predR <- exp(lbh(best$estimate, tigers$Spawn))
result <- cbind(tigers, predR, tigers$Recruit/predR)

# Examine the fit of these best fit parameters
plot1(tigers$Spawn, predR, xlab = "Spawning Biomass", ylab = "Recruitment", maxy = getmax(c(predR, tigers$Recruit)), lwd = 2)
points(tigers$Spawn, tigers$Recruit, pch = 16, cex = 1.1, col = 2)


### 4.7.4 Fitting a Dynamic Model Using Log-Normal Errors ----

# Let's to the same process as above but with dynamic (time series) catch data
data(abdat)
plotspmdat(abdat)

param <- log(c(r = 0.42, K = 9400, Binit = 3400, sigma = 0.05)) # Use log-transformed parameters for increased stability in the surplus production model
obslog <- log(abdat$cpue) # Input log-transformed observation data

bestmod <- nlm(f = negLL, p = param, funk = simpspm, indat = as.matrix(abdat), logobs = obslog) # find the best fit parameters
outfit(bestmod, backtran = TRUE, title = "abdat") # back-transform estimates (why??)

predce <- simpspm(bestmod$estimate, abdat)
ymax <- getmax(c(predce, obslog))
plot1(abdat$year, obslog, type = "p", maxy = ymax, ylab = "Log(CPUE)", xlab = "Year", cex = 0.9)
lines(abdat$year, predce, lwd = 2, col = 2)
# Okay this actually blows my mind that that just happened but: Using the function simpspm (which gets the log of predicted values for the recruitment model (?)) we're getting optimum parameters to fit abdata - this gives us a curve


## 4.8 Likelihoods from the Binomial Distribution ----

# So far we've dealt with continuous variables but sometimes what you observe is binomial (y/n, 0/1)

### 4.8.1 An Example Using Binomial Likelihoods ----

# In an example of a fishery targeting male crabs, you want to know whether targeting only males negatively affects sex ratios
# In a sample we get 60 animals, 40 of which are female. Let's see the relative likelihood of observing such a sample (and answer the q of whether the male target is bad)

n <- 60 # N of individuals in the sample
p <- 0.5 # Assuming a 1:1 M:F ratio, the probability of getting a male
m <- 1:60 # How likely each trial is
binom <- dbinom(m, n, p) # Individual likelihoods
cumbin <- pbinom(m, n, p) # Cumulated distribution

plot1(m, binom, type = "h", xlab = "Number of Males", ylab = "Probability")
abline(v = which.closest(0.025, cumbin), col = 2, lwd = 2) #
# This is a plot of how likely it is to get x number of males in a 60 individual sample
# Red line is the lower bound of the 95% confidence interval (CI = we're 95% sure we should fall in this)
# Because 20 males is outside of that 95%CI, observing 20 males is unlikely to be random (so there is probably an effect of targeting males)

# Because a sex ratio of 1:1 gives us a very low likelihood of getting 20 males out of 60 individuals, Let's try to find the sex ratio that would make finding 20 males likely
# You'd expect a sex ratio of 20/60 (0.333...) but it's good to know what range of values we can have

n <- 60
m <- 20
p <- seq(0.1, 0.6, 0.001) # Changing the sex ratio between 0.1 to 0.6
lik <- dbinom(m, n, p)
plot1(p, lik, lype = "l", xlab = "Prob of 20 males", ylab = "Prob.")
abline(v = p[which.max(lik)], col = 2, lwd = 2)

# Now we want to find the optimised parameter p for this
n <- 60; m <- 20
p <- c(0.1, 0.6)
optimise(function(p){ dbinom(m, n, p)}, interval = p, maximum = TRUE) # function "optimum" is better than nlm for single parameters
# and boom, the sex ratio that makes finding 20 males the most likely is 0.3333168


### 4.8.2 Open Bay Juvenile Fur Seal Population Size ----

# Let's use a real-life example: "what's the fur seal pup population on the island?"
# Tag recapture experiment, 151 pups tagged, then on subsequent days, tagged pups were recounted (32 tagged out of 222 counted, 31 out of 181, and 29 out of 185, see column 1 of matrix "furs")
# Remember that n. tagged/population size = n. recaptured/sample size


furseal <- c(32, 222, 1020, 704, 1337, 161.53, 31, 181, 859, 593, 1125, 135.72, 29, 185, 936, 634, 1238, 154.99)
columns <- c("tagged(m)", "Sample(n)", "Population", "95%Lower", "95%Upper", "StErr")
furs <- matrix(furseal, nrow = 3, ncol = 6, dimnames = list(NULL, columns), byrow = TRUE)

# "What's the original population size ?" Let's examine potential population sizes, and how likely resighting these pups would be with them.
optsol <- matrix(0, nrow = 2, ncol = 2, dimnames = list(furs[1:2, 2], c("p", "Likelihood")))
X <- seq(525, 1850, 1) # Range of potential population sizes
p <- 151/X # A range of proportion tagged
m1 <- furs[1, 1] + 1 # Tags observed one day, with Bailey's adjustment
m2 <- furs[2, 1] + 1 # Tags observed another day, with Bailey's adjustment

n1 <- furs[1, 2] + 1 # Sample size one day with Bailey's adjustment
n2 <- furs[2, 2] + 1 # Sample size another day with Bailey's adjustment

# Bailey's adjustment accounts for the fact that we're dealing with discrete events (idk what that means)

lik1 <- dbinom(m1, n1, p) # Likelihood of the population size with the first resighting sample
optsol[1, ] <- unlist(optimise(function(p) {dbinom(m, n, p)}, p, maximum = TRUE)) # Find the best fit of parameters m, n and p to have the maximum value of p (the proportion?)
m <- furs[2, 1] + 1; n <- furs[2, 2] + 1
lik2 <- dbinom(m2, n2, p) # Likelihood of the population size with the second resighting sample
totlik <- lik1 * lik2
optsol[2, ] <- unlist(optimise(function(p) {dbinom(m, n, p)}, p, maximum = TRUE))


plot1(X, lik1, type = "l", xlab = "Total Pup Numbers", ylab = "Probability", maxy = 0.085, lwd = 2)
abline(v =  X[which.max(lik1)], lwd = 1, col = 1, lty = 1)
lines(X, lik2, lwd = 2, col = 2, lty = 3)
abline(v = X[which.max(lik2)], col = 2, lwd = 1)
# Okay so. We had two samples of resightings, and wanted to know what was the original population size that was most likely to give us these resighting numbers.
# Because the two resighting numbers were different, it gave us two different curves, the overlap of which is likely where the actual population size is.


### 4.8.3 Using Multiple Independent Samples ----

# The likelihood of a set of independent observations is the product of the likelihoods of the individual observations.
# Just above we have totlik, which is exactly that.

totlik <- totlik/sum(totlik)
cumlik <- cumsum(totlik)
plot1(X, totlik, type = "l", lwd = 2, xlab = "Total Pup Numbers", ylab = "Posterior Joint Probability")
percs <- c(X[which.closest(0.025, cumlik)], X[which.max(totlik)], X[which.closest(0.975, cumlik)])
abline(v = percs, lwd = c(1, 2, 1), col = c(2, 1, 2))

# Adding a 6-sample average of resighting observations to the plot. Apparently 'furs' last row is an average of 6 observations
m <- furs[3, 1]; n <- furs[3, 2]
lik3 <- dbinom(m, n, p)
lik4 <- lik3/sum(lik3)
lines(X, lik4, lwd = 2, col = 3, lty = 2)


### 4.8.3 Analytical Approaches ----

# Processes like maturity (mature/not mature) follow a Binomial distribution and over time, a cohort would become 100% mature.
# We can use the logistic curve for this, and use GLM with binomial distribution


## 4.9 Other Distributions ----

# Base R has some distributions that are useful in fisheries
?Distributions

## 4.10 Likelihoods from the Multinomial Distribution ----

# Binomial distribution is good for 1/0 type outcomes of an observation.
# Multinomial distribution is what we use for observations that have more than 2 discrete outcomes possible
# There's a bunch of equations for the probability density function, and the log likelihood of the multinomial distribution

### 4.10.1 Using the Multinomial Distribution ----

# Multinomial dist. is a common one to use when fitting age- or size-composition data
# Let's use an example from abalone-

cw <- 2 # 2mm size classes
mids <- seq(8, 54, cw) # Size class sequence
obs <- c(0, 0, 6, 12, 35, 40, 29, 23, 13, 7, 10, 14, 11, 16, 11, 11, 9, 8, 5, 2, 0, 0, 0, 0) # Counts of each size class
dat <- as.matrix(cbind(mids, obs))
parset()
inthist(dat, col = 2, border = 3, width = 1.8, xlabel = "Shell Length mm", ylabel = "Frequency", xmin = 7, xmax = 55)
# There are two obvious modes, which we assume belong to two distinct cohorts
# We have 5 parameters to estimate: mean and sd for each cohort, and the proportion of the total number of observations that belong to the first cohort (for the other cohort just do 1-proportion)

av <- c(18, 34.5) # A guess of the means of the two cohorts
stdev <- c(2.75, 5.75) # A guess of the sd of the two cohorts
prop1 <- 0.55 # Proportion of observations in cohort 1
n <- sum(obs)  # Total number of shells measured
cohort1 <- (n*prop1*cw)*dnorm(mids, av[1], stdev[1]) # Likelihood of cohort1 (the first set of brackets scale for the 2mm size classes)
cohort2 <- (n*(1-prop1)*cw)*dnorm(mids, av[2], stdev[2]) # Likelihood of cohort2 (the first set of brackets scale for the 2mm size classes)

lines(mids, cohort1, lwd = 2, col = 1)
lines(mids, cohort2, lwd = 2, col = 4)

# Okay the guesses of parameters seem reasonable, but the left cohort is a little off, and the allocation of observations seems biased to the left cohort too.
# Let's... optimise!!! Yay!!!

# We need two functions: one that will generate predicted numbers-per-size-class, and one that will give the likelihood of that prediction

wrapper <- function(pars, obs, sizecl, midval = TRUE) { # This will calculate likelihood
  freqf <- predfreq(pars, sum(obs), sizecl = sizecl, midval = midval) # this will predict values
  return(mnnegLL(obs, freqf))
} # function "wrapper" will need : parameters, observations and sizeclasses, to predict some values - It will then optimise the parameters to reduce multinomial neg. log likelihood

pars <- c(av, stdev, prop1)
wrapper(pars, obs = obs, sizecl = mids) # Calculate total negative log likelihood

# Okay so below we're trying to get the best parameters- We're using two methods, midpoints of size classes (2mm size classes, if the size class is length 7-9, bin it to 8) or bounds of the size class (actually using 7-9)
# I think the latter is meant to be more sensitive

bestmod <- nlm(f = wrapper, p = pars, obs = obs, sizecl = mids, midval = TRUE, typsize = magnitude(pars))
outfit(bestmod, backtran = FALSE, title = "Using Midpts"); cat("\n") # Okay so this gives us the best parameters (in order, mean of cohort1, mean of cohort2, sd of cohort1, sd of cohort 2, proportion of observations belonging to cohort 1)
# Using Midpoints of the size classes

X <- seq((mids[1] - cw/2), (tail(mids, 1) + cw/2), cw)
bestmodb <- nlm(f = wrapper, p = bestmod$estimate, obs = obs, sizecl = X, midval = FALSE, typsize = magnitude(pars))
outfit(bestmodb, backtran = FALSE, title = "Using size-class bounds") # Okay so this gives us the best parameters (in order, mean of cohort1, mean of cohort2, sd of cohort1, sd of cohort 2, proportion of observations belonging to cohort 1)
# Using size-class bounds

# Okay now plot the best fit estimate to the original data
pars <- bestmod$estimate
cohort1 <- (n*pars[5]*cw)*dnorm(mids, pars[1], pars[3])
cohort2 <- (n*(1-pars[5])*cw)*dnorm(mids, pars[2], pars[4])
parsb <- bestmodb$estimate
nedge <- length(mids) + 1
cump1 <- (n*pars[5])*pnorm(X, pars[1], pars[3])
cohort1b <- (cump1[2:nedge] - cump1[1:(nedge - 1)])
cump2 <- (n*(1 - pars[5]))*pnorm(X, pars[2], pars[4])
cohort2b <- (cump2[2:nedge] - cump2[1:(nedge - 1)])

parset()
pick <- which(mids < 28)
inthist(dat[pick,], col = 0, border = 8, width = 1.8, xmin = 5, xmax = 28, xlabel = "Shell Length mm", ylabel = "Frequency", lwd = 3)
lines(mids, cohort1, lwd = 3, col = 1, lty = 2)
lines(mids, cohort1b, lwd = 2, col = 4)
# A tiny difference between using bounds and midpoints in this case

predmid <- rowSums(cbind(cohort1, cohort2))
predbnd <- rowSums(cbind(cohort1b, cohort2b))

result <- as.matrix(cbind(mids, obs, predmid, predbnd, predbnd-predmid))
colnames(result) <- c("mids", "Obs", "Predmid","Predbnd", "Difference")
result <- rbind(result, c(NA, colSums(result, na.rm = TRUE)[2:5]))
# This is just a table of the difference between the midpoint and bounds estimates


## 4.11 Likelihoods from the Gamma Distribution ----

# Gamma distribution is less known, but increasingly used in fisheries
# It has 2 parameters, a scale parameter b and shape parameter c.
# The distribution exists only for positive values, so 0 =< x =< +inf
# There's a bunch of equations of the gamma distribution and the density function to calculate -LL for it

X <- seq(0, 10, 0.1)
dg <- dgamma(X, shape = 1, scale = 1)
plot1(X, dg, xlab = "Quantile", "Probability Density")
lines(X, dgamma(X, shape = 1.5, scale = 1), lwd = 2, col = 2, lty = 2)
lines(X, dgamma(X, shape = 2, scale = 1), lwd = 2, col = 3, lty = 3)
lines(X, dgamma(X, shape = 4, scale = 1), lwd = 2, col = 4, lty = 4)
# Different Gamma distributions, scale parameter = 1, shape parameter changes

## 4.12 Likelihoods from the Beta Distribution ----

# Beta distribution is only defined between 0 and 1, with no possibility of getting values beyond that

x <- seq(0, 1, length = 1000)
parset()
plot(x, dbeta(x, shape1 = 3, shape2 = 1), type = "l", lwd = 2, ylim = c(0,4), xlim = c(0,1), yaxs = "i", panel.first = grid(), xlab = "Variable 0 - 1", ylab = "Beta Probability Density - Scale1 = 3")
bval <- c(1.25, 2, 4, 10)

for (i in 1:length(bval))
  lines(x, dbeta(x, shape1 = 3, shape2 = bval[i]), lwd = 2, col = (i + 1), lty = c(i + 1))
# Varying Beta distributions with shape1 value of 3 and shape2 values ranging from 1 to 10

## 4.13 Bayes' Theorem ----

### 4.13.1 Introduction ----

# Bayesian statistics (based on Bayes Theorem, which describes the probability of an event based on prior knowledge of conditions that might be related to what you're trying to model) are increasingly used in fisheries
# Let's try to compare Bayesians stats and max likelihood methods

# Conditional probabilities (which Bayesian stats are) are used to describe situations where we're interested in the prob. of one thing, given that a previous thing has happened
# Probability of an event Bi happens given that a set of events A has happened before: P(Bi|A)
# P(Bi|A) = (P(A|Bi)*P(Bi))/P(A)
# But P(Bi|A) is basically the likelihood of dataset A given the model plus parameters Bi (the hypothesis), which is where likelihood and Bayesian stats meet.
# P(Bi) is new to us, as it's the probability of Bi before analysis = the prior probability of hypothesis Bi
# P(A) considers all other possible outcomes of data, which is good in a closed system, but ambitious in the open world - see book for explanation


### 4.13.2 Bayesian Methods ----

# In fisheries and ecology, to use Bayes' theorem to generate the required posterior distribution we need three things:
# 1. A list of hypotheses to be considered with the model we're trying (combination and ranges of parameters and models we're trying)
# 2. A likelihood function to calculate the probability density of the observed data given each hypothesis
# 3. A prior probability for each hypothesis (normalised so that the sum of all prior probabilities = 1)

# Identical to determining max likelihood, with the addition of 3.

# The aim of a Bayesian analysis is not to find the optimum model fit, but rather to characterise the relative probability of different possible outcomes from an analysis
# == to characterise the uncertainty about each parameter and model output.
# There may be a more probable result, but it's presented inthe context of the distribution of probabilities for all other possibilities.


### 4.13.3 Prior Probabilities ----

# There is no constraint placed on how you figure out prior probabilities. previous work on the stock or the species, or at least some useful parameter constraints like survivorship can't be >1 - are all usefil
# If there's not enough info to produce informative prior probabilities, then you set uniform prior probabilities (all assigned equal probability)
# The good thing about Bayesian stats is that in life sciences, things can be deduced fairly easily (ex. you'd reasonably expect a deep-living species to be slow-growing, because temperature etc.) so prior probabilities can be given a relative likelihoods based on that
# *but* make sure you document the basis of your priors, otherwise you risk falling into subjective opinion territory and that's no good.
# Some say non-informative priors should be the default, but they depend on what scale you use- a uniform Prob. Density on a linear scale won't be linear on the log-scale

## 4.14 Concluding Remarks ----

# There is no clear winner with what method to use, it's best to use multiple. (that way if big differences arise, you'll know to investigate further)




# 5 Static Models ----

## 5.1 Introduction ----

# This book works on simple population models, which are needed to understand more complex models (age-structured models)
# Static models (the "simple" models) are used to describe functional form of processes like how maturity changes with age or size, stock-recruitment relationship, and selectivity of fishing gear for fished stocks.
# These relationships are given the big assumption of remaining fixed through time

# Dynamic models though try to describe processes like population size, with value at time t is estimated using the value at t-1


## 5.2 Productivity Parameters ----

# A population is made of a collection of individuals, and thus individuals partially summarise the properties of the population


## 5.3 Growth ----

# Ignoring immigration/emigration, stock production = individual growth + recruitment of new individuals
# In this chapter we'll focus on seasonal growth models, and estimating individual growth from tagging data

### 5.3.1 Seasonal Growth Curves ----

# Growth can be seasonal with environmental variables affecting metabolism.
# Pitcher & Macdonald (1973) tried to add a sine wave into the growth rate element of the vB curve to take that into consideration

data(minnow) ; week <- minnow$week; length <- minnow$length
pars <- c(75, 0.1, -10, 3.5); label = c("Linf", "K", "t0", "sigma")
bestvB <- nlm(f = negNLL, p = pars, funk = vB, ages = week, observed = length, typsize = magnitude(pars))
predL <- vB(bestvB$estimate, 0:160)
outfit(bestvB, backtran = FALSE, title = "Non-Seasonal vB", parnames = label)
# In this best fit estimation of parameters, we're getting really big sigma - this is a remnant of the fact that we told R to fit a (straight) line through seasonal (oscillating) growth data, so residuals are huge.

parset(plots = c(2, 1), margin = c(0.35, 0.45, 0.02, 0.05))
plot1(week, length, type = "p", cex = 1, col = 2, xlab = "Weeks", pch = 16, ylab = "Length", defpar = F)
lines(0:160, predL, lwd = 2, col = 1)

resids <- length - vB(bestvB$estimate, week)
plot1(week, resids, type = "l", col = "darkgrey", cex = 0.9, lwd = 2, defpar = F)
points(week, resids, pch = 16, cex = 1.1, col = "red")
abline(h = 0, col = 1, lwd = 1)

# We can reduce that by modifying the vB function to add the sine wave
svb <- function(p, ages, inc = 52) {
  return(p[1]*(1-exp(-(p[4] * sin(2*pi*(ages - p[5])/inc) + p[2] * (ages - p[3])))))
}

spars <- c(bestvB$estimate[1:3], 0.1, 5, 2)
bestsvb <- nlm(f = negNLL, p = spars, funk = svb, ages = week, observed = length, typsize = magnitude(spars))
predLs <- svb(bestsvb$estimate, 0:160)
outfit(bestsvb, backtran = F, title = "Seasonal Growth", parnames = c("Linf", "K", "t0", "C", "s", "sigma"))
# Seasonal adjustment has minor effects on Linf and K, but a bigger effect on t0 and definitely the sigma parameter
# Model fit is also improved (-LL from 150 to 105), reflected by the reduction in residual values (see below)

plot1(week, length, type = "p", cew = 0.9, col = 2, xlab = "Weeks", pch = 16, defpar = F)
lines(0:160, predLs, lwd = 2, col = 1)

resids <- length - svb(bestsvb$estimate, week)
plot1(week, resids, type = "l", col = "darkgray", cex = 0.9, xlab = "Weeks", lty = 3, ylab = "Normal Residuals", defpar = F)
points(week, resids, pch = 16, cex = 1.1, col ="red")
abline(h = 0, col = 1, lwd = 1)
# This is so friggin coooool

# Although usually at population level a non-sine wave curve is good enough for approximation of growth


### 5.3.2 Fabens Method with Tagging Data ----

# Sometimes estimating the age of a species to get length-at-age is very difficult
# One way to bypass that is to tag-recapture and use the difference in length to piece together a function
# Then we optimise to get best fit parameters for the data we have.

data(blackisland); bi <- blackisland
parset()
plot(bi$l1, bi$dl, xlab = "initial length", ylab = "growth increment after recapture (mm)")

### 5.3.3 Fitting Models to Tagging Data ----

# Okay we're trying to find the best curve to describe the above data.
# We'll use2 functions, fabens() and invl() (inverse logistic) to generate predicted growth increment,
# then nlm() to find the optimum parameters using negNLL() as a criteria for fit

linm <- lm(bi$dl ~ bi$l1)
param <- c(170, 0.3, 4); label <- c("Linf", "K", "sigma")
modelvB <- nlm(f = negNLL, p = param, funk = fabens, observed = bi$dl, indat = bi, initL = "l1", delT = "dt")
outfit(modelvB, backtran = F, title = "vB", parnames = label)
predvB <- fabens(modelvB$estimate, bi)
cat("\n")

param2 <- c(25, 130, 35, 3)
label2 = c("MaxDL", "L50", "delta", "sigma")
modelil <- nlm(f = negNLL, p = param2, funk = invl, observed = bi$dl, indat = bi, initL = "l1", delT = "dt")
outfit(modelil, backtran = F, title = "IL", parnames = label2)
predil <- invl(modelil$estimate, bi)

plot(bi$l1, bi$dl, ylim = c(-2, 31), ylab = "Growth increment mm", xlab = "Length mm")
lines(bi$l1, predvB, col = 1)
lines(bi$l1, predil, col = 2)
abline(linm, col = 3, lty = 2)
# Whoohooooo
# The linear regression and the vB are the same, I'm guessing because the vB we gave it was not a curve

parset(plots = c(1, 2))
plot(bi$l1, (bi$dl - predvB)); abline(h = 0, col = 1)
plot(bi$l1, (bi$dl - predil)); abline(h = 0, col = 1)
# Residual plots - the vB residuals are domed whereas the IL is more uniform.


### 5.3.4 A Closer Look at the Fabens Method ----

# What the Fabens transformation does is alter the residual structure and how parameters interact with each other.
# The Fabens model implies different things from the vB model

# The problem with growth models is that it assumes that variation around predicted growth (=residuals) is constant (=constant variance)
# But we can modify this by estimating a coefficient of variation that will modify variance over different lengths for example (coefficient sigma/variance mu)


### 5.3.5 Implementation of Non-Constant Variances ----

# Previously with a constant variance we used negNLL() as a criterion for model fit. Let's maybe change that
?negnormL # Is what we use instead of negNLL within the nlm() call. This will require "funksig" which is the function that is used to calculate sigma (the coefficient for variance)

sigfunk <- function(pars, predobs) return(tail(pars, 1)) # Choosing a sigma that will have no effect on variance
param <- c(170, 0.3, 4); label = c("Linf", "K", "sigma")
modelvb <- nlm(f = negnormL, p = param, funk = fabens, funksig = sigfunk, indat = bi, initL = "l1", delT = "dt")
outfit(modelvb, backtran = FALSE, title = "vB constant sigma")


sigfunk2 <- function(pars, predo) { # Linear with predicted length sigma x predDL, see negnormL, en
  sig <- tail(pars, 1) * predo
  pick <- which(sig <= 0) # No negative sigmas
  sig[pick] <- 0.01 # Possible negative predicted lengths
  return(sig)
}
param <- c(170, 0.3, 1); label = c("Linf", "K", "sigma")
modelvB2 <- nlm(f = negnormL, p = param, funk = fabens, funksig = sigfunk2, indat = bi, initL = "l1", delT = "dt")
outfit(modelvB2, backtran = FALSE, title = "vB inverse DeltaL, sigma < 1")
# Boom now the best fit estimates are different. Because we're giving it different variance structures.

parset(c(1,1))
predvB <- fabens(modelvB$estimate, bi)
predvB2 <- fabens(modelvB2$estimate, bi)
plot(bi$l1, bi$dl)
lines(bi$l1, predvB, col = 1)
lines(bi$l1, predvB2, col = 2)
# it is a little change between the two, a more complex residual structure make model fitting more sensitive to initial conditions
# You have to have good reason to complexify a model, but a changing variance is pretty easy to tell in data.



## 5.4 Objective Model Selection ----

# Another criteria you can base your model fit on is Akaike's Information Criterion


### 5.4.1 Akaike's Information Criterion ----

# If two models give equivalent results, always choose the simplest one.
# Akaike's Information Criterion (AIC) is a likelihood-based criteria for model selection
# The smaller the AIC, the better the model
# AIC basically penalises the model's fit based on how many parameters are being tested.

# From AIC came Schwartz' Bayesian Information Criterion (BIC)
# BIC will penalise complex models more than AIC (only if sample size > 8)


### 5.4.2 Likelihood Ratio Test ----

# Another comparison between different model fits v. complexity is Likelihood Ratio tests
# We want to determine if one of two models (same data, same residual structure, but different model structures or parameters) provide a sig. better  fit to available data
# If these models differ by 1.92 in that ratio test, then they're considered to be sig. different, and one will fit the data better

vb <- modelvb$minimum
il <- modelil$minimum
dof <- 1
round(likeratio(vb, il, dof), 8) # Here the LR test is 28, (wayy over 1.92) so one is much better than the other


### 5.4.3 Caveats on Likelihood Ratio Tests ----

# Always compare comparable models with likelihood ratio tests. duh



## 5.5 Remarks on Growth ----

# Growth curves estimated with existing data, and thus can be weak wherever data is sparse.
# Don't mistake the interpretation of the model's parameters for reality. The curve is only valid where there is data.



## 5.6 Maturity ----

### 5.6.1 Introduction

# Fisheries Management tries to have a target reference point (mature or spawning biomass level) as a desirable state for the stock
# In Aus, apparently the target is B48 (or 48% of the original biomass) as a proxy for MEY

# You also need a limit reference point, beyond which you consider your stock to be overfished

# To produce management advice you need:
# 1. An assessment of the current stock status (spawning biomass depletion)
# 2. Current estimate of unfished spawning biomass
# 3. A Harvest control rule that determines next season's fishing mortality, effort, or catch.

# L50 (length at which ~50% of individuals are mature) is a common measure of maturity
# Life histories will depend on species of course


### 5.6.2 Alternative Maturity Ogives ----

# For maturity we use an S shaped logistic curve that describe the % of mature individuals at/up to a certain length

data(tasab) # This is maturity data for Abalone
properties(tasab)
table(tasab$site, tasab$sex)

propm <- tapply(tasab$mature, tasab$length, mean)
lens <- as.numeric(names(propm))
plot1(lens, propm, type = "p", cex = 0.9, xlab = "Length mm", ylab = "Proportion mature")
# There are 2 sites, and each length is pooled and the % of mature is given a point. Depending on how many observations are in each length, the % may be skewed
# For ex. the single point at 0.5 maturity is just made of 2 observations.. Not super robust.

binglm <- function(x, digits = 6) {
  out <- summary(x)
  print(out$call)
  print(round(out$coefficient, digits))
  cat("\nNull Deviance ", out$null.deviance, "df", out$df.null, "\n")
  cat("Resid.Deviance ", out$deviance, "df", out$df.residual, "\n")
  cat("AIC = ", out$aic, "\n\n")
  return(invisible(out))
}

tasab$site <- as.factor(tasab$site)
smodel <- glm(mature ~ site + length, family = binomial, data = tasab)
outs <- binglm(smodel)

model <- glm(mature ~ length, family = binomial, data = tasab)
outm <- binglm(model)
cof <- outm$coefficients
cat("Lm50 = ", -cof[1, 1]/cof[2, 1], "\n")
cat("IQ = ", 2*log(3)/cof[2, 1], "\n")

# The AIC for the site&length model is higher than the IAC for the length model
# This indicates that site is not an informative variable to include

propm <- tapply(tasab$mature, tasab$length, mean)
lens <- as.numeric(names(propm))
pick <- which((lens > 79) & (lens < 146))
parset()
plot(lens[pick], propm[pick], xlab = "Length mm", ylab = "Proportion mature")
L <- seq(80, 145, 1)
pars <- coef(smodel)
lines(L, mature(pars[1], pars[3], L), col = 2)
lines(L, mature(pars[1] + pars[2], pars[3], L), col = 3)
lines(L, mature(coef(model)[1], coef(model)[2], L), col = 4)
abline(h = c(0.25, 0.5, 0.75), lty = 3, col = "gray")
# It's meant to givve 3 lines, 1 per site and one for both sites, showing a slight but insig. difference between sites


### 5.6.3 The Assumption of Symmetry ----

# The logistic curve assumes symmetry around the L50 point, which may not reflect reality
# But of course they've made asymetrical maturity curves so let's try those

L <- seq(50, 160, 1)
p <- c(a = 0.07, b = 0.2, c = 1, alpha = 100) # You can change the parameters and see how they alter the curve shape
asym <- srug(p = p, sizeage = L) # Shnute & Richards Unified Growth curve
L25 <- linter(bracket(0.25, asym, L))
L50 <- linter(bracket(0.5, asym, L))
L75 <- linter(bracket(0.75, asym, L))
parset()
plot(L, asym, type = "l")
abline(h = c(0.25, 0.5, 0.75))
abline(v = c(L25, L50, L75), lwd = c(1, 2, 1), col = c(1, 2, 1))
# Asymetry is visible by the unequal distance between L25-L50 and L75-L50



## 5.7 Recruitment ----

# 5.7.1 Introduction ----

# Growth and recruitment are the two main contributors of biomass addition in a stock (excluding immigration)
# Recruitment is highly *highly* variable. Whether it's dictated by the spawning stock or environmental variables or both

# here we largely ignore biology and try two of the most common recruitment models.


### 5.7.2 Properties of "Good" Stock Recruitment Relationships ----

# Ricker (1975) lists 4 properties of a stock recruitment relationship that he considered desirable:
# -curve should pass through the origin (if stock = 0 then recruitment = 0)
# -Recruitment should not fall to zero at high stock densities (declines are okay, but not zero)
# -Rate of recruitment (recruit-per-spawner) should decrease continuously with increases in parental stock (density-dependent effect) (this may not always hold)
# -Recruitment must exceed parental stock over some part of the range of possible parental stocks (only true for single-spawning species, but for multi-spawners, read this as "recruitment should be high enough to replace the annual natural mortality deaths")

# Hilborn & Walters added 2 properties to this:
# -The spawning stock curve should be continuous, no sharp changes over small changes of stock size (recruitment should vary smoothly with stock size)
# -The average stock recruitment relationship is constant over time (will likely fail where the ecosystem changes markedly, but in models you can use time-blocks of parameters)


### 5.7.3 Recruitment Overfishing ----

# 2 types of overfishing:
# Growth overfishing: fished before individuals can reach optimal size
# Recruitment overfishing: fished so hard the population can't recruit enough to replace those that die

# A Limit Reference Point is usually used to say 'okay beyond this the population can't replenish'.
# 20% of the original unfished biomass is the one used commonly


### 5.7.4 Beverton & Holt Recruitment ----

B <- 1:3000
bhb <- c(1000, 500, 250, 150, 50) # parameter b
parset()
plot(B, bh(c(1000, bhb[1]), B), ylim = c(0, 1050), xlab = "Spawning Biomass", ylab = "Recruitment")
for (i in 2:5) lines(B, bh(c(1000, bhb[i]), B), lwd = 2, col = i, lty = i)
# boom, Beverton Holt curve for different parameters b, and constant parameter a (1000)


### 5.7.5 Ricker Recruitment ----

# This one exhibits a recruitment decline at higher spawning biomass numbers, the idea is that competitivity or predatory effects (cannibalism of juv. by adults) reduces recruitment at high densities

B <- 1:20000
rickb <- c(0.0002, 0.0003, 0.0004)
parset()
plot(B, ricker(c(10, rickb[1]), B), xlab = "Spawning Biomass", ylab = "Recruitment")
for (i in 2:3) lines(B, ricker(c(10, rickb[i]), B), lwd = 3, col = i, lty = i)


### 5.7.6 Deriso's Generalised Model ----

# Deriso-Schnute proposes a generalised equation of which Beverton-Holt and Ricker are special cases
# Only use it in a simulation model, not in a fitted model.

deriso <- function(p, B) return(p[1] * B * (1 - p[2] * p[3] * B)^(1/p[3]))
B <- 1:10000
plot1(B, deriso(c(10,0.001, -1), B), xlab = "Spawning Biomass", ylab = "Recruitment")  # Beverton-Holt curve, but through Deriso-Schnute
lines(B, deriso(c(10,0.0004, 0.25), B), col = 2) # Deriso-Schnute
lines(B, deriso(c(10,0.0004, 1e-06), B), col = 3) # Ricker, but through Deriso-Schnute
lines(B, deriso(c(10,0.0004, 0.5), B), col = 1) # Odd line, unrealistic outcomes


### 5.7.7 Re-Parameterised Beverton-Holt Equation ----

# Changing things around in Beverton-Holt, for what reason I don't know.
# Using "steepness" of the initial curve as a parameter


### 5.7.8 Re-Parameterised Ricker Equation ----

# Same thing here, using steepness (again idk why)



## 5.8 Selectivity ----

### Introduction ----

# The fishing gear will determine which individuals of the available population will be vulnerable to it.
# Selectivity curves should really be viewed as selevivity/availability curves:
# Estimating selectivity using shallow water data v. deep water data will be different (juveniles may be more abundant in shallows, whereas fishing grounds may be in deeper waters)

# Choice of selectivity curve is important, each gear type will be described by a diff equation.
# Selectivity models in stock assessment can only be fitted if age/size-composition data is available from catches.
# Multinomial likelihoods would usually be used to fit it.

# When you generate predicted catch composition data, you multiply number/size-at-age data by predicted selectivity


### 5.8.2 Logistic Selection ----

# Same curve as for maturity -
# This assumes that vulnerability to gear gradually increases with size/age until 100% vulnerable

# 2 equations used, either logist() or mature()

ages <- seq(0, 50, 1); in50 <- 25
sel1 <- logist(in50, 12, ages)
sel2 <- mature(-3.650425, 0.146017, sizeage = ages)
sel3 <- mature(-6, 0.2, ages)
sel4 <- logist(22, 14, ages, knifeedge = T)
plot1(ages, sel1, xlab = "Age Years", ylab = "Selectivity")
lines(ages, sel2, col = 2) # Same L50 as the black line but diff gradients
lines(ages, sel3, col = 3) # A different parameter b
lines(ages, sel4, col = 4) # Knife edge selection, if age < x, not vulnerable, if age > x, all vulnerable


### 5.8.3 Dome-Shaped Selection ----

# Ascending, plateau and descending curve is also common
# Needs 3 parts, joined together - So the equation would have 5 parts, asc, plat, and desc, + 2 joining functions
# See book for details on equations

L <- seq(1, 30, 1)
p <- c(10,11, 16, 33, -5, -2)
plot1(L, domed(p, L), ylab = "Selectivity", xlab = "Age Yrs")
p1 <- c(8, 12, 16, 33, -5, -1)
lines(L, domed(p1, L), col = 2)
p2 <- c(9, 10, 16, 33, -5, -4)
lines(L, domed(p2, L), col = 3)
# Three selectivity curves, with changed initial and ending age-to-100-selectivity and final age class selectivity



## 5.9 Concluding Remarks for Static Models ----

# It's important to know individual elements (growth, maturity, recruitment) separately before getting to advanced models
# That way any interaction between components can be automatically accounted for.




# 6 On Uncertainty ----

## 6.1 Introduction ----

# Fitting a model = looking for parameter estimates that optimise the relationship between obs. and pred. values
# Exact parameter estimates are not the goal (given that the sample you have is only a part of the real population)
# What we want is repeatability of estimates
# We need to characterise the uncertainty of whatever model we're using. Some uncertainty influence the variability of collected data, other uncertainty influences the type of data available


### 6.1.1 Types of Uncertainty ----

# Uncertainty is often called "error" (as in residual error)
# Process uncertainty; natural random variation in demographic rates or other bio'l properties and processes
# Observation uncertainty: sampling error and measurement error- reflect the fact that samples are representations of a population but can't be perfectly representative
# Model uncertainty: relates to the capacity of the selected model structure to describe the dynamics of the system
  # (different structural models may provide different answers, thus uncertainty exists over which is the best representation of nature)
  # (selection of residual error structure is a special case of model uncertainty, can have implication for parameter estimates)
  # (estimation uncertainty = interactions or correlations between slightly diff. parameters sets, giving out identical log-likelihood)
# Implementation uncertainty: effects of management actions may differ from those intended. Poor definition of management options may lead to implementation uncertainty
  # (institutional uncertainty: bad management objectives leading to unworkable management)
  # (time-lag between making decisions and implementing them can lead to bigger variation)

# Model uncertainty can be both qualitative and quantitative:
# You can have two models, with same structure, and one being better when comparing residual structure
# Or if two models with diff. structures are to be compared, they are different descriptions of the system, you have to decide which one is the more sensical description

# Four different methods for characterising uncertainty around parameter estimates/model outputs
# 1. Bootstrapping: focus on uncertainty inherent in the data samples- functions by examining what would happen to parameter estimates if it were given slightly different samples
# 2. Asymptotic errors: use variance-covariance matrix between parameter estimates to describe uncertainty around parameter values
# 3. Likelihood profiles: on parameters of primary interest, constructed to obtain more specific distributions of each parameter
# 4. Bayesian marginal posteriors: characterise the uncertainty inherent in estimates of model parameters and outputs


## 6.1.2 The Example Model ----

data(abdat); logce <- log(abdat$cpue)
param <- log(c(0.42, 9400, 3400, 0.05))
label <- c("r", "K", "Binit", "sigma")
bestmod <- nlm(f = negLL, p = param, funk = simpspm, indat = abdat, logobs = logce)
outfit(bestmod, title = "SP-Model", parnames = label)

predce <- exp(simpspm(bestmod$estimate, abdat))
optresid <- abdat[, "cpue"]/predce
ymax <- getmax(c(predce, abdat$cpue))
plot1(abdat$year, (predce*optresid), maxy = ymax, ylab = "CPUE", xlab = "year")
points(abdat$year, abdat$cpue, col = 1); lines(abdat$year, predce, col = 2)

# "What is the plausible spread of the predicted CPUE around observed data ?"


## 6.2 Bootstrapping ----

# Bootstrapping = within your total sample, fit your model to multiple random subsections of the sample (can be used to fit standard errors, confidence intervals)


### 6.2.1 Empirical Probability Density Distributions ----

# Assumption is that for a sample from a population, the population's probability density distribution is the sample itself.
# == for n observations, each observation has equal likelihood of occurring (multiple observations may have the same value, but each obs. has the same likelihood)



## 6.3 A Simple Bootstrap Example ----

data(npf)
model <- lm(endeavour ~ tiger, data = npf)
plot1(npf$tiger, npf$endeavour, xlab = "Tiger Prawn (t)", type = "p", ylab = "Endeavour Prawn (t)")
abline(model, col = 1, lwd = 2)
correl <- sqrt(summary(model)$r.squared)
pval <- summary(model)$coefficients[2, 4]
label <- paste("Corr", round(correl, 5), "p ", round(pval, 8)); text(2700, 180, label, pos = 4)
# Prawn catch fishery (x = target spp, y = bycatch spp). quite noisy data
# Very sig. regression, but variability of data means it's hard to know with what confidence we can believe it.
# Let's try bootstrapping; taking 5000 subsamples of this data and calculate correlation coefficient.

set.seed(12321)
N <- 5000 # Number of bootstrap samples
result <- numeric(N)

for (i in 1:N) { # This is the bootstrapping: "reshuffle our 23 datapoints and tell me the correlation between target and bycatch spp"
  pick <- sample(1:23, 23, replace = T) # 23 because data has 23 years of data (datapoints)
  result[i] <- cor(npf$tiger[pick], npf$endeavour[pick])
}

rge <- range(result)
CI <- quants(result)
restrim <- result[result > 0]
parset(cex = 1)
bins <- seq(trunc(range(restrim)[1]*10)/10, 1, 0.01)
outh <- hist(restrim, breaks = bins, main = "", xlab = "Correlation")
abline(v = c(correl, mean(result)), col = c(4, 3), lty = c(1, 2))
abline(v = CI[c(2, 4)], col = 4, lwd = 2)
text(0.48, 400, makelabel("range", rge, sep = " ", sigdig = 4), pos = 4)
label <- makelabel("90% CI", CI[c(2, 4)], sep = " ", sigdig = 4); text(0.48, 300, label, pos = 4)
# Green dashed line is the correlation we observed with the total data, within blue lines (=if our data is in it, we are 95% confident that it's not a fluke)
# Bars are the correlation observations for the 5000 bootstrap samples we tested

# We can be confident that thhe high correlation in the original data is a fair representation of the relationship



## 6.4 Bootstrapping Time-Series Data ----

# When you're looking specifically at time-series, where values at time t depending on the value at t-1, bootstrapping is not valid or sensible (because we'd be subsampling the time-series)
# To bootstrap time-series data, you don't subsample the original data, but we get an optimum model bit to the OG data first, then bootstrap individual residuals at each point

data(abdat); logce <- log(abdat$cpue)
param <- log(c(r = 0.42, K = 9400, Binit = 3400, sigma =0.05))
bestmod <- nlm(f = negLL, p = param, funk = simpspm, indat = abdat, logobs = logce) # Finding best parameters using negLL as a criteria of fit
optpar <- bestmod$estimate # Optimal parameters
predce <- exp(simpspm(optpar, abdat)) # Predicted values using optimal parameters
optres <- abdat[, "cpue"]/predce # Optimal log-normal residuals
optmsy <- exp(optpar[1])*exp(optpar[2])/4
sampn <- length(optres) # number of residuals and number of years

# Let's bootstrap this time-series:
start <- Sys.time()
bootfish <- as.matrix(abdat)
N <- 1000; years <- abdat[, "year"]
columns <- c("r", "K", "Binit", "sigma")
results <- matrix(0, nrow = N, ncol = sampn, dimnames = list(1:N, years))
bootcpue <- matrix(0, nrow = N, ncol = sampn, dimnames = list(1:N, years))
parboot <- matrix(0, nrow = N, ncol = 4, dimnames = list(1:N, columns))

for (i in 1:N) { # fitting the model,
  bootcpue[i,] <- predce * sample(optres, sampn, replace = T) # multiply predicted values by the bootstrapped residuals (because we're suing multiplicative log-normal residuals - idk??)
  bootfish[, "cpue"] <- bootcpue[i,] # calculate the bootstrapped CPUE
  bootmod <- nlm(f = negLL, p = optpar, funk = simpspm, indat = bootfish, logobs = log(bootfish[, "cpue"])) # Iteratively fit the model using all the bootstrapped residuals we created above
  parboot[i,] <- exp(bootmod$estimate) # 1000 sets of optimal parameters, for each bootstrapped residuals
  results[i,] <- exp(simpspm(bootmod$estimate, abdat)) # Bootstrapped predicted values ?
}
cat("total time = ", Sys.time()-start, "seconds \n") # This whole thing took about 8 seconds - keep that in mind for larger/more complex models

plot1(abdat[, "year"], abdat[, "cpue"], type = "n", xlab = "year", ylab = "cpue") # Empty plot
for (i in 1:N)
  lines(abdat[,"year"], results[i, ], lwd = 1, col = "gray") # This adds a line for every bootstrap model we made!! Holy crap
points(abdat[, "year"], abdat[, "cpue"], pch = 16, col = 1) # Observed cpue
lines(abdat[, "year"], predce, lwd = 2, col = 1) # Black line for the optimal model after bootstrappin'

# This is super cool but the star of the show is the optimal parameters, so we need to plot the bootstrapped values of r, K, Binit (the parameters of SPM function simpspm)
# And also bootstrapped MSY calculations - why not.

dohist <- function(invect, nmvar, bins = 30, bootres, avpar) {
  hist(invect[, nmvar], breaks = bins, main = "", xlab = nmvar, col = 0)
  abline(v = c(exp(avpar), bootres[pick, nmvar]), lwd = c(3, 2, 3, 2), col = c(3, 4, 4, 4))
}
msy <- parboot[, "r"] * parboot[, "K"]/4 # calculating MSY with all the bootstrapped values
msyB <- quants(msy) # CI for MSY

parset(plots = c(2, 2), cex = 0.9)
bootres <- apply(parboot, 2, quants); pick <- c(2, 3, 4)
dohist(parboot, nmvar = "r", bootres = bootres, avpar = optpar[1])
dohist(parboot, nmvar = "K", bootres = bootres, avpar = optpar[2])
dohist(parboot, nmvar = "Binit", bootres = bootres, avpar = optpar[3])
hist(msy, breaks = 30, main = "", xlab = "MSY", col = 0)
abline(v = c(optmsy, msyB[pick]), lwd = c(3, 2, 3, 2), col = c(3, 4, 4, 4))
# boom we have the frequency of parameter estimates from the bootstrapping, with 90% CI as blue lines.
# MSY is also calculated with all the bootstrapped values, giving you an idea of what the most likely MSY is.
# COooooooOOl!!


### 6.4.1 Parameter Correlation ----

# Let's look at whether there are correlations between parameters and model outputs

parboot1 <- cbind(parboot, msy)
pairs(parboot1, pch = 16, col = rgb(red = 1, green = 0, blue = 0, alpha = 1/20))
# There's definitely a correlation between r, K and Binit, so in models you might get the same MSY for very different values of r and K, for ex.



## 6.5 Asymptotic Errors ----

# Confidence intervals (90% or 95% CI) is defined as :
# sample mean +/- (t distribution of degrees of freedom)*(sample standard deviation/sqrt(sample size))
# Degrees of freedom: maximum number of logically independent values. Always n - 1 for some reason

# What we want to know is how sure are we of the parameter estimates given by the model?
# In multi-parameter models, you need standard errors for all parameters, which can be done by making a variance-covariance matrix for the parameters in the vicinity of optimum parameter set (here, Hessian matrix)
# In this matrix, you assume that maximum likelihood (the fit criteria) is multivariate Normally distributed
# From the matrix, you can get standard errors and then confidence intervals
# Max. likelihood --> variance-covariance matrix (assumed Normally Distributed) --> parameter standard error --> parameter confidence interval

data(abdat)
param <- log(c(r = 0.42, K = 9400, Binit = 3400, sigma = 0.05))

# Just below is the same optimisation process we've done so far, with the added "hessian = TRUE" which creates a Hessian matrix for our parameters on top of everything else
bestmod <- nlm(f = negLL, p = param, funk = simpspm, indat = abdat, logobs = log(abdat[, "cpue"]), hessian = TRUE)

outfit(bestmod, backtran = TRUE)
bestmod # Here you see the added matrix on the output

vcov <- solve(bestmod$hessian) # this inverts the matrix (^-1)
sterr <- sqrt(diag(vcov)) # extracts the diagonal of the matrix, which is used to get the standard error of each parameter
optpar <- bestmod$estimate # take our optimised model estimates
U95 <- optpar + qt(0.975, 20)*sterr # upper confidence interval
L95 <- optpar - qt(0.975, 20)*sterr # lower confidence interval

# Boom, you got your upper and lower CI for every parameter in the modellllllll wooohoooo
cat("\n               r      K       Binit     sigma \n")
cat("Upper 95%", round(exp(U95), 5), "\n") # Backtransforming
cat("Optimum  ", round(exp(optpar), 5), "\n")
cat("Lower 95%", round(exp(L95), 5), "\n") # Backtransforming


### 6.5.1 Uncertainty about the Model Outputs ----

# Above we have CI for model parameters.
# To get confidence intervals around model outputs.
# Same things as above, you assume that log-likelihood around the optimum solution is symmetrical (multivariate Normal)
# But from hereon, you have to generate random parameter vectors from the estimated Normal distribution (mean of optimum parameter)
# The random parameter vectors are like bootstrap values that you generate CI from
# Using multivariate Normal, parameter-to-parameter correlations are accounted for (not sure how)


### 6.5.2 Sampling from a Multivariate Normal Distribution ----

library(mvtnorm)
N <- 1000 # number of multivariate normal parameter vectors
years <- abdat[, "year"]; sampn <- length(years) # 24 years
# Then create a matrix of 24 columns, 1000 rows
mvncpue <- matrix(0, nrow = N, ncol = sampn, dimnames = list(1:N, years))
columns <- c("r", "K", "Binit", "sigma")

# Then fill another matrix with 1000 random vectors of our 4 parameters from the previous section, picked from a multivariate random distribution
mvnpar <- matrix(exp(rmvnorm(N, mean = optpar, sigma = vcov)), nrow = N, ncol = 4, dimnames = list(1:N, columns))

# Calculate 1000 cpue trajectories from our 1000 vectors of parameters
for (i in 1:N) mvncpue[i,] <- exp(simpspm(log(mvnpar[i,]), abdat))
msy <- mvnpar[, "r"]*mvnpar[, "K"]/4

plot1(abdat[, "year"], abdat[,"cpue"], type = "p", xlab = "Year", ylab = "cpue", cex = 0.9)

for (i in 1:N) lines(abdat[, "year"], mvncpue[i,], col = "gray", lwd = 1)
points(abdat[, "year"], abdat[, "cpue"], pch = 16, cex = 1) # Original data
lines(abdat[, "year"], exp(simpspm(optpar, abdat)), lwd = 2, col = 1) # The optimal parameter set
# HHHOOOOOLLLLYYYYY MOOOOOOOLLLLYYYY
# Okay so we have the 1000 trajectories from the 1000 randomly picked sets of parameters (defined by multivariate Normal distribution around the optimum parameters and variance-covariance matrix)

# Correlations between 1000 parameters
pairs(cbind(mvnpar, msy), pch = 16, col = rgb(red = 1, 0, 0, alpha = 0.1))


mvnres <- apply(mvnpar, 2, quants)
pick <- c(2, 3, 4)
meanmsy <- mean(msy)
msymvn <- quants(msy)

plothist <- function(x, optp, label, resmvn) {
  hist(x, breaks = 30, main = "", xlab = label, col = 0)
  abline(v = c(exp(optp), resmvn), lwd = c(3, 2, 3, 2), col = c(3, 4, 4, 4))
}

par(mfrow = c(2, 2), mai = c(0.45, 0.45, 0.05, 0.05), oma = c(0, 0, 0, 0))
par(cex = 0.85, mgp = c(1.35, 0.35, 0), font.axis = 7, font = 7, font.lab = 7)
plothist(mvnpar[, "r"], optpar[1], "r", mvnres[pick, "r"])
plothist(mvnpar[, "K"], optpar[2], "K", mvnres[pick, "K"])
plothist(mvnpar[, "Binit"], optpar[3], "Binit", mvnres[pick, "Binit"])
plothist(msy, meanmsy, "MSY", msymvn[pick])
# These are plots of 1000 parameter estimates, and the derived MSY. Green line is the arithmetic mean, central blue line is median, outer lines are 90% CI around the median



## 6.6 Likelihood Profiles ----

# To understand the relative contribution of each model parameter on the model fit (negLL) we make likelihood profiles
# What we're doing is fixing one of the parameters a bit far from its optimal value, and refitting the model, finding a new optimum, and seeing how much the negLL is worsened

# Below is a completely normal process we've done before for fitting a model
data(abdat); logce <- log(abdat$cpue)
param <- log(c(r = 0.43, K = 9400, Binit = 3400, sigma = 0.05))
optmod <- nlm(f = negLL, p = param, funk = simpspm, indat = abdat, logobs = logce)
outfit(optmod, parnames = c("r", "K", "Binit", "sigma"))

# Now to fix one parameter and allow all others to be tested, we use negLLP() instead of negLL as a model fit criteria
# If you ignore the initpar and notfixed arguments, you should get the same outputs as negLL
negLLP <- function(pars, funk, indat, logobs, initpar = pars, notfixed = c(1:length(pars)), ...) {
  usepar <- initpar # initial parameters
  usepar[notfixed] <- pars[notfixed] # parameters that are not fixed
  npar <- length(usepar)
  logpred <- funk(usepar, indat, ...)
  pick <- which(is.na(logobs))
  if (length(pick) > 0) {
    LL <- -sum(dnorm(logobs[-pick], logpred[-pick], exp(pars[npar]),
                     log = T))
  } else {
    LL <- -sum(dnorm(logobs, logpred, exp(pars[npar])))
  }
  return(LL)
}

# Let's test is negLLP gives the same outputs as negLL if we ignore initpar and notfixed
param <- log(c(r = 0.42, K = 9400, Binit = 3400, sigma = 0.05))
bestmod <- nlm(f = negLLP, p = param, funk = simpspm, indat = abdat, logobs = logce) # using negLLP but not initpar or notfixed
outfit(bestmod, parnames = c("r", "K", "Binit", "sigma")) # we're all good

# Okay now for the likelihood profile

rval <- seq(0.325, 0.45, 0.001) # test sequence
ntrial <- length(rval)
columns <- c("r", "K", "Binit", "sigma", "-veLL")
result <- matrix(0, nrow = ntrial, ncol = length(columns), dimnames = list(rval, columns)) # empty matrix to fill with the for-loop
bestest <- c(r = 0.32, K = 11000, Binit = 4000, sigma = 0.05)

for (i in 1:ntrial) { # for all test sequence,
  param <- log(c(rval[i], bestest[2:4])) # parameter r varies (rval sequence) and K, Binit and sigma are fixed on bestest values
  parinit <- param
  bestmodP <- nlm(f = negLLP, p = param, funk = simpspm, initpar = parinit, indat = abdat, logobs = log(abdat$cpue), notfixed = c(2:4), typsize = magnitude(param), iterlim = 1000) # optimise the model for K, Binit and sigma, using the fixed sequence of r
  bestest <- exp(bestmodP$estimate)
  result[i,] <- c(bestest, bestmodP$minimum) # store each result
}
result # you can see the negLL getting bigger and bigger (worse and worse model fit) as r goes further from its optimal value (0.32 i think)
minLL <- min(result[, "-veLL"])


### 6.6.1 Likelihood Ratio-Based Confidence Intervals ----

# We can plot the likelihood profile we just made:
plotprofile(result, var = "r", lwd = 2)

# Using almost identical code to above, testing the likelihood profile for K

Kval <- seq(7200, 12000, 10) # test sequence
ntrial <- length(Kval)
columns <- c("r", "K", "Binit", "sigma", "-veLL")
resultK <- matrix(0, nrow = ntrial, ncol = length(columns), dimnames = list(Kval, columns)) # empty matrix to fill with the for-loop
bestest <- c(r = 0.45, K = 7500, Binit = 2800, sigma = 0.05)

for (i in 1:ntrial) { # for all test sequence,
  param <- log(c(bestest[1], Kval[i], bestest[3:4])) # parameter K varies (Kval sequence) and r, Binit and sigma are fixed on bestest values
  parinit <- param
  bestmodP <- nlm(f = negLLP, p = param, funk = simpspm, initpar = parinit, indat = abdat, logobs = log(abdat$cpue), notfixed = c(1, 3, 4), typsize = magnitude(param), iterlim = 1000) # optimise the model for K, Binit and sigma, using the fixed sequence of r
  bestest <- exp(bestmodP$estimate)
  resultK[i,] <- c(bestest, bestmodP$minimum) # store each result
}
resultK # you can see the negLL getting bigger and bigger (worse and worse model fit) as r goes further from its optimal value (0.32 i think)
minLL <- min(resultK[, "-veLL"])

plotprofile(resultK, var = "K", lwd = 2) # a bit more asymetrical than r


### 6.6.2 -ve Log-Likelihoods or Likelihoods ----

# Neg log likelihoods aren't super intuitive, so let's backtransform them
likes <- exp(-resultK[, "-veLL"])/sum(exp(-resultK[, "-veLL"]), na.rm = TRUE)
resK <- cbind(resultK, likes, cumlike = cumsum(likes))
plot(likes)


### 6.6.3 Percentile Likelihood Profiles for Model Outputs ----

# So far it's all been likelihood profiles for parameters, but when you need MSY (kind of an indirect output from the model) likelihood profiles for ex. whaddaya do ?
# You can make these likelihood profiles by adding a penalty term to the negative log-likelihood that constrains likelihood to the optimal target
# See book for equations
# Essentially the penalty is in the form of a new model fit criteria (negLLO) where a weighting factor is involved (the bigger, the more penalty)
# We'll use negLLO in the likelihood profile to find the optimal parameters that give the best output (MSY), with lowest penalty for not being exactly the optimal parameters

negLLO <- function(pars, funk, indat, logobs, wght, optvar, varval) {
  logpred <- funk(pars, indat)
  LL <- -sum(dnorm(logobs, logpred, exp(tail(pars, 1)), log = T)) + wght * ((varval - optvar)/optvar)^2
  return(LL)
}

# We're testing the effect of different MSY values on the neg log likelihood
msyP <- seq(740, 1020, 2.5);
optmsy <- exp(optmod$estimate[1])*exp(optmod$estimate[2])/4
ntrial <- length(msyP)
wait <- 400
columns <- c("r", "K", "Binit", "sigma", "-veLL", "MSY", "pen", "TrialMSY")
resultO <- matrix(0, nrow = ntrial, ncol = length(columns), dimnames = list(msyP, columns))
bestest <- c(r = 0.47, K = 7300, Binit = 2700, sigma = 0.05)

for (i in 1:ntrial) {
  param <- log(bestest)
  bestmodO <- nlm(f = negLLO, p = param, funk = simpspm, indat = abdat, logobs = log(abdat$cpue), wght = wait, optvar = optmsy, varval = msyP[i], iterlim = 1000)
  bestest <- exp(bestmodO$estimate)
  ans <- c(bestest, bestmodO$minimum, bestest[1]*bestest[2]/4, wait*((msyP[i] - optmsy)/optmsy)^2, msyP[i])
  resultO[i,] <- ans
}
resultO
minLLO <- min(resultO[, "-veLL"])
plotprofile(resultO, var = "TrialMSY")

#Adjusting the weighting has to be done manually
# Not sure I completely understand this.


## 6.7 Bayesian Posterior Distributions ----

# If your likelihood profile is very steep around your optimum, you can have fairly high confidence in your optimum
# If it's pretty flat around the optimum, uncertainty is a bit higher.
# When you have loads of parameters' likelihood profiles to compare, it gets a bit more fiddly
# There are many ways to do this, but we'll use Markov Chain Monte Carlo (MCMC), specifically the Gibbs-within-Hastings approach (more flexible apparently)

# Markov Chain describes a process, where each step is determined probabilistically from the previous step. The final state of the chain is a description of a random distribution.
# The aim here is to produce a Markov Chain whose final equilibrium state (stationary distribution) gives a description of the target (posterior distribution of Bayesian stats)

# You give the Markov Chain parameter values, your observed data, and the model you're using, and these make a likelihood space.
# The MCMC process will step through the parameter space (with rules on likelihoods of each new candidate parameter set) to determine which steps become part of the Markov Chain.
# Each step is a question of "which parameter vectors will be accepted and which rejected?"
# Each iteration (step) produces a new set of parameters to test, and in Gibbs-within-Hastings, one parameter at a time
# Each iteration (=set of parameter tested) combined with observed data and model is given a likelihood
# The set of parameters is accepted or rejected based on how much the likelihood changes from the previous step

# This process is likely to result in auto-correlation between sequential parameter sets, which can bias the conclusions you make about variation across parameter space
#(if you get the same path every time, you'll think this is the only good parameter set, which may not be true)
# To solve this, you can thin out the chain so that the final chain contains only every n'th point. (worked example later)


### 6.7.1 Generating the Markov Chain ----

# To generate a Markov Chain you need to define :
# the likelihood of an initial set of parameters given a set of data and the Bayesian prior probability of the parameter set
# You then generate a new canditate parameter set by randomly incrementing at least one parameter, which will alter the implied likelihood
# If the ratio between previous and candidate likelihoods >1, the it's accepted by the Markov Chain


### 6.7.2 The Starting Point ----

# Ideally you would start the chain by selecting a vector of parameters close to but not identical to the optimum solution


### 6.7.3 The Burn-In Period ----

# When the chain starts, it'll stochastically meander before generally going to vectors of parameters with higher likelihoods
# This initial bit of vectors are known as the 'burn-in period' and it's recommended to remove them.
# Removing the first couple hundred is good enough in our case


### 6.7.4 Convergence to the Stationary Distribution ----

# Markov Chains will generally converge to a stable solution
# There have been questions of "how many chains" and "how long should the chains be" to arrive at this stable solution
# We'll look at empirical methods for determining evidence for convergence


### 6.7.5 The Jumping Distribution ----

# When testing a new candidate parameter vector, you need to generate that new set according to a "jumping distribution" (essentially how big a step between current and candidate parameters)
# The Normal random deviate distribution is common, scaled by alpha
# If you're testing parameter sets that are too different, your success rate will be too low and require a lot of iterations to converge,
# But if they're too similar, success rate will be too high and never converge
# Scaling takes trial and error


### 6.7.6 Application of MCMC to the Example ----

# We'll apply a Markov Chain to abdat

data(abdat)
plotspmdat(abdat)


### 6.7.7 Markov Chain Monte Carlo ----

# Pre-requisites for MCMC:
# - functions to calculate the negLL and prior probability or each candidate parameter set
# - the parameter set we'll start the MCMC with, and the burn-in size
# - thinning rate (how often in the chain we accept a result)
# - how many independent chains we're generating, and how long they each are
# - how to select weightings (scaling)


### 6.7.8 A First Example of an MCMC ----

data(abdat); logce <- log(abdat$cpue)
fish <- as.matrix(abdat)
begin <- Sys.time() # Allows you to calculate how long it all takes to run
chains <- 1 # 1 chain per run, usually more
burnin <- 0 # for the first 3 chains, no burn-in
N <- 100 # Number of steps in the chain
step <- 4 # thinning rate of 1 x 4 parameters = 4, so no thinning (???)
priorcalc <- calcprior # prior probability function
scales <- c(0.065, 0.055, 0.065, 0.425) # found by trial and error
set.seed(128900) # to get same results as the book, omit irl

inpar <- log(c(r = 0.4, K = 11000, Binit = 3600, sigma = 0.05))
result1 <- do_MCMC(chains, burnin, N, step, inpar, negLL, calcpred = simpspm, calcdat = fish, obsdat = logce, priorcalc, scales)

inpar <- log(c(r = 0.35, K = 8500, Binit = 3400, sigma = 0.05))
result2 <- do_MCMC(chains, burnin, N, step, inpar, negLL, calcpred = simpspm, calcdat = fish, obsdat = logce, priorcalc, scales)

inpar <- log(c(r = 0.55, K = 9500, Binit = 3200, sigma = 0.05))
result3 <- do_MCMC(chains, burnin, N, step, inpar, negLL, calcpred = simpspm, calcdat = fish, obsdat = logce, priorcalc, scales)

burnin <- 50
step <- 16 # thinning rate of 4 x 4 parameters = 16
N <- 10000 # 16 x 10000 = 160 000 steps + 50 burnin

inpar <- log(c(r = 0.4, K = 9400, Binit = 3400, sigma = 0.05))
result4 <- do_MCMC(chains, burnin, N, step, inpar, negLL, calcpred = simpspm, calcdat = fish, obsdat = logce, priorcalc, scales)

post1 <- result1[[1]][[1]]
post2 <- result2[[1]][[1]]
post3 <- result3[[1]][[1]]
postY <- result4[[1]][[1]]
cat("time   =", Sys.time() -begin, "\n")
cat("Accept =", result4[[2]], "\n")


parset(cex = 0.85)
P <- 75
plot(postY[, "K"], postY[, "r"], type = "p", cex = 0.2, xlim = c(7000, 13000), ylim = c(0.28, 0.47), col = 8, xlab = "K", ylab = "r", panel.first = grid())
lines(post2[1:P, "K"], post2[1:P, "r"], lwd = 1, col = 1)
points(post2[1:P, "K"], post2[1:P, "r"], pch = 15, cex = 1)

lines(post1[1:P, "K"], post1[1:P, "r"], lwd = 1, col = 1)
points(post1[1:P, "K"], post1[1:P, "r"], pch = 1, cex = 1)

lines(post3[1:P, "K"], post3[1:P, "r"], lwd = 1, col = 1)
points(post3[1:P, "K"], post3[1:P, "r"], pch = 2, cex = 1)

# This is a plot of the path that the chains took for the 3 chains we tried.
# They start at different origins (that we set) and have no burn-in period.
# The gray points are from a fourth chain with a 50 point burn-in and a thinning rate of 4, giving an idea of where chains should converge (the stationary distribution)

# Verify parameter correlation
posterior <- result4[[1]][[1]]
msy <- posterior[, 1]*posterior[, 2]/4
pairs(cbind(posterior[, 1:4], msy), pch = 16, col = rgb(1, 0, 0, 1/50), font = 7)

posterior <- result4[[1]][[1]]
par(mfrow = c(4, 2), mai = c(0.4, 0.4, 0.05, 0.05), oma = c(0, 0, 0, 0))
par(cex = 0.8, mgp = c(1.35, 0.35, 0), font.axis = 7, font = 7, font.lab = 7)
label <- colnames(posterior)
N <- dim(posterior)[1]
for (i in 1:4) {
 ymax <- getmax(posterior[, i]); ymin <- getmin(posterior[, i])
 plot(1:N, posterior[, i], type = "l", lwd = 1, ylim = c(ymin, ymax), panel.first = grid(), ylab = label[i], xlab = 'Step')
 plot(density(posterior[, i]), lwd = 2, col = 2, panel.first = grid(), main = "")
} # These are "trace" plots, telling you that for the 4th chain we generated, parameters went up and down as the chain went along (left plots), and their values were mostly around __ (right plots)
# Ideally you would have trace plots exactly like what sigma is doing (left plot is a hairy caterpillar, literally), but what we have for r, K and Binit is a lor of movement together == autocorrelation

# To test autocorrelation we use acf() - not sure what lag is but whatevs

posterior <- result4[[1]][[1]]
label <- colnames(posterior)[1:4]
parset(plots = c(2, 3), cex = 0.85)
for (i in 1:4) auto <- acf(posterior[, i], type = "correlation", lwd = 2, plot = T, ylab = label[i], lag.max = 20)
# This shows that r, K and Binit have high correlation, but not sigma

# If we run MCMC with higher thinning rate, this should reduce correlation

begin <- gettime()
scales <- c(0.06, 0.05, 0.06, 0.4)
inpar <- log(c(r = 0.4, K = 9400, Binit = 3400, sigma = 0.05))
result <- do_MCMC(chains = 1, burnin = 100, N = 1000, thinstep = 512, inpar, negLL, calcpred = simpspm, calcdat = fish, obsdat = logce, priorcalc, scales, schaefer = T)
posterior <- result[[1]][[1]]
label <- colnames(posterior)[1:4]
parset(plots = c(2, 2), cex = 0.85)
for (i in 1:4) auto <- acf(posterior[, i], type = "correlation", lwd = 2, plot = T, ylab = label[i], lag.max = 20)
# Woohoo correlation much reduced!!! from thinning rate increasing from 4 (16) to 128 (512)
# But it does take longer to run


### 6.7.9 Marginal Distributions ----

dohist <- function(x, xlab) {
  return(hist(x, main = "", breaks = 50, col = 0, xlab = xlab, ylab = "", panel.first = grid()))
}

param <- log(c(r = 0.42, K = 9400, Binit = 3400, sigma = 0.05))
bestmod <- nlm(f = negLL, p = param, funk = simpspm, indat = abdat, logobs = log(abdat$cpue))
optval <- exp(bestmod$estimate)
posterior <- result[[1]][[1]]
par(mfrow = c(5, 1), mai = c(0.4, 0.3, 0.025, 0.05), oma = c(0, 1, 0, 0))

np <- length(param)
for (i in 1:np) {
  outH <- dohist(posterior[, i], xlab = colnames(posterior)[i])
  abline(v = optval[i], lwd = 3, col = 4)
  tmp <- density(posterior[, i])
  scaler <- sum(outH$counts)*(outH$mids[2]-outH$mids[1])
  tmp$y <- tmp$y * scaler
  lines(tmp, lwd = 2, col = 2)
}
msy <- posterior[, "r"]*posterior[, "K"]/4
mout <- dohist(msy, xlab = "MSY")
tmp <- density(msy)
tmp$y <- tmp$y * (sum(mout$counts)*(mout$mids[2]-mout$mids[1]))
lines(tmp, lwd = 2, col = 2)
abline(v = (optval[1]*optval[2]/4))
# These plots tell us which values of each parameter are the most accepted in the chain



## 6.8 The Use of Rcpp ----

# Running steps of the MCMC can take a long time to run (depending on the computer)
# You can check what function is taking the longest to run with Rprof()
library(MQMF)

# (below takes forever to run, essentially the output is a table with the top time-intensive functions of the line of code. funk is the top one)
data(abdat); logce <- log(abdat$cpue); fish <- as.matrix(abdat)
param <- log(c(r = 0.39, K = 9200, Binit = 3400, sigma = 0.05))
Rprof(append = T)
result <- do_MCMC(chains = 1, burnin = 100, N = 20000, thinstep = 16, inpar = param, infunk = negLL1, calcpred = simpspm, calcdat = fish, obsdat = logce, priorcalc = calcprior, scales = c(0.07, 0.06, 0.07, 0.45))
Rprof(NULL)
outprof <- summaryRprof()

# To speed up the execution of MCMC (and many other things) we can combine R code with C++ (another programming language)
# I've downloaded a C++ compiler and the Rcpp package to use for this
# For compiler: CRAN repository website > download R for Windows > Rtools > download


### 6.8.1 Addressing Vectors and Matrices ----

# In a vector in R, you would address elements in the vector as index 1, 2, 3... for the first, second and third element of the vector.
# (ex. in c(r, K, Binit, sigma), r is indexed as 1, K as 2 etc.)
# In C++ they're indexed starting with 0, so r would be 0, K as 1, Binit as 2 and sigma as 3.
# Be very careful when you use R and C++, these mistakes can happen


### 6.8.2 Replacement for simpspm() ----

library(Rcpp)
# Essentially what we cant is a function in C++ that does exactly the same things as simpspm but quicker
# Below we're giving this long character string (in '...') with the C++ function
cppFunction('NumericVector simpspmC(NumericVector pars, NumericMatrix indat, LogicalVector schaefer) {
            int nyrs = indat.nrow();
            NumericVector predce(nyrs);
            NumericVector biom(nyrs+1);
            double Bt, qval;
            double sumq = 0.0;
            double p = 0.00000001;
            if (schaefer(0) == TRUE) {
              p = 1.0;
            }
            NumericVector ep = exp(pars);
            biom[0] = ep[2];
            for (int i = 0; i < nyrs; i++) {
              Bt = biom[i];
              biom[(i+1)] = Bt+(ep[0]/p)*Bt*(1-pow((Bt/ep[1]),p))-indat(i,1);
            if (biom[(i+1)] < 40.0) biom[(i+1)] = 40.0;
            sumq += log(indat(i,2)/biom[i]);
            }
            qval = exp(sumq/nyrs);
            for (int i = 0; i < nyrs; i++) {
              predce[i] = log(biom[i] * qval);
            }
            return predce;
}')

# Using microbenchmark we can compare the speed difference between the two functions, simpspm and simpspmC
library(microbenchmark)
data(abdat)
fishC <- as.matrix(abdat) # simpspmC needs a matrix to work with - faster than a dataframe

# Parameters and functions
inpar <- log(c(r = 0.389, K = 9200, Binit = 3300, sigma = 0.05))
spmR <- exp(simpspm(inpar, fishC))
spmC <- exp(simpspmC(inpar, fishC, schaefer = TRUE))

rows <- 1:length(spmR)
columns <- c("spmR", "spmC")

comp <- matrix(0, nrow = length(spmR), ncol = length(c("spmR", "spmC")), dimnames = list(rows, columns)) # empty matrix to fill with the for-loop
comp[, 1] <- spmR; comp[, 2] <- spmC # The two functions give the exact same outputs. yay!!

# Speed comparison
out <- microbenchmark(
  simpspm(inpar, fishC, schaefer = TRUE),
  simpspmC(inpar, fishC, schaefer = TRUE),
  times = 1000
);
out2 <- summary(out)[, 2:8]
out2 <- rbind(out2, out2[2,]/out2[1,])
rownames(out2) <- c("simpspm", "simpspmC", "TimeRatio")
out2 # The C++ function simpspmC is wayyy facter

# Below is another way to see how much faster C++ is
set.seed(167423)
beginR <- gettime()
setscale <- c(0.07, 0.06, 0.07, 0.45)
reps <- 2000
param <- log(c(R = 0.39, K = 9200, Binit = 3400, sigma = 0.05))
resultR <- do_MCMC(chains = 1, burnin = 100, N = reps, thinstep = 128, inpar = param, infunk = negLL1, calcpred = simpspm, calcdat = fishC, obsdat = log(abdat$cpue), schaefer = TRUE, priorcalc = calcprior, scales = setscale)
timeR <- gettime() - beginR
cat("time = ", timeR, "\n")
cat("acceptance rate = ", resultR$arate, " \n")
postR <- resultR[[1]][[1]]

set.seed(167423)
beginC <- gettime()
param <- log(c(R = 0.39, K = 9200, Binit = 3400, sigma = 0.05))
resultC <- do_MCMC(chains = 1, burnin = 100, N = reps, thinstep = 128, inpar = param, infunk = negLL1, calcpred = simpspmC, calcdat = fishC, obsdat = log(abdat$cpue), schaefer = TRUE, priorcalc = calcprior, scales = setscale)
timeC <- gettime() - beginC
cat("time = ", timeC, "\n")
cat("acceptance rate = ", resultC$arate, " \n")
postC <- resultC[[1]][[1]]
cat("Time Ratio = ", timeC/timeR) # Boom, much faster

# Because we set the same seed for both MCMC chains (R and C++ functions) we get the exact same density of parameters accepted (see below).
par(mfrow = c(1, 1), mai = c(0.45, 0.45, 0.05, 0.05), oma = c(0, 0, 0, 0))
maxy <- getmax(c(density(postR[, "K"])$y, density(postC[, "K"])$y))
plot(density(postR[, "K"]), lwd = 2, col = 1, xlab = "K", ylab = "Density", main ="", ylim = c(0, maxy), panel.first = grid())
lines(density(postC[, "K"]), lwd = 3, col = 3, lty = 2)


### 6.8.3 Multiple Independent Chains ----

# Ideally you'd have more than one chain. But it's a tradeoff between how much time you have and having many chains
# 3 minimum, the more the better, it provides proof of convergence.
# We'll run some example simple models, with 3 chains

setscale <- c(0.07, 0.06, 0.07, 0.45)
set.seed(9393074) # for same results as the book, omit irl
reps <- 10000
beginC <- gettime()
resultC <- do_MCMC(chains = 3, burnin = 100, N = reps, thinstep = 256, inpar = param, infunk = negLL1, calcpred = simpspmC, calcdat = fishC, obsdat = log(fishC[, "cpue"]), schaefer = TRUE, priorcalc = calcprior, scales = setscale)
cat("time = ", gettime() - beginC, " secs  \n")

par(mfrow=c(2, 2))
label <- c("r", "K", "Binit", "sigma")
for (i in 1:4) {
  plot(density(resultC$result[[2]][, i]), lwd = 2, col = 1,
       xlab = label[i], ylab = "Density", panel.first = grid())
  lines(density(resultC$result[[1]][, i]), lwd = 2, col = 2)
  lines(density(resultC$result[[3]][, i]), lwd = 2, col = 3)
}# These plots show 3 lines (black, green and red) of the densities of parameters accepted by the MCMC, showing a good level of convergence
# The process took 150 SECONDS THAT'S VERY LONG

# You can also generate summary stats, but I'm too lazy


### 6.8.4 Replicates Required to Avoid Serial Correlation ----

# If the thinning rate is too low, there'll be auto-correlation in the trace of each parameter
# How big the thinning rate needs to be is another question we need to answer

# There's a balance to be struck between high thinning rate and full variation of the different parameters

# Let's compare 2 thinning rates

setscale <- c(0.07, 0.06, 0.07, 0.45)
param <- log(c(R = 0.39, K = 9200, Binit = 3400, sigma = 0.05))

result1 <- do_MCMC(chains = 1, burnin = 100, N = 2000, thinstep = 1024, inpar = param, infunk = negLL1, calcpred = simpspmC, calcdat = fishC, obsdat = log(abdat$cpue), schaefer = TRUE, priorcalc = calcprior, scales = setscale)
result2 <- do_MCMC(chains = 1, burnin = 100, N = 1000, thinstep = 2048, inpar = param, infunk = negLL1, calcpred = simpspmC, calcdat = fishC, obsdat = log(abdat$cpue), schaefer = TRUE, priorcalc = calcprior, scales = setscale)

# The reason we need to remove within-sequence correlation is that it can interfere with stationary distribution convergence

# Below is a comparison of the two thinning rates, see how higher rate reduces the auto-correlation
posterior1 <- result1$result[[1]]
posterior2 <- result2$result[[1]]
label <- colnames(posterior1)[1:4]
par(mfrow = c(4, 2), mai = c(0.25, 0.45, 0.05, 0.05), oma = c(1, 0, 1, 0))
par(cex = 0.85, mgp = c(1.35, 0.35, 0), font.axis = 7, font = 7, font.lab = 7)
for (i in 1:4) {
  auto <- acf(posterior1[, i], type = "correlation", plot = TRUE, ylab = label[i], lag.max = 20, xlab = "", ylim = c(0, 0.3), lwd = 2)
  if (i == 1) mtext(1024, side = 3, line = -0.1, outer = FALSE, cex = 1.2)
  auto <- acf(posterior2[, i], type = "correlation", plot = TRUE, ylab = label[i], lag.max = 20, xlab = "", ylim = c(0, 0.3), lwd = 2)
  if (i == 1) mtext(2048, side = 3, line = -0.1, outer = FALSE, cex = 1.2)
}

parset(plots = c(2, 2), cex = 0.85)
label <- c("r", "K", "Binit", "sigma")
for (i in 1:4) {
  plot(density(result1$result[[1]][, i]), lwd = 4, col = 1, xlab = label[i], ylab = "Density", main = "", panel.first = grid())
  lines(density(result2$result[[1]][, i]), lwd = 2, col = 5, lty = 2)
} # Aaand a comparison of the two chains' (one has low thinning but high replicates, the other vice versa) acceptance densities for the 4 parameters
# They're pretty comparable, showing the balance of thinning rate and replicates.
# In general we would trust the chain with a higher thinning rate

# Another summary stats table I'm too lazy to copy



## 6.9 Concluding Remarks ----

# The whole point of this chapter was to highlight the amount of uncertainty in analysis
# So many little things you can adjust to increase or destroy the accuracy of things




# 7 Surplus Production Model ----

## 7.1 Introduction ----

# So far we've looked at static models (stable through time) and surplus production models (like the schaefer model)
# Surplus production models are dynamic models that use time-series data and can be used for stock assessment

# SPMs pool all aspects of production (recruitment, growth, mortality) into one function, dealing with undifferenciated biomass (age, size, sex and other differences are ignored)
# For formal stock assessments you need to model dynamic behaviour and productivity of an exploited stock - basically the stock's response to varied fishing pressure through time
# From this response to fishing intensity you can assess the stock's productivity.


### 7.1.1 Data Needs ----

# The minimum we need to estimate parameters is at least one time-series of relative abundance and an associated time-series of catch-data.
# Catch data can cover more years than the index data.
# Relative abundance can be CPUE, trawl survey, acoustic surveys etc.


### 7.1.2 The Need for Contrast ----

# SPMs assume that the stocks are in equilibrium which isn't the case really, so outcomes were too optimistic.
# They lacked contrast in effort levels and thus were uninformative about dynamics of populations.
# "lacking contrast" means that catch and effort is only available for a limited range of stock abundance levels and limited fishing intensity levels.
# Essentially, the less effort/catch data you have, the less range of responses you have to different fishing intensities,
# Meaning your stock could vary because of environmental factors, but it appears that the stock responds to the fishery in unexpected ways despite no changes in effort/catch. (And I'm guessing the model doesn't know how to deal with that)
# Basically more data = more contrast = more robust model. Makes sense

# A big assumption is that your index of relative abundance is an informative measure of stock abundance through time
# Very bery stable CPUEs can exist despite stocks declining, or highly variable CPUEs can be observed despite a stock not varying that much. Tech creep leading to effort creep can contribute


### 7.1.3 When are Catch-Rates Informative? ----

# Essentially you want to be sure that your index of relative abundance (ex. CPUE) and actual abundance have a logical relationship
# Example: if you allow effort to increase, CPUE over time will decrease (fishing above productivity), and vice versa
# If you can find such a relationship in the data it means there's probably some degree of contrast.

data("schaef") # In this tuna fishery example, let's see if CPUE and catch have enough of a relationship for CPUE to be used as a relative abundance index for an SPM

# This shows a correlation between catches and cpue for this dataset with a lag of 2 years (if you physically lag CPUE by 2 years, correlation is more apparent)
ccf(x = schaef[, "catch"], y = schaef[, "cpue"], type = "correlation", ylab = "Correlation", plot = TRUE)
# This suggests that there's enough contrast to inform an SPM

parset(plots = c(3, 1), margin = c(0.35, 0.4, 0.05, 0.05))
plot1(schaef[1:20, "year"], schaef[3:22, "catch"], ylab = "catch", xlab = "Year", defpar = F, lwd = 2)
plot1(schaef[3:22, "year"], schaef[3:22, "cpue"], ylab = "CPUE", xlab = "Year", defpar = F, lwd = 2)
plot1(schaef[1:20, "catch"], schaef[1:20, "cpue"], ylab = "CPUE", xlab = "catch", defpar = F, lwd = 2, type = "p")
model2 <- lm(schaef[3:22, "cpue"] ~schaef[1:20, "catch"])
abline(model2, lwd = 2, col = 2)

summary(model2)



## 7.2 Some Equations ----

# The relative abundance index we're using is assumed to reflect biomass available
# And the biomass is assumed to be affected by catches removed by the fishery.

# So if you use CPUE as a relative abundance index, you're dealing with exploitable biomass, not spawning biomass (??)

# See equations in book for calculating biomass and relative abundance
# This all assumes that any catch of the fishery will have an effect across the whole stock within a time step.
# So if you have big spatial structure in a stock or fishery, this is a big assumption...


### 7.2.1 Production Functions ----

# Different stocks will have different productivity shapes.
# Below is the Schaefer model and the Fox model, which may describe different stocks
# Their equations involve parameter p, which alters the shape of the dome

prodfun <- function(r, Bt, K, p) return((r*Bt/p)*(1-(Bt/K)^p))
densdep <- function(Bt, K, p) return((1/p)*(1-(Bt/K)^p))
r <- 0.75; K <- 1000 ; Bt <- 1:1000
sp <- prodfun(r, Bt, K, 1) # Schaefer production equivalent, p = 1
sp0 <- prodfun(r, Bt, K, p = 1e-08) # Fox production equivalent, p = 1e-08
sp3 <- prodfun(r, Bt, K, 3) # Left-skew, p = 3

parset(plots = c(2, 1))
plot1(Bt, sp, type = "l", lwd = 2, xlab = "Stock Size", ylab = "Surplus Production", maxy = 200, defpar = F)
lines(Bt, sp0*(max(sp)/max(sp0)), lwd = 2, col = 2, lty = 2)
lines(Bt, sp3*(max(sp)/max(sp3)), lwd = 3, col = 3, lty = 3)

plot1(Bt, densdep(Bt, K, p = 1), xlab = "Stock Size", defpar = FALSE, ylab = "Density-Dependence", maxy = 2.5, lwd = 2)
lines(Bt, densdep(Bt, K, 1e-08), lwd = 2, col = 2, lty = 2)
lines(Bt, densdep(Bt, K, 3), lwd = 3, col = 3, lty = 3)

# the general Schaefer model is symmetrical in production on either side of 0.5K (the max. of the curve being MSY)
# While the Fox model is skewed to one side, and the left-skew model to the right
# The red line (Fox) would be more suited to species that thrive at low population levels (like... sardines ? idk) = low density-dependence at low stock size
# Whereas the green line is more suited to species like marine mammals, whose stock growth would only decline at high stock sizes = high density dependence at low stock sizes (??)


### 7.2.2 The Schaefer Model ----

# For the schaefer model you'll need a carrying capacity/unfished biomass K (if the stock has been depleted at the start of the dataset, you'll also need initial biomass Binit,), reproductive rate r, and catchability coefficient q.
# Catchability can have one estimate (no big changes in CPUE) or multiple (if big changes in CPUE over the time series)


### 7.2.3 Sum of Squared Residuals ----

# We would fit the model using least squares (sum of squared residual errors)


### 7.2.4 Estimating Management Statistics ----

# As seen in the plots, MSY will be different depending on whether you use Schaefer or Fox equations.
# Below it the MSY for both
param <- c(r = 1.1, K = 1000, Binit = 800, sigma = 0.075)
cat("MSY Schaefer = ", getMSY(param, p = 1), "\n")
cat("MSY Fox      = ", getMSY(param, p = 1e-08), "\n")

# You can also obtain MSY-derived management advice like "the effort that would yield MSY"


### 7.2.5 The Trouble with Equilibria ----

# The idea of MSY and related stats revolve around the idea of equilibrium
# This rarely happens in real life. If the effort needed for MSY is applied, MSY can only be attained if Bmsy was the starting point, otherwise there's too much fishing leaving no biomass to replenish
# So watch out when using equilibrim statistics



## 7.3 Model Fitting ----

### 7.3.1 A Possible Workflow for Stock Assessment ----

# In order of things to do for stock assessment you might:
#  1. read in time series of catch and relative abundance data (check completeness, missing values etc.)
#  2. use ccf() to determine whether CPUE data relative to catch data is informative (sig. negative correlation would provide strong proof)
#  3. define/guess initial parameters
#  4. use plotsmmod() to plot the implications of assumed initial parameter sets for the dynamics - useful to search for initial parameter sets
#  5. use nlm() or fitSPM() to search for optimal parameter sets
#  6. use plotsmmod() again to see relative fit of optimal parameters (especially residual plot)
#  7. examine robustness of the model fit, using multiple initial parameter set starting points, see robustSPM()
#  8. once satisfied with the robustness, use spmphaseplot() to plot the phase diagram of biomass v. harvest rate to determine stock status visually
#  9. use spmboot() or other methods to assess uncertainty in the model fit and outputs
#  10. document and defend reached conclusions.

# Below is the initial model 'fit' to initial parameters
data(schaef); shaef <- as.matrix(schaef) # 1.
param <- log(c(r = 0.1, K = 2250000, sigma = 0.5)) # 3.
negatL <- negLL(param, simpspm, schaef, logobs = log(schaef[, "cpue"]))
ans <- plotspmmod(inp = param, indat = schaef, schaefer = T, addrmse = T, plotprod = F) # 6.

# Top right plot: green dash is simple loess fit, black line is from the guessed input parameters
# Bottom left plot: red line is predicted MSY
# Bottom right: strong residual pattern up to 1950ish, number on the top left is the root mean sq error of log normal residuals

# Then we get to actual model fitting, using two different methods (nlm and optim) (5.)
pnams <- c("r", "K", "Binit", "sigma")
best <- optim(par = param, fn = negLL, funk = simpspm, indat = schaef, logobs = log(schaef[, "cpue"]), method = "BFGS")
outfit(best, digits = 4, title = "Optim", parnames = pnams)
best2 <- nlm(negLL, best$par, funk = simpspm, indat = schaef, logobs = log(schaef[, "cpue"]))
outfit(best2, digits = 4, title = "nlm", parnames = pnams)
# The outputs give the same optimal parameters for both methods, which is a good precaution to take

# Now we take optimum parameters and visualise the fit (6.)
ans <- plotspmmod(inp = best2$estimate, indat = schaef, addrmse = T, plotprod = T)
# The output of plotspmmod is plots, but also invisible object we'll want to use later.

# Comparing parametric MSY with numerical MSY (you can obtain MSY from the optimal parameters -parametric MSY- or from the productivity curce -numerical MSY-)
round(ans$Dynamics$sumout, 3)
cat("\n Productivity Statistics \n")
summspm(ans)
# idk what this is


### 7.3.2 Is the Analysis Robust ? ----

# Numerical methods (which we are using to find parameter/model fit) can find false minima, which is why we used 2 methods above
# To test the robustness of the model fit we examine the influence initial parameters on the model fit using robustSPM()

data(schaef); schaef <- as.matrix(schaef); reps <- 12
param <- log(c(r = 0.15, K = 2250000, Binit = 2250000, sigma = 0.5))
ansS <- fitSPM(pars = param, fish = schaef, schaefer = T, maxiter = 1000, funk = simpspm, funkone = F)
set.seed(777852)
robout <- robustSPM(inpar = ansS$estimate, fish = schaef, N = reps, scaler = 40, verbose = F, schaefer = T, funk = simpspm, funkone = F)
# robout creates N vectors of parameters, each slightly different (scaler adjusts how different I think), each assigned a likelihood, then all N vectors are sorted by likelihood
# Let's repeat this 100 times

set.seed(777854)
robout2 <- robustSPM(inpar = ansS$estimate, fish = schaef, N = 100, scaler = 25, verbose = F, schaefer = T, funk = simpspm, funkone = F, steptol = 1e-06)
lastbits <- tail(robout2$results[, 6:11], 10)
# Essentially by increasing the number of reps, you can see how similar the likelihood is no matter how much you change vectors of parameters

result <- robout2$results
parset(plots = c(2, 2), margin = c(0.35, 0.45, 0.05, 0.05))
hist(result[, "r"], breaks = 15, col = 2, xlab = "r")
hist(result[, "K"], breaks = 15, col = 2, xlab = "K")
hist(result[, "Binit"], breaks = 15, col = 2, xlab = "Binit")
hist(result[, "MSY"], breaks = 15, col = 2, xlab = "MSY")
# Here plots don't work like in the book but it's meant to show parameters are close together but there's still variation.


### 7.3.3 Using Different Data ----

# Same process as before, different data
set.seed(777854)
data("dataspm"); fish <- dataspm
param <- log(c(r = 0.24, K = 5174, Binit = 2846, sigma = 0.164))
ans <- fitSPM(pars = param, fish = fish, schaefer = T, maxiter = 1000, funkone = T)
out <- robustSPM(ans$estimate, fish, N = 100, scaler = 15, verbose = F, funkone = T)
result <- tail(out$results[, 6:11], 10) # Notice how the last 6 trials differ quite a bit (see -veLL) from the others, with parameters quite different
# I'm guessing this means that this model isn't a robust as the one we had before ?



## 7.4 Uncertainty ----

# We've had a look at how different parameter values can give us the same-ish numerical fit.
# We now need to know the precision and uncertainty of those parameter sets.


## 7.4.1 Likelihood Profiles ----

# Likelihood profiles provide insight into how model fit quality changes with slightly different parameters
# A model is optimally fitted with likelihood methods, then some parameters are fixed to pre-determined values, and we fit the model with only one parameter changing.

data(abdat); fish <- as.matrix(abdat)
colnames(fish) <- tolower(colnames(fish))
pars <- log(c(r = 0.4, K = 9400, Binit = 3400, sigma = 0.05))
ans <- fitSPM(pars, fish, schaefer = T)
answer <- plotspmmod(ans$estimate, abdat, schaefer = T, addrmse = T)

doprofile <- function(val, loc, startest, indat, notfix = c(2:4)) {
  pname <- c("r", "K", "Binit", "sigma", "-veLL")
  numv <- length(val)
  outpar <- matrix(NA, nrow = numv, ncol = 5, dimnames = list(val, pname))
  for (i in 1:numv) {
    param <- log(startest)
    param[loc] <- log(val[i])
    parinit <- param
    bestmod <- nlm(f = negLLP, p = param, funk = simpspm, initpar = parinit, indat = indat, logobs = log(indat[, "cpue"]), notfixed = notfix)
    outpar[i, ] <- c(exp(bestmod$estimate), bestmod$minimum)
  }
  return(outpar)
}

rval <- seq(0.32, 0.46, 0.001)
outr <- doprofile(rval, loc = 1, startest = c(rval[1], 11500, 5000, 0.25), indat = fish, notfix = c(2:4))
Kval <- seq(7200, 11500, 200)
outk <- doprofile(Kval, loc = 2, c(0.4, 7200, 6500, 0.3), indat = fish, notfix = c(1, 3, 4))

parset(plots = c(2, 1), cex = 0.85, outmargin = c(0.5, 0.5, 0, 0))
plotprofile(outr, var = "r", defpar = F, lwd = 2)
plotprofile(outk, var = "K", defpar = F, lwd = 2)
# Likelihood profiles of parameters r and K of the Schaefer model fit to abdat data. middle green line is the -veLL minimum (most likely value) and thick green lines are the 95% CI.

# When you look at single parameters you ignore the correlation between parameters (which Schaefer is known to have) so be careful
# When you do likelihood profiles, change the initial parameters and re-run the profile to see how consisten the model fit is with the fixed parameters


### 7.4.2 Bootstrap Confidence Intervals ----

# One way to characterise uncertainty in a model fit is to put CI around parameters and model outputs by bootstrapping residuals:
# Essentially you take bootstrap samples of the log-normal residuals associated with CPUE and use those to generate new bootstrap CPUE samples to replace original CPUE time-series.
# Each new bootstrap sample is re-fit and solutions stored for further analysis

data(dataspm)
fish <- as.matrix(dataspm)
colnames(fish) <- tolower(colnames(fish))
pars <- log(c(r = 0.25, K = 5500, Binit = 3000, sigma = 0.25))
ans <- fitSPM(pars, fish, schaefer =T, maxiter = 1000)
answer <- plotspmmod(ans$estimate, fish, schaefer = T, addrmse = T)
# Okay that's the initial optimum fit, let's bootstrap

set.seed(210368)
reps <- 1000
boots <- spmboot(ans$estimate, fishery = fish, iter = reps)
str(boots, max.level = 1)
# Output is the dynamics of each bootstrap run, with predicted model biomass, bootstrap CPUE sample, predicted CPUE for each bootstrap sample, depletion time-series, annual harvest rate time-series

bootpar <- boots$bootpar
rows <- colnames(bootpar)
columns <- c(c(0.025, 0.05, 0.5, 0.95, 0.975), "Mean")
bootCI <- matrix(NA, nrow = length(rows), ncol = length(columns), dimnames = list(rows, columns))
for (i in 1:length(rows)) {
  tmp <- bootpar[, i]
  qtil <- quantile(tmp, probs = c(0.025, 0.05, 0.5, 0.95, 0.975), na.rm = T)
  bootCI[i, ] <- c(qtil, mean(tmp, na.rm = T))
}
bootCI # The summary of each parameter's confidence interval, using bootstrap CI

colf <- c(1, 1, 1, 4); lwdf <- c(1, 3, 1, 3); ltyf <- c(1, 1, 1, 2)
colsf <- c(2, 3, 4, 6); parset(plots = c(3, 2))

hist(bootpar[, "r"], breaks = 25, main = "", xlab = "r")
abline(v = c(bootCI["r", colsf]), col = colf, lwd = lwdf, lty = ltyf)

uphist(bootpar[, "K"], maxval = 14000, breaks = 25, main = "", xlab = "K")
abline(v = c(bootCI["K", colsf]), col = colf, lwd = lwdf, lty = ltyf)

hist(bootpar[, "Binit"], breaks = 25, main = "", xlab = "Binit")
abline(v = c(bootCI["Binit", colsf]), col = colf, lwd = lwdf, lty = ltyf)

uphist(bootpar[, "MSY"], maxval = 14000, breaks = 25, main = "", xlab = "MSY")
abline(v = c(bootCI["MSY", colsf]), col = colf, lwd = lwdf, lty = ltyf)

hist(bootpar[, "Depl"], breaks = 25, main = "", xlab = "Final Depletion")
abline(v = c(bootCI["Depl", colsf]), col = colf, lwd = lwdf, lty = ltyf)

hist(bootpar[, "Harv"], breaks = 25, main = "", xlab = "End Harvest Rate")
abline(v = c(bootCI["Harv", colsf]), col = colf, lwd = lwdf, lty = ltyf)
# Boom, plots of the CI around each parameter and output

dynam <- boots$dynam
years <- fish[, "year"]
nyrs <- length(years)
parset()
ymax <- getmax(c(dynam[,, "predCE"], fish[, "cpue"]))
plot(fish[, "year"], fish[, "cpue"], type = "n", ylim = c(0, ymax), xlab = "Year", ylab = "CPUE", yaxs = "i", panel.first = grid())
for (i in 1:reps) lines(years, dynam[i,, "predCE"], lwd = 1, col = 8)
lines(years, fish[, "cpue"], cex = 1.2, pch = 16, col = 1)
points(years, fish[, "cpue"], cex = 1.2, pch = 16, col = 1)
percs <- apply(dynam[,, "predCE"], 2, quants)
arrows(x0 = years, y0 = percs["5%",], y1 = percs["95%",], length = 0.03, angle = 90, code = 3, col = 1)
# Above is a plot of: original observed CPUE (black dots), optimum predicted CPUE (solid line), 1000 bootstrap predicted CPUE (gray lines), and 90th percentile CI around the predicted values

# We've used Schaefer so far, let's try fitting it to Fox

pars <- log(c(r = 0.15, K = 6500, Binit = 3000, sigma = 0.2))
ansF <- fitSPM(pars, fish, schaefer = F, maxiter = 1000)
bootsF <- spmboot(ansF$estimate, fishery = fish, iter = reps, schaefer = F)
dynamF <- bootsF$dynam

parset()
ymax <- getmax(c(dynam[,, "predCE"], fish[, "cpue"]))
plot(fish[, "year"], fish[, "cpue"], type = "n", ylim = c(0, ymax), xlab = "Year", ylab = "CPUE", yaxs = "i", panel.first = grid())
for (i in 1:reps) lines(years, dynamF[i,, "predCE"], lwd = 1, col = 1)
for (i in 1:reps) lines(years, dynam[i,, "predCE"], lwd = 1, col = 8)
lines(years, fish[, "cpue"], cex = 1.2, pch = 16, col = 1)
points(years, fish[, "cpue"], cex = 1.2, pch = 16, col = 1)
percs <- apply(dynam[,, "predCE"], 2, quants)
arrows(x0 = years, y0 = percs["5%",], y1 = percs["95%",], length = 0.03, angle = 90, code = 3, col = 1)
# Same plot as above, with thin black lines being Fox bootstrap predicted CPUE
# You could argue Fox captures the variability of the data more than Schaefer, but it is more uncertain


### 7.4.3 Parameter Correlations ----

# Parameters may be correlated, plot them against each other to test

pairs(boots$bootpar[, c(1:4, 6, 7)], lower.panel = panel.smooth, upper.panel = panel.cor)
# Strong relationship between r and K


### 7.4.4 Asymptotic Errors ----

# Asymptotic errors are derived from cariance-covariance matrix, can be used to describe variability and interactions between parameters of a model

data(dataspm)
fish <- as.matrix(dataspm)
colnames(fish) <- tolower(colnames(fish))
pars <- log(c(r = 0.25, K = 5200, Binit = 2900, sigma = 0.2))
ans <- fitSPM(pars, fish, schaefer = T, maxiter = 1000, hess = T)
outfit(ans)

vcov <- solve(ans$hessian)
label <- c("r", "K", "Binit", "sigma")
colnames(vcov) <- label; rownames(vcov) <- label
outvcov <- rbind(vcov, sqrt(diag(vcov)))
rownames(outvcov) <- c(label, "StErr")

library(mvtnorm)
N <- 1000
mvn <- length(fish[, "year"])
mvncpue <- matrix(0, nrow = N, ncol = mvn, dimnames = list(1:N, fish[, "year"]))
columns <- c("r", "K", "Binit", "sigma")
optpar <- ans$estimate
mvnpar <- matrix(exp(rmvnorm(N, mean = optpar, sigma = vcov)), nrow = N, ncol = 4, dimnames = list(1:N, columns))
msy <- mvnpar[, "r"]*mvnpar[, "K"]/4
nyr <- length(fish[, "year"])
depletion <- numeric(N)
for (i in 1:N) {
  dynamA <- spm(log(mvnpar[i, 1:4]), fish)
  mvncpue[i, ] <- dynamA$outmat[1:nyr, "predCE"]
  depletion[i] <- dynamA$outmat["2016", "Depletion"]
}
mvnpar <- cbind(mvnpar, msy, depletion)

plot1(fish[, "year"], fish[, "cpue"], type = "p", xlab = "Year", ylab = "CPUE", maxy = 2)
for (i in 1:N) lines(fish[, "year"], mvncpue[i ,], col = "gray", lwd = 1)
points(fish[, "year"], fish[, "cpue"], cex = 1.2, col = 1)
lines(fish[, "year"], exp(simpspm(optpar, fish)), cex = 1.2, pch = 16, col = 1)
percs <- apply(mvncpue, 2, quants)
arrows(x0 = fish[, "year"], y0 = percs["5%",], y1 = percs["95%",], length = 0.03, angle = 90, code = 3, col = 1)
msy <- mvnpar[, "r"]*mvnpar[, "K"]/4
# Similar CI to bootstrap, but a lot more skewed than bootstrap.

# If we look for parameter sets that led to these extreme results:
pickd <- which(mvncpue[, "2016"] < 0.4)
plot1(fish[, "year"], fish[, "cpue"], type = "n", xlab = "Year", ylab = "CPUE", maxy = 6.25)
for (i in 1:length(pickd))
  lines(fish[, "year"], mvncpue[pickd[i],], col = 1, lwd = 1)
points(fish[, "year"], fish[, "cpue"], pch = 16, cex = 1.25, col = 4)
lines(fish[, "year"], exp(simpspm(optpar, fish)), lwd = 3, col = 2, lty = 2)

# We can compare the parameter sets

parset(plots = c(2, 2), cex = 0.85)
outplot <- function(var1, var2, pickdev) {
  plot1(mvnpar[, var1], mvnpar[, var2], type = "p", pch = 16, cex = 1, defpar = F, xlab = var1, ylab = var2, col = 8)
  points(mvnpar[pickdev, var1], mvnpar[pickdev, var2], pch = 16, cex = 1)
}
outplot("r", "K", pickd)
outplot("sigma", "Binit", pickd)
outplot("r", "Binit", pickd)
outplot("K", "Binit", pickd)
# In these plots, the cloud of points is all the predicted CPUEs, and black dots are the ones that are weird.
# For all these runs, parameters are plotted against each other, and the common denominator for all of these plots is that low Binit values = weird predicted CPUEs

# We can compare bootstrap and asymptotic error parameter sets
bt <- apply(bootpar, 2, range)[, c(1:4, 6, 7)]
ay <- apply(mvnpar, 2, range)
out <- rbind(bt, ay)
rownames(out) <- c("MinBoot", "MaxBoot", "MinAsym", "MaxAsym")
out # Quite different


### 7.4.5 Sometimes Asymptotic Errors Work ----

# Sometimes AE (asym. errors) give very similar results to BS (bootstrap)

data(abdat)
fish <- as.matrix(abdat)
pars <- log(c(r = 0.4, K = 9400, Binit = 3400, sigma = 0.05))
ansA <- fitSPM(pars, fish, schaefer = T, maxiter = 1000, hess = T)
vcovA <- solve(ansA$hessian)
mvn <- length(fish[, "year"])
N <- 1000
mvncpueA <- matrix(0, nrow = N, ncol = mvn, dimnames = list(1:N, fish[, "year"]))
columns <- c("r", "K", "Binit", "sigma")
optparA <- ansA$estimate
mvnparA <- matrix(exp(rmvnorm(N, mean = optparA, sigma = vcovA)), nrow = N, ncol = 4, dimnames = list(1:N, columns))
msy <- mvnparA[, "r"]*mvnparA[, "K"]/4
for (i in 1:N) mvncpueA[i, ] <- exp(simpspm(log(mvnparA[i, ]), fish))
mvnparA <- cbind(mvnparA, msy)

plot1(fish[, "year"], fish[, "cpue"], type = "p", xlab = "Year", ylab = "CPUE", maxy = 2.5)
for (i in 1:N) lines(fish[, "year"], mvncpueA[i, ], col = 2, lwd = 1)
points(fish[, "year"], fish[, "cpue"], pch = 16, cex = 1)
lines(fish[, "year"], exp(simpspm(optparA, fish)), lwd = 2, col = 0)
# AE to generate plausible parameter sets, and their implied CPUE trajectories. Optimum is white line


### 7.4.6 Bayesian Posteriors ----

data(abdat); fish <- as.matrix(abdat)
param <- log(c(r = 0.3, K = 11500, Binit = 3300, sigma = 0.05))
foxmod <- nlm(f = negLL1, p = param, funk = simpspm, indat = fish, logobs = log(fish[, "cpue"]), iterlim = 1000, schaefer = F)
optpar <- exp(foxmod$estimate)
ans <- plotspmmod(inp = foxmod$estimate, indat = fish, schaefer = F, addrmse = T, plotprod = T)

# Let's try MCMC
library(Rcpp)
cppFunction('NumericVector simpspmC(NumericVector pars, NumericMatrix indat, LogicalVector schaefer) {
            int nyrs = indat.nrow();
            NumericVector predce(nyrs);
            NumericVector biom(nyrs+1);
            double Bt, qval;
            double sumq = 0.0;
            double p = 0.00000001;
            if (schaefer(0) == TRUE) {
              p = 1.0;
            }
            NumericVector ep = exp(pars);
            biom[0] = ep[2];
            for (int i = 0; i < nyrs; i++) {
              Bt = biom[i];
              biom[(i+1)] = Bt+(ep[0]/p)*Bt*(1-pow((Bt/ep[1]),p))-indat(i,1);
            if (biom[(i+1)] < 40.0) biom[(i+1)] = 40.0;
            sumq += log(indat(i,2)/biom[i]);
            }
            qval = exp(sumq/nyrs);
            for (int i = 0; i < nyrs; i++) {
              predce[i] = log(biom[i] * qval);
            }
            return predce;
}')
set.seed(698381) # to have same results as book , omit irl
begin <- gettime()
inscale <- c(0.07, 0.05, 0.09, 0.45)
pars <- log(c(r = 0.205, K = 11300, Binit = 3200, sigma = 0.044))
result <- do_MCMC(chains = 1, burnin = 50, N = 2000, thinstep = 512, infunk = negLL, calcpred = simpspmC, obsdat = log(fish[, "cpue"]), inpar = pars, calcdat = fish, priorcalc = calcprior, scales = inscale, schaefer = F)
cat("acceptance rate", result$arate, "\n")
cat("time", gettime() -begin, "\n")
post1 <- result[[1]][[1]]
p <- 1e-08
msy <- post1[, "r"]*post1[, "K"]/((p + 1)^((p+1)/p))

# Marginal distribution
parset(plots = c(2, 2), cex = 0.85)
plot(density(post1[, "r"]), lwd = 2, main = "", xlab = "r")
plot(density(post1[, "K"]), lwd = 2, main = "", xlab = "K")
plot(density(post1[, "Binit"]), lwd = 2, main = "", xlab = "Binit")
plot(density(msy), lwd = 2, main = "", xlab = "MSY")
# Lumpiness indicates more than 2000 iterations are needed

puttxt <- function(xs, xvar, ys, yvar, lvar, lab = "", sigd = 0) {
  text(xs*xvar[2], ys*yvar[2], makelabel(lab, lvar, sep = "  ", sigdig = sigd), cex = 1.2, font = 7, pos = 4)
}
kran <- range(post1[, "K"]); rran <- range(post1[, "r"])
mran <- range(msy)

parset(plots = c(1, 2), margin = c(0.35, 0.35, 0.05, 0.1))
plot(post1[, "K"], post1[, "r"], type = "p", cex = 0.5, xlim = kran, ylim = rran, col = "gray", xlab = "K", ylab = "r", panel.first = grid())
points(optpar[2], optpar[1], pch = 16, col = 1, cex = 1.75)
addcontours(post1[, "K"], post1[, "r"], kran, rran, contval = c(0.5, 0.9), lwd = 2, col = 1)

plot(post1[, "K"], msy, type = "p", cex = 0.5, xlim = kran, ylim = mran, col = "gray", xlab = "K", ylab = "MSY", panel.first = grid())
points(optpar[2], getMSY(optpar, p), pch = 16, col = 1, cex = 1.75)
addcontours(post1[, "K"], msy, kran, rran, mran, contval = c(0.5, 0.9), lwd = 2, col = 1) # Doesn't work fsr

# Anyway there are MCMC marginal distribution outputs as scattergrams of r and K parameters, and K and MSY values.
# Contours are 50th and 90th percentiles. Dots are successful parameter vectors

parset(plots = c(4, 1), margin = c(0.3, 0.45, 0.05, 0.05), outmargin = c(1, 0, 0, 0))
label <- colnames(post1)
N <- dim(post1)[1]
for (i in 1:3) {
  plot(1:N, post1[, i], type = "l", ylab = label[i])
  abline(h = median(post1[, i]), col = 2)
}

msy <- post1[, 1]*post1[, 2]/4
plot(1:N, msy, type = "l", lwd = 1)
abline(h = median(msy), col = 2)
mtext("Step", side = 1, outer = T, line = 0, font = 7, cex = 1)
# Traces of the MCMC, 3 model parameters sets, and MSY estimates

set.seed(6396679)
inscale <- c(0.07, 0.05, 0.09, 0.45)
pars <- log(c(r = 0.205, K = 11300, Binit = 3200, sigma = 0.044))
result <- do_MCMC(chains = 5, burnin = 50, N = 2000, thinstep = 512, inpar = pars, infunk = negLL1, calcpred = simpspmC, obsdat = log(fish[, "cpue"]), calcdat = fish, priorcalc = calcprior, scales = inscale, schaefer = F)
cat("acceptance rate", result$arate, "\n")

parset(plots = c(2, 1))
post <- result[[1]][[1]]
plot(density(post[, "K"]), lwd = 2, col = 1, ylim = c(0, 4.4e-04))
for(i in 2:5) lines(density(result$result[[i]][, "K"]), lwd = 2, col = i)
p <- 1e-08
post <- result$result[[1]]
msy <- post[, "r"]*post[, "K"]/((p+1)^((p+1)/p))
plot(density(msy), lwd = 2, col = 1, ylim = c(0, 0.0175))
for(i in 2:5) {
  post <- result$result[[i]]
  msy <- post[, "r"]*post[, "K"]/((p+1)^((p+1)/p))
  lines(density(msy), lwd = 2, col = i)
}
# Marginal posterior for the K parameter and implied MSY from the 5 chains, 2000 reps



## 7.5 Management Advice ----

### 7.5.1 Two Views of Risk ----

# Once you have a fit model, how do you advise management?
# We can model future dynamics.
# If you have a fishery goal, target reference point and limit reference point you can solve effort or catch for that goal


### 7.5.2 Harvest Strategies ----

# Harvest strategy defines decision framework. Three components:
# 1. Means of monitoring and collecting data for each fishery of interest
# 2. Defined manner in which each fishery is to be assessed (relative to biological ref points like mortality etc.)
# 3. Pre-defined harvest control rules to translate stock status into management advice



## 7.6 Risk Assessment Projections ----

# There are big assumptions in stock assessment models, including that:
# You're successful at capturing the important parts of the dynamics controlling the stock
# You assume that stock productivity will remain the same int the future

# If an assessment is highly uncertain then future projections will also be highly uncertain


### 7.6.1 Deterministic Projections ----

library(MQMF)
data(abdat); fish <- as.matrix(abdat)
param <- log(c(r = 0.3, K = 11500, Binit = 3300, sigma = 0.05))
bestmod <- nlm(f = negLL1, p = param, funk = simpspm, schaefer = F, logobs = log(fish[, "cpue"]), indat = fish, hessian = T)
optpar <- exp(bestmod$estimate)
ans <- plotspmmod(inp = bestmod$estimate, indat = fish, schaefer = F, target = 0.4, addrmse = T, plotprod = F)
# Okay, best fit model, all good, let's project future biomass

out <- spm(bestmod$estimate, indat = fish, schaefer = F)
str(out, width = 65, strict.width = "cut")

catches <- seq(700, 1000, 50) # Our scenario cathces for projection
projans <- spmprojDet(spmobj = out, projcatch = catches, plotout = T)


### 7.6.2 Accounting for Uncertainty ----

# The projections are deterministic, they fail to account for uncertainty


### 7.6.3 Using Asymptotic Errors ----

# To get asymptotic errors we need multiple plausible parameter vectors from optimum model fit
library(mvtnorm)
marpar <- parasympt(bestmod, N = 1000)
projs <- spmproj(marpar, fish, projyr = 10, constC = 900)
outp <- plotproj(projs, out, qprob = c(0.1, 0.5), refpts = c(0.2, 0.4))
# This is 1000 projections derived from inverse hessian and mean parameter estimates to generate 1000 plausible parameter vectors and projecting each vector forward
# Blue line is fisheries data limit, thick black is optimum fit, dashed are target and limit ref points
# What the plot says is that with sustained catches, you can say on average that stock will decline, but still within target


### 7.6.4 Using Bootstrap Parameter Vectors ----

reps <- 1000
boots <- spmboot(bestmod$estimate, fishery = fish, iter = reps, schaefer = F)
matparb <- boots$bootpar[, 1:4]
projb <- spmproj(matparb, fish, projyr = 10, constC = 900)
outb <- plotproj(projb, out, qprob = c(0.1, 0.5), refpts = c(0.2, 0.4))
# Slight difference from AEs, but still pretty similar


### 7.6.5 Using Samples from a Bayesian Posterior ----

param <- log(c(r = 0.3, K = 11500, Binit = 3300, sigma = 0.05))
set.seed(444608)
N <- 1000
result <- do_MCMC(chains = 1, burnin = 100, N = N, thinstep = 2048, inpar = param, infunk = negLL, calcpred = simpspmC, calcdat = fish, obsdat = log(fish[, "cpue"]), priorcalc = calcprior, schaefer = F, scales = c(0.065, 0.055, 0.1, 0.475))
parB <- result[[1]][[1]]
cat("Acceptance rate = ", result[[2]], "\n")
parset(plots = c(2, 1))
acf(parB[, 2], lwd = 2)
plot(1:N, parB[, 2], type = "l", ylab = "K", ylim = c(8000, 19000))
# Lack of autocorrelation, only few spikes on lower plot

matparB <- as.matrix(parB[, 1:4])
projs <- spmproj(matparB, fish, constC = 900, projyr = 10)
plotproj(projs, out, qprob = c(0.1, 0.5), refpts = c(0.2, 0.4))
# 1000 projections of constant 900t catch, derived from using 1000 samples from the Bayesian posterior



## 7.7 Concluding Remarks ----

# ...
# Holy shit I finished the book

