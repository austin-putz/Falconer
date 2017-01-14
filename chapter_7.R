
#=========================================================================================#
# 
# Chapter 7 Pygmy example
# 
#=========================================================================================#

# a and d defined
# alpha = a + d(q - p)
# mean = p^2*a + 2pqd + q^2*-a
# mean = a(p - q) + 2dpq, simplified

#=========================================================================================#
# Libraries
#=========================================================================================#

  library(ggplot2)
  library(lsmeans)
  library(car)
  library(dplyr)
  library(pander)
    panderOptions('table.split.table', 150)

#=========================================================================================#
# Function
#=========================================================================================#

  createData <- function(n=1000, q=0.5, values, sd=1) {
  	
  	if (q < 0 | q > 1) stop("q needs to be between 0 and 1")
  	
  	# allele freqs
  	p  <- 1 - q
  	
  	# genotype freqs
  	AA <- p^2
  	Aa <- 2 * p * q 
  	aa <- q^2
  	
  	print(AA); print(Aa); print(aa)
  	
  	# proportions of N
  	AA_n <- n * AA
  	Aa_n <- n * Aa
  	aa_n <- n * aa

  	print(AA_n); print(Aa_n); print(aa_n)
  	
  	# sample from random distn to get values
  	AA_samples <- rnorm(values[3], sd=sd, n=AA_n)
  	Aa_samples <- rnorm(values[2], sd=sd, n=Aa_n)
  	aa_samples <- rnorm(values[1], sd=sd, n=aa_n)
  	
  	# create vectors
  	values <- c(AA_samples, Aa_samples, aa_samples)
  	copies  <- c(rep(2, AA_n), rep(1, Aa_n), rep(0, aa_n))
  	
  	# create dataset
  	data <- data.frame(copies, values)
  	
  	# return data frame
  	return(data)
  	
  }

#=========================================================================================#
# Data
#=========================================================================================#

# create dataset
  pygmy <- createData(n=1000, q=0.1, c(6, 12, 14), sd=2)
  pander(some(pygmy))

# center phenotypic values
  pygmy$centered <- pygmy$values - mean(pygmy$values)

# plot of data
  plot(pygmy$copies, pygmy$values, xlab="Copies", ylab="Phenotype")
  plot(pygmy$copies, pygmy$centered, xlab="Copies", ylab="Phenotype - Mean")

# true mean values (see book)
  true_values <- c(6, 12, 14)

# aggregate to get estimated means
  pygmy.agg <- pygmy %>% group_by(copies) %>% summarise(mean(values))
  names(pygmy.agg) <- c("copies", "mean")
  pygmy.agg

# add true_values
  pygmy.agg$mean_true <- true_values

# allele freqs
  p <- mean(pygmy$copies) / 2
  print(p)
  q <- 1 - p
  print(q)

# true parameters
  m <- mean(true_values[c(1,3)])
  m
  a <- true_values[3] - m
  a
  d <- true_values[2] - m
  d

# estimated parameters
  m_est <- mean(pygmy.agg$mean[c(1,3)])
  m_est
  a_est <- pygmy.agg$mean[3] - mean(pygmy.agg$mean[c(1,3)])
  a_est
  d_est <- pygmy.agg$mean[2] - mean(pygmy.agg$mean[c(1,3)])
  d_est

# another centered column, but not from mean of homozygotes
  pygmy$centered_2 <- pygmy$values - mean(pygmy.agg$mean[c(1,3)])

# genotypic values, true and estimated
  pygmy.agg$G_true <- c(true_values[1] - m, true_values[2] - m, true_values[3] -m)
  pygmy.agg$G_est  <- c(pygmy.agg$mean[1] - m_est, pygmy.agg$mean[2] - m_est, pygmy.agg$mean[3] - m_est)
  pygmy.agg

# means of population
  M <- a * (p - q) + (2 * d * p * q)
  M
  M_est <- a_est * (p - q) + (2 * d_est * p * q)
  M_est
  pop_mean <- m + M
  pop_mean

# add freq
  pygmy.agg$freq <- c(q^2, 2*p*q, p^2)
  pygmy.agg

# alpha of true values
	alpha <- a + (d * (q - p))  
  alpha

# alpha from estimated values
  alpha_est <- a_est + (d_est * (q - p)) # estimate of that reg coef
  alpha_est

#=========================================================================================#
# Model
#=========================================================================================#

# estimate a with true and estimates (DISREGARDS FREQUENCIES!!!!!)
  pygmy.true.lm <- lm(G_true ~ copies, data=pygmy.agg)
  pygmy.true.lm
  pygmy.est.lm  <- lm(G_est  ~ copies, data=pygmy.agg)
  pygmy.est.lm

# estimate alpha from dataset (frequencies included)
# this is the BV line
  pygmy.lm <- lm(centered ~ copies, data=pygmy)
  pygmy.lm

#=========================================================================================#
# Other
#=========================================================================================#

# add breeding values with true and estimated alpha
  pygmy.agg$bv_true <- c(-2*p*alpha, (q-p)*alpha, 2*q*alpha)
  pygmy.agg$bv_est  <- c(-2*p*alpha_est, (q-p)*alpha_est, 2*q*alpha_est)
  pygmy.agg

# Dominance deviation
  pygmy.agg$D_true <- pygmy.agg$G_true - pygmy.agg$bv_true
  pygmy.agg$D_est  <- pygmy.agg$G_est  - pygmy.agg$bv_est
  pygmy.agg

# now get alpha from bv estimates
  pygmy.bv_true.lm <- lm(bv_true ~ copies, data=pygmy.agg)
  pygmy.bv_true.lm
  pygmy.bv_est.lm  <- lm(bv_est ~ copies, data=pygmy.agg)
  pygmy.bv_est.lm

#=========================================================================================#
# Plots
#=========================================================================================#

# plot centered phenotypic values vs copies
  plot(pygmy$copies, pygmy$centered, xlab="Copies", ylab="Phenotype - Mean")
  abline(a=pygmy.lm[1], b=pygmy.lm[2], 
  			 col="red", cex=1)
  text(x=0.5, y=0, "red is alpha line")

# plot
  plot(pygmy.agg$copies, pygmy.agg$G_true, type="p", 
  		 col="blue", pch=20, cex=2, ylim=c(-6, 10))        # true G values
  points(pygmy.agg$copies, pygmy.agg$G_est, type="p",
  			 col="steelblue", pch=20, cex=2, ylim=c(-6,6))  # est G values
  abline(a=pygmy.true.lm[1], b=pygmy.true.lm[2], 
  			 col="blue", cex=2)                 # true genotype line
  abline(a=pygmy.est.lm[1], b=pygmy.est.lm[2], 
  			 col="steelblue", cex=2)            # est genotype line
  abline(a=pygmy.lm[1], b=pygmy.lm[2], 
  			 col="red", cex=1)                  # esimated alpha from dataset
  abline(a=pygmy.bv_true.lm[1], b=pygmy.bv_true.lm[2],
  			 col="purple", cex=1)               # alpha from true values
  text(x=0.5, y=5.8, "blue is true genotype line (a)")
  text(x=0.5, y=4,  "steelblue is estimated genotype line (a)")
  text(x=0.5, y=-5,  "red is estimated alpha line")
  text(x=0.5, y=-5.8, "purple is true alpha line")


#   abline(a=pygmy.agg.lm$coef[1], b=pygmy.agg.lm$coef[2], col="blue", cex=2)
#   points(pygmy.agg$copies, pygmy.agg$bv, col="red", pch=20, cex=2)
#   abline(a=pygmy.agg.bv.lm$coef[1], b=pygmy.agg.bv.lm$coef[2], col="red", cex=2)
#   points(pygmy.agg$copies, pygmy.agg$D, col="orange", pch=20, cex=2)

# plot of genotypic and breeding values
  ggplot(data=pygmy.agg, aes(x=copies, y=G_true)) +
	geom_point(size=5, color="blue",   aes(y= G_est)) +
	geom_point(size=5, color="red",    aes(y=bv_true)) +
	geom_point(size=5, color="orange", aes(y=D_true)) +
	geom_smooth(data=pygmy.agg, method="lm", color="blue", aes(y=G_est), se=FALSE) +
	geom_smooth(data=pygmy.agg, method="lm", color="red", aes(y=bv_true), se=FALSE) +
	geom_smooth(method="lm", color="orange", aes(y=D_true), se=FALSE) +	
	scale_x_continuous(breaks=c(0,1,2), labels=c("0","1","2"), limits=c(0, 2)) +
	annotate("text", x=1.5, y=0.7,  label=paste("Population Mean =", pop_mean), size=5) +
    annotate("text", x=1.5, y=0,    label=paste("q =", q), size=5) +
	annotate("text", x=1.5, y=-0.7, label=paste("Slope =", pygmy.lm$coef[2]),    size=5) +
	annotate("text", x=1.5, y=-1.4, label=paste("allele substitution =", alpha), size=5) +
	theme(text=element_text(size=18))

# true values
  ggplot(data=pygmy.agg, aes(x=copies, y=G_est)) +
	geom_point(size=5, color="blue", aes(size=freq)) +
	geom_point(size=5, aes(y=bv_est), color="red") +
	geom_point(size=5, aes(y=D_est), color="orange") 



#=========================================================================================#
# Do permutation of data to get 
#=========================================================================================#

  pygmy.10.lm <- lm(values ~ copies, data=pygmy)
  summary(pygmy.10.lm)
  tstat <- as.data.frame(summary(pygmy.10.lm)[4])[2,3]
  
  permuteG <- function(data, x, y, n.rounds, plot=FALSE){
  	
  	# subset data
  	data = data[, c(y,x)]
  	
  	print(head(data))
  	
  	# set up vector
  	tvalue_vec = c(rep(-999, n.rounds))
  	
  	for (i in 1:n.rounds){
  		
  	  	# permute
  		newdata = data[sample(1:nrow(data), nrow(data), replace=FALSE), ]
  		
  		# cbind them
  		dataPerm = as.data.frame(cbind(data[, y], newdata[, x]))
  		names(dataPerm) = c(y, x)
  		
  		# fit lm
  		xy.lm = lm(as.numeric(get(y)) ~ as.numeric(get(x)), data= dataPerm)
  		
  		# remove table
  		table = as.data.frame(summary(xy.lm)[4])
  		
  		# remove t-stat
  		tstat = table[2, 3]
  		
  		# put into vector
  		tvalue_vec[i] = tstat
  		
  	}
  	
  	# histogram
	if (plot == TRUE) hist(tvalue_vec)
  	
  	# return vector
  	return(tvalue_vec)
  	
  }

# do permutation
  tstats <-  permuteG(data=pygmy, x="copies", y="values", plot=TRUE, n.rounds=1000)

# histogram of t-statistics from permuted data
  ggplot(as.data.frame(tstats), aes(x=tstats)) +
  geom_histogram(fill="white", color="steelblue") +
  geom_vline(xintercept = tstat, color="red")
  
# table of t-stats
  table(tstat < tstats)

# calculate p-value
  table(tstat < tstats)[2] / (table(tstat < tstats)[1] + table(tstat < tstats)[2])


#=========================================================================================#
# Jeremy’s bullshit
#=========================================================================================#

#  phenotype.means <- c(14,12,6)

# frequency of allele
#  p <- 0.9

# other measures
#  q <- 1 - p
#  m <- (phenotype.means[3] + phenotype.means[1]) / 2
#  a <- (phenotype.means[3] - m) ## or can use (m - phenotype.means[1])
#  d <- phenotype.means[2] - m

## Generate Phenotypes based on allele frequency and phenotypic means
## Assumed know genotype is phenotype without any error
## can change just change to some number and add a residual using rnorm
#  phenotype.sd <- 0

## generate to samples (i.e. draw from each parent to make genotype)
#  DF <- data.frame(Copies = rbinom(10000,2,p))
#  DF$Pheno <- ifelse(DF$Copies == 0,phenotype.means[1],NA)
#  DF$Pheno <- ifelse(DF$Copies == 1,phenotype.means[2],DF$Pheno)
#  DF$Pheno <- ifelse(DF$Copies == 2,phenotype.means[3],DF$Pheno)
#  DF <- DF[ ,c(2,1)]

## Check Genotype and allele Frequencies
#  table(DF$Copies)/10000
#  (sum(DF$Copies)/(nrow(DF)*2))






































































