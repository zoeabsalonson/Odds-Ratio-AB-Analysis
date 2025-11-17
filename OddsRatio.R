rm(list=ls())
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Downloads/School/WWU/Winter 2025/MATH 445 - Computational Statistics/W25MATH445Project")
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)

#See https://stats.stackexchange.com/questions/588428/confidence-intervals-for-odds-ratio

#The OR_CI() function computes the odds ratio confidence interval.

#Input argument:
#n11: A numeric specifying the cell count for (1,1).
#n1d: A numeric specifying the total cell count for the first row.
#n21: A numeric specifying the cell count for (2,1).
#n2d: A numeric specifying the total cell count for the second row.
#dist: A character string specifying the distribution to be used. 
#"z" for standard normal, and "t" for t with Welch's adjustment
#adjust: A character string specifying the cell count adjustment. 
#"Agresti" for the independent-smoothed adjustment of Agresti.
#"Gart" for the Gart adjustment.
#"Woolf" for no adjustment, corresponding to the original Woolf logit interval
#alpha: A numeric related to confidence level. Note that the confidence level
#is given by 1-alpha. It is set to 0.05 by default.

OR_CI <- function(n11, n1d, n21, n2d, dist=c("z","t"), 
                  adjust=c("Agresti","Gart","Woolf"), alpha=0.05)
{
  dist <- match.arg(dist)
  adjust <- match.arg(adjust)
  #The second column counts
  n12 <- n1d - n11
  n22 <- n2d - n21
  #Total column counts
  nd1 <- n11 + n21
  nd2 <- n12 + n22
  #Grand total count
  N <- n11 + n12 + n21 + n22
  
  #Cell count adjustments
  if(adjust == "Gart")
  {
    n11a <- n11 + 0.5
    n12a <- n12 + 0.5
    n21a <- n21 + 0.5
    n22a <- n22 + 0.5
  }
  else if(adjust == "Agresti")
  {
    c11 <- 2*n1d*nd1/(N^2)
    c12 <- 2*n1d*nd2/(N^2)
    c21 <- 2*n2d*nd1/(N^2)
    c22 <- 2*n2d*nd2/(N^2)
    
    n11a <- n11 + c11
    n12a <- n12 + c12
    n21a <- n21 + c21
    n22a <- n22 + c22
  }
  else #adjust == "Woolf"
  {
    n11a <- n11
    n12a <- n12
    n21a <- n21
    n22a <- n22
  }
  
  #theta is the odds ratio (OR) estimate
  theta  <- (n11a*n22a)/(n12a*n21a)
  
  #Standard error
  se <- sqrt(1/n11a + 1/n12a + 1/n21a + 1/n22a)
  
  #Critical value calculation
  if(dist=="z")
  {
    crit <- qnorm(alpha/2, lower.tail=FALSE)
  }
  if(dist=="t")
  {
    Na <- n11a + n12a + n21a + n22a
    p11a <- n11a/Na
    p12a <- n12a/Na
    p21a <- n21a/Na
    p22a <- n22a/Na
    
    f1 <- 1/p11a + 1/p12a + 1/p21a + 1/p22a
    f3 <- 1/(p11a^3) + 1/(p12a^3) + 1/(p21a^3) + 1/(p22a^3)
    
    #Welch's degrees of freedom
    nu <- (2*Na*f1^2)/(f3 - f1^2)
    
    crit <- qt(alpha/2, df=nu, lower.tail=FALSE)
  }
  
  #Lower and upper bound of the confidence interval on the log scale
  loglower <- log(theta) - crit*se
  logupper <- log(theta) + crit*se
  
  #Lower and upper bound of the confidence inerval on the original scale
  lower <- exp(loglower)
  upper <- exp(logupper)
  
  #Handling NAs by assigning 0 and infinity for the lower and upper bound.
  lowerNAs <- which(is.na(lower))
  upperNAs <- which(is.na(upper))
  
  if(length(lowerNAs) > 0) lower[lowerNAs] <- 0
  if(length(upperNAs) > 0) upper[upperNAs] <- Inf
  
  results <- rbind(lower, upper)
  rownames(results) <- c("lower", "upper")
  
  return(results)
}

#The simulations matrix below gives all the different simulation scenarios to consider.
#n1d and n2d gives the total number of trials for the first and second population.
#p1 is the probability of success for the first population.
#OR is the odds ratio.
#alpha is the value so that 100(1-alpha)% represents the confidence level.

#Be sure to change these vectors below to match your research plan.
ORex <- c(1, 2, 5, 12, 50)
n1dex <- c(3, 10, 20, 30, 35, 75, 100, 600, 1000, 1005)
n2dex <- c(3, 10, 20, 30, 35, 75, 100, 600, 1000, 1005)
p1ex <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)
simulations <- expand.grid(p1=p1ex, n1d=n1dex, n2d=n2dex, 
                           OR=ORex, alpha=0.05)
simulations

#numsimcases gives the number of simulation cases.
numsimcases <- dim(simulations)[1]
numsimcases

#Number of Monte Carlo simulations 
#(10000 or more is recommended, but use 1000 for preliminary study)
MC <- 10000

#Initializing some vectors
empirical_CR_za <- c()
empirical_CR_ta <- c()
empirical_CR_zg <- c()
empirical_CR_tg <- c()
empirical_CR_zw <- c()
empirical_CR_tw <- c()
p2 <- c()
for(i in 1:numsimcases)
{
  print(paste("Running simulation ", i, sep=""))
  #Retrieving information about the i-th simulation scenario.
  curr_sim <- simulations[i,]
  curr_n1d <- curr_sim$n1d
  curr_n2d <- curr_sim$n2d
  curr_p1 <- curr_sim$p1
  curr_OR <- curr_sim$OR
  curr_p2 <- curr_p1/(curr_OR*(1-curr_p1)+curr_p1)
  p2[i] <- curr_p2
  curr_alpha <- curr_sim$alpha
  
  #Generating random n11's and n21's.
  n11s <- rbinom(n=MC, size=curr_n1d, prob=curr_p1)
  n21s <- rbinom(n=MC, size=curr_n2d, prob=curr_p2)
  
  #dist = "z" and adjust = "Agresti"
  case_za <- OR_CI(n11=n11s, n1d=curr_n1d, n21=n21s, n2d=curr_n2d, 
                 dist="z", adjust="Agresti", alpha=curr_alpha)

  #dist = "t" and adjust = "Agresti"
  case_ta <- OR_CI(n11=n11s, n1d=curr_n1d, n21=n21s, n2d=curr_n2d, 
                   dist="t", adjust="Agresti", alpha=curr_alpha)
  
  #dist = "z" and adjust = "Gart"
  case_zg <- OR_CI(n11=n11s, n1d=curr_n1d, n21=n21s, n2d=curr_n2d, 
                   dist="z", adjust="Gart", alpha=curr_alpha)
  
  #dist = "t" and adjust = "Gart"
  case_tg <- OR_CI(n11=n11s, n1d=curr_n1d, n21=n21s, n2d=curr_n2d, 
                   dist="t", adjust="Gart", alpha=curr_alpha)  
  
  #dist = "z" and adjust = "Woolf"
  case_zw <- OR_CI(n11=n11s, n1d=curr_n1d, n21=n21s, n2d=curr_n2d, 
                   dist="z", adjust="Woolf", alpha=curr_alpha)
  
  #dist = "t" and adjust = "Woolf"
  case_tw <- OR_CI(n11=n11s, n1d=curr_n1d, n21=n21s, n2d=curr_n2d, 
                   dist="t", adjust="Woolf", alpha=curr_alpha)    
  
  #Empirical coverage rate calculations
  empirical_CR_za[i] <- mean((case_za[1,] <= curr_OR)*(curr_OR <= case_za[2,])==1)
  empirical_CR_ta[i] <- mean((case_ta[1,] <= curr_OR)*(curr_OR <= case_ta[2,])==1)  
  empirical_CR_zg[i] <- mean((case_zg[1,] <= curr_OR)*(curr_OR <= case_zg[2,])==1)
  empirical_CR_tg[i] <- mean((case_tg[1,] <= curr_OR)*(curr_OR <= case_tg[2,])==1)   
  empirical_CR_zw[i] <- mean((case_zw[1,] <= curr_OR)*(curr_OR <= case_zw[2,])==1)
  empirical_CR_tw[i] <- mean((case_tw[1,] <= curr_OR)*(curr_OR <= case_tw[2,])==1)   
}

#Empirical coverage rates
empirical_CR_za
empirical_CR_ta
empirical_CR_zg
empirical_CR_tg
empirical_CR_zw
empirical_CR_tw

#All the results can be found here.
all_results <- cbind(simulations,
                 p2,
                 empirical_CR_za,
                 empirical_CR_ta,
                 empirical_CR_zg,
                 empirical_CR_tg,
                 empirical_CR_zw,
                 empirical_CR_tw                 
                 )
all_results

# Defining and testing values of interest

pdf("plots4.pdf", onefile = TRUE, width = 10, height = 8)  # Open PDF file

plot_list <- list() # Creates a plot list to hold all the plots (to print after)
index <- 1 # Initialize indices for plot list

for (i in c()) # VALUES MUST BE INSERTED, iterates through all values of interest for n1
{
  for (j in c()) # VALUES MUST BE INSERTED, iterates through all values of interest for OR
  {  
    for (k in c()) # VALUES MUST BE INSERTED, iterates through all values of interest for n2
    { 
      results_or <- filter(all_results, OR == j, n1d == i, n2d == k) # Filters results to give plots for values of interest
      
      plot_list[[index]] <- ggplot() + 
        geom_line(data=results_or, aes(x=p1, y=empirical_CR_za), color='red') + # Connects points with a line
          #geom_point(data=results_or, aes(x=p1, y=empirical_CR_za)) +  # Adds dots for each data point
        
        geom_line(data=results_or, aes(x=p1, y=empirical_CR_ta), color='orange') +
          #geom_point(data=results_or, aes(x=p1, y=empirical_CR_ta)) + 
        
        geom_line(data=results_or, aes(x=p1, y=empirical_CR_zg), color='yellow') +
          #geom_point(data=results_or, aes(x=p1, y=empirical_CR_zg)) + 
        
        geom_line(data=results_or, aes(x=p1, y=empirical_CR_tg), color='green') +
          #geom_point(data=results_or, aes(x=p1, y=empirical_CR_tg)) + 
        
        geom_line(data=results_or, aes(x=p1, y=empirical_CR_zw), color='blue') +
          #geom_point(data=results_or, aes(x=p1, y=empirical_CR_zw)) + 
        
        geom_line(data=results_or, aes(x=p1, y=empirical_CR_tw), color='purple') +
          #geom_point(data=results_or, aes(x=p1, y=empirical_CR_tw)) + 
        
        geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +  # Line at 0.95
        labs(y = "Empirical Coverage Rate", x = expression(p[1])) +
        labs(title = paste("OR =", j, ", n1d =", i, ", n2d =", k)) # Titles each plot with which values are looked at
        
        index <- index + 1  # Increase index after each plot creation
    }
  }
}

if (length(plot_list) > 0) # Make sure plot_list isn't empty
  {
  print(ggarrange(plotlist = plot_list, ncol = 3, nrow = 3)) # Print plots to pdf arranged 3x3 plots per page
  }

dev.off() # Close pdf

pdf("legend.pdf", onefile = TRUE, width = 10, height = 8)
library(ggplot2)

# Dummy data
legend_data <- data.frame(
  category = c("za", "ta", "zg", "tg", "zw", "tw"),
  x = 1,  # Dummy x-values
  y = 1   # Dummy y-values
)

# Create a legend using points (you can change to lines, bars, etc.)
legend_plot <- ggplot(legend_data, aes(x, y, color = category)) +
  geom_point(size = 5) +  # Use points, but can change to other geoms
  scale_color_manual(values = c("za" = "red", "ta" = "orange", "zg" = "yellow", "tg" = "green", "zw" = "blue", "tw" = "purple")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(color = "Confidence Interval") +
  theme_void() +  # Remove axes
  theme(legend.position = "left")  # Adjust legend placement

# Print the legend
legend_plot

dev.off()


#Save simulation results as a CSV file.
write.csv(all_results, file="all_results.csv", quote=FALSE, row.names=FALSE)