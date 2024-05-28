# Load necessary libraries
library(ggplot2)
library(dplyr)
library(fdANOVA)
library(zoo)
library(plotly)
library(reshape2)
library(data.table)
library(forecast)
library(fda)
library(fda.usc)
library(janitor)

# Load the smoothed functional data object from the file
smoothed_data <- readRDS("C:/Users/ACER/Desktop/studies/FDA/smoothed_fd_stock_prices.rds")

# Load sector information
sector_info <- read.csv("C:/Users/ACER/Desktop/studies/FDA/eda/unique_sp500_stocks_and_sectors.csv", header = TRUE)
sector_info <- sector_info %>%
  mutate(Stock = toupper(Stock))

# Extract stock names from the smoothed functional data
stock_names <- colnames(eval.fd(1:dim(smoothed_data$coefs)[1], smoothed_data))

# Validate that all stock headers in smoothed_data are listed in sector_info
missing_stocks <- setdiff(stock_names, sector_info$Stock)
if (length(missing_stocks) > 0) {
  stop("Error: Not all stock columns in the price data have corresponding sector info.")
}

# Generate group factors
group_factors <- factor(sector_info$Sector[match(stock_names, sector_info$Stock)])

# Define the time sequence
t.seq <- 1:126

# Function to perform Z-test
Ztwosample <- function(x, y, t.seq, alpha = 0.05) {
  if (class(x) != "fd") stop("X must be fd object")
  if (class(y) != "fd") stop("Y must be fd object")
  k <- length(t.seq)
  
  mu.x <- mean.fd(x)
  mu.y <- mean.fd(y)
  
  n <- dim(x$coef)[2]
  m <- dim(y$coef)[2]
  
  delta <- (mu.x - mu.y)
  delta.t <- eval.fd(t.seq, delta)
  
  z.x <- center.fd(x)
  z.y <- center.fd(y)
  
  z.x.t <- eval.fd(t.seq, z.x)
  z.y.t <- eval.fd(t.seq, z.y)
  z.t <- cbind(z.x.t, z.y.t)
  
  if (n > k) {
    Sigma <- (t(z.t) %*% z.t) / (n + m - 2)
  } else {
    Sigma <- (z.t %*% t(z.t)) / (n + m - 2)
  }
  
  gamma.t <- diag(Sigma)
  Zpointwise <- sqrt((n * m) / (n + m)) * delta.t / sqrt(gamma.t)
  
  crit <- qt(1 - alpha / 2, n - 2)
  crit.val <- rep(crit, k)
  params <- list(critical.value = crit)
  
  mx <- max(cbind(Zpointwise, crit.val))
  mn <- min(cbind(Zpointwise, -crit.val))
  
  plot(t.seq, Zpointwise, type = "l", xlab = 'Time', ylab = "Z statistics",
       main = "Two samples t-test", ylim = c(mn - 0.5, mx + 0.5))
  lines(t.seq, crit.val, lty = 2, lwd = 2, col = "blue")
  lines(t.seq, -crit.val, lty = 2, lwd = 2, col = "blue")
  
  return(list(statistics.pointwise = Zpointwise,
              params = params))
}

# Function to perform L2-norm based two sample test
L2.stat.twosample <- function(x, y, t.seq, alpha = 0.05, method = 1, replications = 100) {
  if (class(x) != "fd") stop("X must be fd object")
  if (class(y) != "fd") stop("Y must be fd object")
  
  mu.x <- mean.fd(x)
  mu.y <- mean.fd(y)
  
  n <- dim(x$coefs)[2]
  m <- dim(y$coefs)[2]
  
  k <- length(t.seq)
  
  cn <- (n * m) / (n + m)
  delta <- (mu.x - mu.y)
  delta.t <- eval.fd(t.seq, delta)
  
  z.x <- center.fd(x)
  z.y <- center.fd(y)
  
  z.x.t <- eval.fd(t.seq, z.x)
  z.y.t <- eval.fd(t.seq, z.y)
  z.t <- cbind(z.x.t, z.y.t)
  
  if (n > k | m > k) {
    Sigma <- (t(z.t) %*% z.t) / (n + m - 2)
  } else {
    Sigma <- (z.t %*% t(z.t)) / (n + m - 2)
  }
  
  A <- sum(diag(Sigma))
  B <- sum(diag(Sigma %*% Sigma))
  
  L2stat <- cn * t(delta.t) %*% delta.t
  L2stat <- L2stat[1]
  
  btL2stat <- numeric(replications)
  
  if (method == 1) { # naive method
    A2 <- A^2
    B2 <- B
    alp <- B2 / A
    df <- A2 / B2
    pvalue <- 1 - pchisq(L2stat / alp, df)
    params <- list(alpha = alp, df = df)
  } else if (method == 2) {  # bootstrapping method
    for (i in 1:replications) {
      rep1 <- sample(1:n, n, replace = TRUE)
      xstar.coefs <- x$coefs[, rep1]
      xstar.names <- x$fdnames
      xstar.names$reps <- rep1
      xstar.fd <- fd(xstar.coefs, x$basis, xstar.names)
      
      rep2 <- sample(1:m, m, replace = TRUE)
      ystar.coefs <- y$coefs[, rep2]
      ystar.names <- y$fdnames
      ystar.names$reps <- rep2
      ystar.fd <- fd(ystar.coefs, y$basis, ystar.names)
      
      mu.x.star <- mean.fd(xstar.fd)
      mu.y.star <- mean.fd(ystar.fd)
      delta.star <- (mu.x.star - mu.y.star)
      delta.star.t <- eval.fd(t.seq, delta.star)
      
      btmu <- apply(delta.star.t, 1, mean) - delta.t
      btL2stat[i] <- cn * t(btmu) %*% btmu
    }
    pvalue <- mean(btL2stat >= L2stat)
    params <- list(boot.stat = btL2stat)
  }
  return(list(statistics = L2stat, pvalue = pvalue, params = params))
}

# Function to perform F-type two sample test
F.stat.twosample <- function(x, y, t.seq, alpha = 0.05, method = 1, replications = 100) {
  if (class(x) != "fd") stop("X must be fd object")
  if (class(y) != "fd") stop("Y must be fd object")
  
  mu.x <- mean.fd(x)
  mu.y <- mean.fd(y)
  
  n <- dim(x$coefs)[2]
  m <- dim(y$coefs)[2]
  
  k <- length(t.seq)
  
  cn <- (n * m) / (n + m)
  delta <- (mu.x - mu.y)
  delta.t <- eval.fd(t.seq, delta)
  
  z.x <- center.fd(x)
  z.y <- center.fd(y)
  
  z.x.t <- eval.fd(t.seq, z.x)
  z.y.t <- eval.fd(t.seq, z.y)
  z.t <- cbind(z.x.t, z.y.t)
  
  if (n > k | m > k) {
    Sigma <- (t(z.t) %*% z.t) / (n + m - 2)
  } else {
    Sigma <- (z.t %*% t(z.t)) / (n + m - 2)
  }
  
  A <- sum(diag(Sigma))
  B <- sum(diag(Sigma %*% Sigma))
  
  Fstat <- (cn * t(delta.t) %*% delta.t) / A
  Fstat <- Fstat[1]
  
  btFstat <- numeric(replications)
  
  if (method == 1) { # naive method
    kappa <- A^2 / B
    pvalue <- 1 - pf(Fstat, kappa, (n + m - 2) * kappa)
    params <- list(df1 = kappa, df2 = (n + m - 2) * kappa)
  } else if (method == 2) {  # bootstrapping method
    for (i in 1:replications) {
      rep1 <- sample(1:n, n, replace = TRUE)
      xstar.coefs <- x$coefs[, rep1]
      xstar.names <- x$fdnames
      xstar.names$reps <- rep1
      xstar.fd <- fd(xstar.coefs, x$basis, xstar.names)
      
      rep2 <- sample(1:m, m, replace = TRUE)
      ystar.coefs <- y$coefs[, rep2]
      ystar.names <- y$fdnames
      ystar.names$reps <- rep2
      ystar.fd <- fd(ystar.coefs, y$basis, ystar.names)
      
      mu.x.star <- mean.fd(xstar.fd)
      mu.y.star <- mean.fd(ystar.fd)
      delta.star <- (mu.x.star - mu.y.star)
      delta.star.t <- eval.fd(t.seq, delta.star)
      
      bt.z.x <- center.fd(xstar.fd)
      bt.z.y <- center.fd(ystar.fd)
      
      bt.z.x.t <- eval.fd(t.seq, bt.z.x)
      bt.z.y.t <- eval.fd(t.seq, bt.z.y)
      z.star.t <- cbind(bt.z.x.t, bt.z.y.t)
      
      if (n > k | m > k) {
        btSigma <- (t(z.star.t) %*% z.star.t) / (n + m - 2)
      } else {
        btSigma <- (z.star.t %*% t(z.star.t)) / (n + m - 2)
      }
      
      btmu <- apply(delta.star.t, 1, mean) - delta.t
      btFstat[i] <- (cn * t(btmu) %*% btmu) / sum(diag(btSigma))
    }
    pvalue <- mean(btFstat >= Fstat)
    params <- list(btFstat = btFstat)
  }
  return(list(statistics = Fstat, pvalue = pvalue, params = params))
}

# Define the F.stat function for one sample test
F.stat <- function(fdobj, t.seq, mu0, alpha = 0.05) {
  if (class(fdobj) != "fd") stop("fdobj must be fd object")
  if (class(mu0) != "fd") stop("mu0 must be fd object")
  
  k <- length(t.seq)
  n <- dim(fdobj$coefs)[2]
  
  delta <- (mean.fd(fdobj) - mu0)
  delta.t <- eval.fd(t.seq, delta)
  
  z <- center.fd(fdobj)
  z.t <- eval.fd(t.seq, z)
  
  if (n > k) {
    Sigma <- (t(z.t) %*% z.t) / (n - 1)
  } else {
    Sigma <- (z.t %*% t(z.t)) / (n - 1)
  }
  
  A <- sum(diag(Sigma))
  B <- sum(diag(Sigma %*% Sigma))
  
  Fstat <- (n * t(delta.t) %*% delta.t) / A
  Fstat <- Fstat[1]
  
  if (A^2 / B > 0) {
    kappa <- A^2 / B
    pvalue <- 1 - pf(Fstat, kappa, (n - 1) * kappa)
    params <- list(df1 = kappa, df2 = (n - 1) * kappa)
  } else {
    pvalue <- NA
    params <- NULL
  }
  
  return(list(statistics = Fstat, pvalue = pvalue, params = params))
}

# Define the L2.stat function for one sample test
L2.stat <- function(fdobj, t.seq, mu0, alpha = 0.05, method = 1, replications = 100) {
  if (class(fdobj) != "fd") stop("fdobj must be fd object")
  if (class(mu0) != "fd") stop("mu0 must be fd object")
  
  k <- length(t.seq)
  n <- dim(fdobj$coefs)[2]
  
  delta <- (mean.fd(fdobj) - mu0)
  delta.t <- eval.fd(t.seq, delta)
  
  z <- center.fd(fdobj)
  z.t <- eval.fd(t.seq, z)
  
  if (n > k) {
    Sigma <- (t(z.t) %*% z.t) / (n - 1)
  } else {
    Sigma <- (z.t %*% t(z.t)) / (n - 1)
  }
  
  A <- sum(diag(Sigma))
  B <- sum(diag(Sigma %*% Sigma))
  
  L2stat <- n * sum(delta.t^2)
  
  if (method == 1) { # naive method
    A2 <- A^2
    B2 <- B
    alp <- B2 / A
    df <- A2 / B2
    pvalue <- 1 - pchisq(L2stat / alp, df)
    params <- list(alpha = alp, df = df)
  } else if (method == 2) { # bootstrapping method
    btL2stat <- numeric(replications)
    for (i in 1:replications) {
      rep <- sample(1:n, n, replace = TRUE)
      xstar.coefs <- fdobj$coefs[, rep]
      xstar.names <- fdobj$fdnames
      xstar.names$reps <- rep
      xstar.fd <- fd(xstar.coefs, fdobj$basis, xstar.names)
      
      mu.star <- mean.fd(xstar.fd)
      delta.star <- (mu.star - mu0)
      delta.star.t <- eval.fd(t.seq, delta.star)
      
      btmu <- apply(delta.star.t, 1, mean) - delta.t
      btL2stat[i] <- n * sum(btmu^2)
    }
    pvalue <- mean(btL2stat >= L2stat)
    params <- list(boot.stat = btL2stat)
  } else {
    stop("method should be either 1 or 2")
  }
  
  return(list(statistics = L2stat, pvalue = pvalue, params = params))
}

# Define the sectors
sectors <- c("Communication Services", "Consumer Discretionary", "Consumer Staples", "Energy", 
             "Financials", "Health Care", "Industrials", "Information Technology", "Materials", 
             "Real Estate", "Utilities")

# Perform tests for each pair of sectors
for (i in 1:(length(sectors) - 1)) {
  for (j in (i + 1):length(sectors)) {
    sector_1_name <- sectors[i]
    sector_2_name <- sectors[j]
    
    sector_1 <- smoothed_data[which(group_factors == sector_1_name), ]
    sector_2 <- smoothed_data[which(group_factors == sector_2_name), ]
    
    cat("Testing sectors:", sector_1_name, "vs", sector_2_name, "\n")
    
    # Perform Z two-sample test
    Ztwosample_result <- tryCatch({
      Ztwosample(sector_1, sector_2, t.seq)
    }, error = function(e) {
      cat("Error in Ztwosample:", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(Ztwosample_result)) {
      cat("Z-twosample test: Statistically Significant =", any(Ztwosample_result$statistics.pointwise > Ztwosample_result$params$critical.value), "\n")
    }
    
    # Perform two-sample F-stat test
    F_stat_twosample_result <- tryCatch({
      F.stat.twosample(sector_1, sector_2, t.seq)
    }, error = function(e) {
      cat("Error in F.stat.twosample:", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(F_stat_twosample_result)) {
      cat("F-stat twosample test: p-value =", F_stat_twosample_result$pvalue, "Statistically Significant =", F_stat_twosample_result$pvalue < 0.05, "\n")
    }
    
    # Perform two-sample L2-stat test
    L2_stat_twosample_result <- tryCatch({
      L2.stat.twosample(sector_1, sector_2, t.seq)
    }, error = function(e) {
      cat("Error in L2.stat.twosample:", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(L2_stat_twosample_result)) {
      cat("L2-stat twosample test: p-value =", L2_stat_twosample_result$pvalue, "Statistically Significant =", L2_stat_twosample_result$pvalue < 0.05, "\n")
    }
    
    cat("\n")
  }
}

# Hypothesis testing function using smoothed Data
fANOVA.pointwise <- function(fdobj, groups, t.seq, alpha = 0.05) {
  data_matrix <- eval.fd(t.seq, fdobj)
  n <- nrow(data_matrix)
  pvals <- numeric(n)
  mean.p <- matrix(NA, ncol = length(levels(groups)), nrow = n)
  Tukey.posthoc <- vector("list", n)
  
  for (i in 1:n) {
    dt <- data.frame(values = data_matrix[i, ], group = groups, stringsAsFactors = FALSE)
    
    # Skip if any group has less than two observations
    if (any(table(dt$group) < 2)) {
      pvals[i] <- NA
      mean.p[i, ] <- NA
      Tukey.posthoc[[i]] <- NA
      next
    }
    
    av <- aov(values ~ group, data = dt)
    pvals[i] <- summary(av)[[1]][["Pr(>F)"]][1]
    mean.p[i, ] <- tapply(dt$values, dt$group, mean, na.rm = TRUE)
    
    if (!is.null(av) && "group" %in% names(TukeyHSD(av))) {
      Tukey.posthoc[[i]] <- TukeyHSD(av)$group[, "p adj"]
    } else {
      Tukey.posthoc[[i]] <- rep(NA, length(levels(groups)))
    }
  }
  
  overall_mean <- rowMeans(data_matrix, na.rm = TRUE)
  
  par(mfrow = c(2, 1))
  plot(t.seq, pvals, type = "l", main = "Pointwise ANOVA p-values", xlab = "Time", ylab = "p-value", ylim = c(0, 1))
  abline(h = alpha, col = "red", lty = 2)
  plot(t.seq, overall_mean, type = "l", main = "Overall Mean per Time Point", xlab = "Time", ylab = "Mean Value")
  par(mfrow = c(1, 1))
  
  return(list(p.values = pvals, TukeyHSD = Tukey.posthoc, group.means = mean.p, overall.mean = overall_mean, data_matrix = data_matrix))
}

# Run the hypothesis testing function
results <- fANOVA.pointwise(smoothed_data, group_factors, t.seq)

# Access the returned data_matrix from the results
data_matrix <- results$data_matrix

# Print the overall p-values
cat("P-values:\n")
print(results$p.values)

# Print the overall mean values
cat("Overall Mean per Time Point:\n")
print(results$overall.mean)

# Print the group means for each time point
cat("Group Means:\n")
print(results$group.means)

# Print the Tukey HSD post hoc test results
cat("Tukey HSD Results:\n")
print(results$TukeyHSD)

# Create significant matrix for heatmap visualization
num_pairs <- length(levels(group_factors)) * (length(levels(group_factors)) - 1) / 2
sig_matrix <- matrix(0, nrow = length(t.seq), ncol = num_pairs)

# Check the structure of Tukey HSD results
pair_names <- NULL
if (!is.null(results$TukeyHSD) && length(results$TukeyHSD) > 0 && is.list(results$TukeyHSD[[1]])) {
  pair_names <- rownames(results$TukeyHSD[[1]])
  colnames(sig_matrix) <- pair_names
  
  for (i in 1:length(t.seq)) {
    if (!is.na(results$p.values[i]) && results$p.values[i] < 0.05) {
      tukey_results <- results$TukeyHSD[[i]]
      if (!is.null(tukey_results) && length(tukey_results) > 1) {
        significant_pairs <- rownames(tukey_results[tukey_results < 0.05])
        for (pair in significant_pairs) {
          sig_matrix[i, pair] <- 1
        }
      }
    }
  }
}

# Convert to a dataframe for easier manipulation
sig_df <- as.data.frame(sig_matrix)
sig_df$time_point <- t.seq

# Summarize the frequency of significant differences for each pair
sig_summary <- colSums(sig_df[, -ncol(sig_df)])
sig_summary <- sort(sig_summary, decreasing = TRUE)

# Create direction matrix
direction_matrix <- matrix(0, nrow = length(t.seq), ncol = num_pairs)

if (!is.null(results$TukeyHSD) && length(results$TukeyHSD) > 0 && is.list(results$TukeyHSD[[1]])) {
  pair_names <- rownames(results$TukeyHSD[[1]])
  colnames(direction_matrix) <- pair_names
  
  for (i in 1:length(t.seq)) {
    if (!is.na(results$p.values[i]) && results$p.values[i] < 0.05) {
      tukey_results <- results$TukeyHSD[[i]]
      if (!is.null(tukey_results) && length(tukey_results) > 1) {
        significant_pairs <- rownames(tukey_results[tukey_results$p.adj < 0.05, ])
        for (pair in significant_pairs) {
          diff <- tukey_results[pair, "diff"]
          direction_matrix[i, pair] <- sign(diff)
        }
      }
    }
  }
}

# Convert to a dataframe for easier manipulation
direction_df <- as.data.frame(direction_matrix)
colnames(direction_df) <- pair_names
direction_df$time_point <- t.seq

# Plot sector performance over time
sector_performance <- data.frame(time_point = t.seq)

# Calculate the mean performance for each sector at each time point
for (sector in levels(group_factors)) {
  sector_indices <- which(group_factors == sector)
  sector_means <- rowMeans(data_matrix[, sector_indices], na.rm = TRUE)
  sector_performance[[sector]] <- sector_means
}

# Melt the data for plotting
sector_performance_melted <- melt(sector_performance, id.vars = "time_point", variable.name = "Sector", value.name = "MeanPerformance")

# Plot sector performance over time
p <- ggplot(sector_performance_melted, aes(x = time_point, y = MeanPerformance, color = Sector)) +
  geom_line(size = 1.2, alpha = 0.8) +
  labs(title = "Sector Performance Over Time", subtitle = "Mean Performance of Each Sector",
       x = "Time Point", y = "Mean Performance", color = "Sector") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~Sector, scales = "free_y")

# Convert to an interactive plotly object
p <- ggplotly(p)
p
