# Load libraries that are required
library(fda)
library(fda.usc)
library(foreign)
library(fdaoutlier)

# Load the dataset
stock_prices <- read.csv("C:/Users/ACER/Desktop/studies/FDA/eda/normalized_sp500_daily_closing_prices_2nd_half_2023.csv", header = TRUE)

# Ensure dates are correctly formatted
stock_prices$Date <- as.Date(stock_prices$Date)

# Prepare data for functional data analysis
argvals <- 1:nrow(stock_prices)  # date: days
rangeval <- c(min(argvals), max(argvals))
stock_data_matrix <- as.matrix(stock_prices[, -1])  # Exclude the Date column

# LOOCV function to assess model fit, focusing on predictive performance
loocv_error <- function(nbasis, stock_data_matrix, argvals) {
  total_error <- 0
  n <- nrow(stock_data_matrix)  # number of rows
  
  for (i in 1:n) {
    train_indices <- setdiff(1:n, i)
    
    # Create basis object for the entire dataset range
    basis_obj <- create.bspline.basis(rangeval = c(min(argvals), max(argvals)), nbasis = nbasis)
    
    # Create functional data object for the training subset
    if(length(train_indices) > 1) {  # Check to avoid creating fd object with a single observation
      fd_train <- Data2fd(argvals = argvals[train_indices], y = stock_data_matrix[train_indices, ], basisobj = basis_obj)
      
      # Compute the mean of functional data
      if(is.fd(fd_train)) {  # Check if fd_train is indeed a functional data object
        mean_fd <- mean.fd(fd_train)
        predicted <- eval.fd(argvals[i], mean_fd)
        
        # Compute error for the excluded observation
        actual <- stock_data_matrix[i, ]
        error <- sum((actual - predicted)^2)
        
        total_error <- total_error + error
      } 
    }
  }
  
  # Calculate and return mean error
  mean_error <- total_error / n
  return(mean_error)
}

# Test the LOOCV function with a range of basis functions
nbasis_values <- c(10, 15, 20, 25, 30, 40)
errors <- sapply(nbasis_values, function(nb) loocv_error(nb, stock_data_matrix, argvals))

# Test the LOOCV function with a range of basis functions
errors <- unlist(lapply(nbasis_values, function(nb) loocv_error(nb, stock_data_matrix, argvals)))


# Plotting the errors to find the optimal number of basis functions 
plot(nbasis_values, errors, type = "b", xlab = "Number of Basis Functions", ylab = "Average LOOCV Error")

# Explore different `nbasis` values by plotting functional data objects
par(mfrow=c(3,2))
for(nb in nbasis_values) {
  basis_obj <- create.bspline.basis(rangeval, nb)
  fd_stock_prices <- Data2fd(argvals, stock_data_matrix, basisobj = basis_obj)
  plot(fd_stock_prices, main = paste("nbasis =", nb))
}

# Choosed nbasis using the graph
optimal_nbasis <- 15

# Function to compute and plot plot GCV, df and SSE values with different lambda values
par(mfrow=c(2,2)) 
compute_smoothed_data <- function(lambdas, basis_obj, argvals, stock_data_matrix) {
  gcv <- rep(0, length(lambdas))
  df <- rep(0, length(lambdas))
  sse <- rep(0, length(lambdas))
  
  for (i in seq_along(lambdas)) {
    fdPar_obj <- fdPar(basis_obj, Lfdobj = int2Lfd(2), lambda = lambdas[i])
    smoothed_data <- smooth.basis(argvals, stock_data_matrix, fdPar_obj)
    gcv[i] <- sum(smoothed_data$gcv)
    df[i] <- smoothed_data$df
    sse[i] <- smoothed_data$SSE
  }
  
  return(list(gcv = gcv, df = df, sse = sse))
}

# Integrate lambda and plot GCV, df, and SSE values
lambda_values <- 10^seq(-10, -1, length.out = 21) # from 1e-10 to 0.1
smoothed_data_results <- compute_smoothed_data(lambda_values, create.bspline.basis(rangeval, 15), argvals, stock_data_matrix)

# Plotting GCV, df, and SSE values
plot(log10(lambda_values), smoothed_data_results$df, type = 'b', xlab = 'log lambda', ylab = 'Degrees of Freedom (df)')
plot(log10(lambda_values), smoothed_data_results$sse, type = 'b', xlab = 'log lambda', ylab = 'Sum of Squared Errors (SSE)')
plot(log10(lambda_values), smoothed_data_results$gcv, type = 'b', xlab = 'log lambda', ylab = 'Generalized Cross-Validation (GCV)')

# Choosed lambda
optimal_lambda = 0.001

# Bspline parameters
optimal_norder = 4
basis_obj <- create.bspline.basis(rangeval, nbasis = optimal_nbasis, norder = optimal_norder)

# Perform Functional Principal Component Analysis (FPCA)
fd_stock_prices <- Data2fd(argvals = argvals, y = stock_data_matrix, basisobj = basis_obj)
pca_result <- pca.fd(fd_stock_prices)

# Plot cumulative percentage of variance explained by principal components
plot(cumsum(pca_result$values) / sum(pca_result$values),type = 'b', xlab = 'Number of components', ylab = 'Cumulative Variance Explained')
abline(h=0.95, col = "red",lty=2)


# Define a fdPar object for smoothing
fdParObj <- fdPar(basis_obj, Lfdobj = int2Lfd(2), lambda = optimal_lambda)

# Smooth the data using the fdparobj
smooth_list <- smooth.basis(argvals, stock_data_matrix, fdParObj)
smoothed_curves <- eval.fd(argvals, smooth_list$fd)

# Plot original and smoothed data
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
matplot(argvals, stock_data_matrix, type = 'l', lty = 1, col = 1:ncol(stock_data_matrix), ylab = "Stock Price", main = "Original Data", xlab = "Time")
matplot(argvals, smoothed_curves, type = 'l', lty = 1, col = 1:ncol(smoothed_curves), ylab = "Smoothed Stock Price", main = "Smoothed Data", xlab = "Time")
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)  # Resetting to default layout


# Create fd object for depth analysis
fdataobj_stock_prices <- fdata(smoothed_curves, argvals)

# Plot to visually check the smoothed functional data
plot(fdataobj_stock_prices)

# Fraiman-Muniz Depth
out_FM_stock <- depth.FM(fdataobj_stock_prices, trim = 0.1, draw = TRUE)
# Modal Depth
out_mode_stock <- depth.mode(fdataobj_stock_prices, trim = 0.1, draw = TRUE)
# Random Projection Depth
out_RP_stock <- depth.RP(fdataobj_stock_prices, trim = 0.1, draw = TRUE)

# Calculate band depth
stock_bd <- band_depth(dt = t(smoothed_curves))
names(stock_bd) <- colnames(smoothed_curves)

# Plot band depth
plot(stock_bd, type = "l", main = "Band Depth for Stock Prices", xlab = "Index", ylab = "Depth")

# Calculate modified band depth
stock_mbd <- modified_band_depth(t(smoothed_curves))
names(stock_mbd) <- colnames(smoothed_curves)

# Plot modified band depth
plot(stock_mbd, type = "l", main = "Modified Band Depth for Stock Prices", xlab = "Index", ylab = "Depth")

# Perform functional boxplots to detect outliers with band depth
fbplot_obj_bd <- functional_boxplot(t(smoothed_curves), depth_method = "bd")
outliers_bd <- fbplot_obj_bd$outliers

# Perform functional boxplots to detect outliers with modified band depth
fbplot_obj_mbd <- functional_boxplot(t(smoothed_curves), depth_method = "mbd")
outliers_mbd <- fbplot_obj_mbd$outliers

# Use muod function for outlier detection
outliers_muod <- muod(t(smoothed_curves), cut_method = c("boxplot"))

# Print outliers detected by MUOD for each component
cat("Shape outliers: ", outliers_muod$outliers$shape, "\n")
cat("Amplitude outliers: ", outliers_muod$outliers$amplitude, "\n")
cat("Magnitude outliers: ", outliers_muod$outliers$magnitude, "\n")

# Visualize functional boxplots
fbplot(t(smoothed_curves), method = "BD2")
fbplot(t(smoothed_curves), method = "MBD")

# Combine all outlier indices
all_outlier_columns_ind <- unique(c(outliers_muod$outliers$shape, outliers_muod$outliers$amplitude, outliers_muod$outliers$magnitude))
outlier_cols <- colnames(smoothed_curves[, all_outlier_columns_ind])
cat("Number of outliers stocks: ", length(outlier_cols), "\n")

# Remove outliers from smoothed_curves by columns
cleaned_smoothed_curves <- smoothed_curves[, -all_outlier_columns_ind]

# Check the dimensions of the cleaned data
dim(cleaned_smoothed_curves)

# Plot with original data and outliers highlighted
matplot(argvals, smoothed_curves, type='l', col=ifelse(colnames(smoothed_curves) %in% outlier_cols, 'red', 'gray'), lty=1, main="Original Data with Outliers Highlighted", xlab="Time", ylab="Stock Price")
legend("topright", legend=c("Outlier", "Non-Outlier"), col=c("red", "gray"), lty=1, cex=0.8)

# Calculate central tendency measures
fd_stock_prices <- Data2fd(argvals, y = cleaned_smoothed_curves, basisobj = basis_obj)

# Calculate the mean function
fd_mean <- mean.fd(fd_stock_prices)

data_matrix <- eval.fd(argvals, fd_stock_prices)

# Calculate the median for each time point
pointwise_median <- apply(data_matrix, 1, median)

# Correctly evaluate the mean function over 'argvals'
mean_values <- eval.fd(argvals, fd_mean)

# Calculate the standard deviation for each time point
sd_values <- apply(data_matrix, 1, sd)

par(mfrow = c(1, 2)) # Split the plot area into 1 row and 2 columns 

# Central Tendency Measures Plot
plot(argvals, mean_values, type = 'l', lwd = 2, col = 'blue', ylim = range(c(mean_values - sd_values, mean_values + sd_values, pointwise_median)), xlab = 'Time', ylab = 'Values', main = 'Mean, Median, and SD of Functional Data')
lines(argvals, pointwise_median, col = 'grey', lwd = 2) # Plot Median
lines(argvals, mean_values + sd_values, col = "red", lty = 2) # +1 SD
lines(argvals, mean_values - sd_values, col = "green", lty = 2) # -1 SD
legend("topright", legend = c("Mean", "Median", "+1 SD", "-1 SD"), col = c("blue", "grey", "red", "green"), lty = 1:2, lwd = 2)

# Dispersion Measures Plot
# For the dispersion measures plot, let's focus on plotting the standard deviation itself to illustrate variability
plot(argvals, sd_values, type = 'l', col = 'purple', lwd = 2, xlab = 'Time', ylab = 'Standard Deviation', main = 'Variability Over Time')

# Save the cleaned_smoothed_curves as an RDS file
saveRDS(fd_stock_prices, file = "C:/Users/ACER/Desktop/studies/FDA/smoothed_fd_stock_prices.rds")
























## Canadian Weather PCA

# Precipitation data
logprecav <- CanadianWeather$dailyAv[, , 'log10precip']
yearRng <- c(0, 365)
daybasis <- create.fourier.basis(yearRng, 365)
Lcoef <- c(0, (2 * pi / diff(yearRng))^2, 0)
harmaccelLfd <- vec2Lfd(Lcoef, yearRng)

lambda <- 1e6
fdParobj <- fdPar(daybasis, harmaccelLfd, lambda)
logprec.fd <- smooth.basis(day.5, logprecav, fdParobj)$fd

nharm <- 4
pcalist <- pca.fd(logprec.fd, nharm, centerfns = FALSE)
plot(pcalist)
plot(pcalist$harmonics)

# plot(pcalist)
# legend("topright", legend = paste("PC", 1:nharm), col = 1:nharm, lty = 1:nharm, cex = 0.8)
# 
# # Plot the harmonics
# plot(pcalist$harmonics)
# legend("topright", legend = paste("Harmonic", 1:nharm), col = 1:nharm, lty = 1:nharm, cex = 0.8)


#plotscores(pcalist, loc = 5)  ###*** check latency and why code get stuck

# Function to restore the original curves using PCA scores and harmonics
restore_curves <- function(pcalist, n_components) {
  fd_pca_list <- list()
  for (i in 1:n_components) {
    fd_pca <- mean.fd(logprec.fd)
    for (j in 1:i) {
      fd_pca <- fd_pca + pcalist$scores[i, j] * pcalist$harmonics[j]
    }
    fd_pca_list[[i]] <- fd_pca
  }
  return(fd_pca_list)
}

# Restore original curves using PCA
fd_pca_list <- restore_curves(pcalist, nharm)

print(length(fd_pca_list))

# Plot restored curves
opar <- par(mfrow=c(2,2), ask = TRUE)
for (i in 1:length(fd_pca_list)) {
  plot(fd_pca_list[[i]], ylim=c(-1, 1), ylab = paste(i, "PC"))
  lines(logprec.fd[i], col = 2)
}
par(opar)

#### Rotation
varmx <- varmx.pca.fd(pcalist)
plot(varmx)
plot(varmx$harmonics)
#plotscores(varmx, loc = 5)   ###*** check latency and why code get stuck

# Restore original curves using rotated PCA
fd_vrm_list <- restore_curves(varmx, nharm)

# Check the length of fd_vrm_list
print(length(fd_vrm_list))

# Plot restored curves using rotated PCA
opar <- par(mfrow=c(2,2), ask = TRUE)
for (i in 1:length(fd_vrm_list)) {
  plot(fd_vrm_list[[i]], ylim=c(-1, 1), ylab = paste(i, "PC"))
  lines(logprec.fd[i], col = 2)
}
par(opar)










