rm(list=ls(all=TRUE)) # Clear all
setwd("//uoa.auckland.ac.nz/engdfs/home/mzha612/Desktop/PartIV/PartIVProject/Sample Distribution")

# distribution type
d_type <- 3
v_type <- 2  # Vessel type - arteriole, capillary, venule
c_type <- 1  # Characteristic type - diameter, length,

d_name <- c("norm","lognorm",'gamma')
v_name <- c("art","cap","ven")
c_name <- c("dia","len")

# sample size
n<-100000

mydata = read.csv("Digitised_Data.csv")  # read from first sheet



if(c_type == 1){ # Diameter
  num_values <- 20
  x <- mydata$Diameter_Freq[0:num_values]
  
  if(v_type == 1){ # Arteriole
    
    p_x <- mydata$Art_Dia_Freq[0:num_values]
    
  }else if(v_type == 2){ # Capillary
    p_x <- mydata$Cap_Dia_Freq[0:num_values]
    
  }else if(v_type == 3){ # Venule
    p_x <- mydata$Ven_Dia_Freq[0:num_values]
    
  }
  
}else if(c_type == 2){ # Length
  num_values <- 10
  x <- mydata$Length_Freq[0:num_values]
  

  if(v_type == 1){ # Arteriole
    
    p_x <- mydata$Art_Len_Freq[0:num_values]
    
  }else if(v_type == 2){ # Capillary
    p_x <- mydata$Cap_Len_Freq[0:num_values]
    
  }else if(v_type == 3){ # Venule
    p_x <- mydata$Ven_Len_Freq[0:num_values]
    
  }
  
}
filename <- paste(v_name[v_type],c_name[c_type],d_name[d_type],sep = "_")
filename <- paste(filename,".csv")
  
# Create variables



# Check
plot(x,p_x)

# Standardize the cumulative probabilities
cum_x <- cumsum(p_x); 
cum_x <- cum_x / cum_x[num_values]

# Compute sum of squared residuals to a fit
if (d_type == 1){
  f <- function(q) {
  res <- pnorm(x, q[1], q[2]) - cum_x
  sum(res * res)
  }
}else if(d_type == 2){
  f <- function(q) {
    res <- pnorm(log(x), q[1], q[2]) - cum_x
    sum(res * res)
  }
}else{
  f <- function(q) {
    res <- pgamma(x, q[1], q[2]) - cum_x
    sum(res * res)
  }
}

# Find the least squares fit
coeff <-(fit <- nlm(f, c(3, 2)))$estimate

# Plot the fit
if(d_type == 1){
  plot(x, cum_x)
  curve(pnorm(x, coeff[1], coeff[2]), add=TRUE)
  
}else if(d_type == 2){
  plot(log(x), cum_x)
  curve(pnorm(x, coeff[1], coeff[2]), add=TRUE)
}else{
  plot(x, cum_x)
  ccc <- curve(pgamma(x, coeff[1], coeff[2]), add=TRUE)
}

# Sample


if(d_type == 1){
  x<-rnorm(n, coeff[1], coeff[2])
  h <- hist(x,plot=FALSE)
  h$counts = h$counts/sum(h$counts)
  plot(h)
  line(x,dnorm(x, coeff[1], coeff[2]))
  # # Write to file
  # write.csv(x,file='x_norm.csv',row.names=FALSE)
  
}else if(d_type == 2){
  x<-exp(rnorm(n, coeff[1], coeff[2]))
  h <- hist(x,plot=FALSE)
  h$counts = h$counts/sum(h$counts)
  plot(h)
  line(x,dnorm(x, coeff[1], coeff[2]))
  
  # # Write to file
  # write.csv(x,file='x_lognorm.csv',row.names=FALSE)
  
}else {
  x<-rgamma(n, coeff[1], coeff[2])
  h <- hist(x,plot=FALSE)
  # breaks = seq(-1.25,56.25,2.5)
  h$counts = h$counts/sum(h$counts)
  plot(h)
  line(x,dgamma(x, coeff[1], coeff[2]))
  

}

# Write to file
write.csv(x,file=filename,row.names=FALSE)



