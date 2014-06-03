par(mfrow = c(1, 2))
#metrics <- c("bsor", "bsim", "bnes")
metrics <- c("bsim")
for (k in unique(metrics)){
  for(m in unique(c("tax"))){
    if (m == "tax"){
      dissim <- comm.dissim
      partial.mantels.comp <- taxonomic.partial.mantels.comp
    }
    else {
      dissim <- fun.dissim
      partial.mantels.comp <- functional.partial.mantels.comp
    }
    
    
    if(k == "bnes" && m == "fun"){
      k <- "bsne"
    }
    else {
      k <- k
    }
    
    
    for (j in unique(comm$seasn)){
      ylabel <- paste("Community Dissimilarity Residuals:", k, sep = " ")
      plot(x = c(-500, 500), y = c(-0.25, 0.25), ylab = ylabel, xlab = "Elevational Distance Residuals", type = "n",
           main = paste())
      
      for (i in unique(comm$year)){
        c.dist <- dissim[[paste(i, j, k, sep = "")]]
        e.dist <- ele.dissim[[paste(i, j, sep = "")]]
        s.dist <- spa.dissim[[paste(i, j, sep = "")]]
        
        residual_vectors <- cbind(c.dist, e.dist, s.dist)
        resid_vecs <- as.data.frame(residual_vectors)
        
        comm_residualsG <- resid(lm(c.dist ~ e.dist, data = resid_vecs))
        eco_residuals <- resid(lm(s.dist ~ e.dist, data = resid_vecs))
        comm_residualsE <- resid(lm(c.dist ~ s.dist, data = resid_vecs))
        geo_residuals <- resid(lm(e.dist ~ s.dist, data = resid_vecs))
        
        ind <- substring(k, 2, 4)
        
        
        
        
        
        elep <- partial.mantels.comp[partial.mantels.comp$index == ind & 
                                       partial.mantels.comp$gradient == "ele" & 
                                       partial.mantels.comp$year == i & 
                                       partial.mantels.comp$seasn == j, "pval3"]
        
        
        
        
        linecol <- ifelse(elep < 0.05, "red", "blue")
        eline <- ifelse(elep < 0.05, 1, 3)
        abline(lm(comm_residualsE~-1+geo_residuals), col = linecol, lty = eline)
        
        
      }
    }
  }
}




















partial_plot2F <- function(i, j, k){
  c.dist <- fun.dissim[[paste(i, j, k, sep = "")]]
  e.dist <- ele.dissim[[paste(i, j, sep = "")]]
  s.dist <- spa.dissim[[paste(i, j, sep = "")]]
  
  residual_vectors <- cbind(c.dist, e.dist, s.dist)
  resid_vecs <- as.data.frame(residual_vectors)
  
  comm_residualsG <- resid(lm(c.dist ~ e.dist, data = resid_vecs))
  eco_residuals <- resid(lm(s.dist ~ e.dist, data = resid_vecs))
  comm_residualsE <- resid(lm(c.dist ~ s.dist, data = resid_vecs))
  geo_residuals <- resid(lm(e.dist ~ s.dist, data = resid_vecs))
  
  ind <- substring(k, 2, 4)
  
  elep <- functional.partial.mantels.comp[functional.partial.mantels.comp$index == ind & 
                                           functional.partial.mantels.comp$gradient == "ele" & 
                                           functional.partial.mantels.comp$year == i & 
                                           functional.partial.mantels.comp$seasn == j, "pval3"]
  spap <- functional.partial.mantels.comp[functional.partial.mantels.comp$index == ind & 
                                           functional.partial.mantels.comp$gradient == "spa" & 
                                           functional.partial.mantels.comp$year == i & 
                                           functional.partial.mantels.comp$seasn == j, "pval3"]
  
  year <- i
  seasn <- ifelse(j == "J", "Wet Season", "Dry Season")
  
  #par(mfrow = c(1, 1))
  
  #ylabel <- paste("Community Dissimilarity Residuals:", k, sep = " ")
  
  ylabel <- ifelse(k == "bsor", "Functional Beta Diversity", 
         ifelse(k == "bsim", "Functional Turnover", "Functional Nestedness"))
  pcol <- ifelse(j == "J", "red", "skyblue")
  
  plot(geo_residuals, comm_residualsE, ylab = ylabel, xlab = "Elevational Distance (Residuals)", 
       pch = 21, bg = pcol, col = "white")
  eline <- ifelse(elep < 0.05, 1, 0)
  abline(lm(comm_residualsE~-1+geo_residuals), col = "black", lty = eline)
  #plot(eco_residuals, comm_residualsG, ylab = ylabel, xlab = "Geographical Distance Residuals", 
  #     pch = 21, bg = "black", col = "white")
  #sline <- ifelse(spap < 0.05, 1, 0)
  #abline(lm(comm_residualsG~-1+eco_residuals), col = "red", lty = sline)
}


partial_plot2 <- function(i, j, k){
  c.dist <- comm.dissim[[paste(i, j, k, sep = "")]]
  e.dist <- ele.dissim[[paste(i, j, sep = "")]]
  s.dist <- spa.dissim[[paste(i, j, sep = "")]]
  
  residual_vectors <- cbind(c.dist, e.dist, s.dist)
  resid_vecs <- as.data.frame(residual_vectors)
  
  comm_residualsG <- resid(lm(c.dist ~ e.dist, data = resid_vecs))
  eco_residuals <- resid(lm(s.dist ~ e.dist, data = resid_vecs))
  comm_residualsE <- resid(lm(c.dist ~ s.dist, data = resid_vecs))
  geo_residuals <- resid(lm(e.dist ~ s.dist, data = resid_vecs))
  
  ind <- substring(k, 2, 4)
  
  elep <- taxonomic.partial.mantels.comp[taxonomic.partial.mantels.comp$index == ind & 
                                           taxonomic.partial.mantels.comp$gradient == "ele" & 
                                           taxonomic.partial.mantels.comp$year == i & 
                                           taxonomic.partial.mantels.comp$seasn == j, "pval3"]
  spap <- taxonomic.partial.mantels.comp[taxonomic.partial.mantels.comp$index == ind & 
                                           taxonomic.partial.mantels.comp$gradient == "spa" & 
                                           taxonomic.partial.mantels.comp$year == i & 
                                           taxonomic.partial.mantels.comp$seasn == j, "pval3"]
  
  year <- i
  seasn <- ifelse(j == "J", "Wet Season", "Dry Season")
  
 # par(mfrow = c(1, 1))
  
  #ylabel <- paste("Community Dissimilarity Residuals:", k, sep = " ")
  
  ylabel <- ifelse(k == "bsor", "Taxonomic Beta Diversity", 
                   ifelse(k == "bsim", "Taxonomic Turnover", "Taxonomic Nestedness"))
  pcol <- ifelse(j == "J", "red", "skyblue")
  plot(geo_residuals, comm_residualsE, ylab = ylabel, xlab = "Elevational Distance (Residuals)", 
       pch = 21, bg = pcol, col = "white")
  eline <- ifelse(elep < 0.05, 1, 0)
  abline(lm(comm_residualsE~-1+geo_residuals), col = "black", lty = eline)
  #plot(eco_residuals, comm_residualsG, ylab = ylabel, xlab = "Geographical Distance Residuals", 
  #     pch = 21, bg = "black", col = "white")
  #sline <- ifelse(spap < 0.05, 1, 0)
  #abline(lm(comm_residualsG~-1+eco_residuals), col = "red", lty = sline)
}







par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))


partial_plot2(2012, "J", "bsim")
partial_plot2(2012, "J", "bnes")
partial_plot2(2012, "S", "bsim")
partial_plot2(2012, "S", "bnes")



partial_plot2F(2010, "J", "bsim")
partial_plot2F(2010, "J", "bsne")
partial_plot2F(2010, "S", "bsim")
partial_plot2F(2010, "S", "bsne")






