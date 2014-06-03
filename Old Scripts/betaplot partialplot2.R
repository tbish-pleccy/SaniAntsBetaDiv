comm.dissim
fun.dissim
ele.dissim
spa.dissim



beta.plot <- function(year, season, metric, composition, gradient){
  
  if (composition == "taxonomic"){
    comp.dist.list <- comm.dissim
  }
  else {
    comp.dist.list <- fun.dissim
  }
  
  if (gradient == "elevation"){
    grad.dist.list <- ele.dissim
  }
  else {
    grad.dist.list <- spa.dissim
  }
  
  comm.dist <- comp.dist.list[[paste(year, season, metric, sep = "")]]
  grad.dist <- grad.dist.list[[paste(year, season, sep = "")]]
  
  ylabel <- metric
  xlabel <- gradient
  mlabel <- paste(year, season)
  plot(comm.dist ~ grad.dist, pch = 21, col = "white", bg = "black", 
       ylab = ylabel, xlab = xlabel, main = mlabel, ylim = c(0, 1))
  if (mantel(comm.dist ~ grad.dist)[4] < 0.05) {
    abline((lm(comm.dist ~ grad.dist)))}
  
}



par(mfrow = c(2, 4))
beta.plot(2006, "J", "bsor", "taxonomic", "elevation")
beta.plot(2007, "J", "bsor", "taxonomic", "elevation")
beta.plot(2008, "J", "bsor", "taxonomic", "elevation")
beta.plot(2009, "J", "bsor", "taxonomic", "elevation")
beta.plot(2010, "J", "bsor", "taxonomic", "elevation")
beta.plot(2011, "J", "bsor", "taxonomic", "elevation")
beta.plot(2012, "J", "bsor", "taxonomic", "elevation")

par(mfrow = c(2, 4))
beta.plot(2006, "S", "bsor", "taxonomic", "elevation")
beta.plot(2007, "S", "bsor", "taxonomic", "elevation")
beta.plot(2008, "S", "bsor", "taxonomic", "elevation")
beta.plot(2009, "S", "bsor", "taxonomic", "elevation")
beta.plot(2010, "S", "bsor", "taxonomic", "elevation")
beta.plot(2011, "S", "bsor", "taxonomic", "elevation")
beta.plot(2012, "S", "bsor", "taxonomic", "elevation")

par(mfrow = c(2, 4))
beta.plot(2006, "J", "bsim", "taxonomic", "elevation")
beta.plot(2007, "J", "bsim", "taxonomic", "elevation")
beta.plot(2008, "J", "bsim", "taxonomic", "elevation")
beta.plot(2009, "J", "bsim", "taxonomic", "elevation")
beta.plot(2010, "J", "bsim", "taxonomic", "elevation")
beta.plot(2011, "J", "bsim", "taxonomic", "elevation")
beta.plot(2012, "J", "bsim", "taxonomic", "elevation")

par(mfrow = c(2, 4))
beta.plot(2006, "S", "bsim", "taxonomic", "elevation")
beta.plot(2007, "S", "bsim", "taxonomic", "elevation")
beta.plot(2008, "S", "bsim", "taxonomic", "elevation")
beta.plot(2009, "S", "bsim", "taxonomic", "elevation")
beta.plot(2010, "S", "bsim", "taxonomic", "elevation")
beta.plot(2011, "S", "bsim", "taxonomic", "elevation")
beta.plot(2012, "S", "bsim", "taxonomic", "elevation")

par(mfrow = c(2, 4))
beta.plot(2006, "J", "bnes", "taxonomic", "elevation")
beta.plot(2007, "J", "bnes", "taxonomic", "elevation")
beta.plot(2008, "J", "bnes", "taxonomic", "elevation")
beta.plot(2009, "J", "bnes", "taxonomic", "elevation")
beta.plot(2010, "J", "bnes", "taxonomic", "elevation")
beta.plot(2011, "J", "bnes", "taxonomic", "elevation")
beta.plot(2012, "J", "bnes", "taxonomic", "elevation")

par(mfrow = c(2, 4))
beta.plot(2006, "S", "bnes", "taxonomic", "elevation")
beta.plot(2007, "S", "bnes", "taxonomic", "elevation")
beta.plot(2008, "S", "bnes", "taxonomic", "elevation")
beta.plot(2009, "S", "bnes", "taxonomic", "elevation")
beta.plot(2010, "S", "bnes", "taxonomic", "elevation")
beta.plot(2011, "S", "bnes", "taxonomic", "elevation")
beta.plot(2012, "S", "bnes", "taxonomic", "elevation")




















partial_plot2(2012, "J", "bsim")

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
