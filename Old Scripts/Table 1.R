###
# Create Table 1
###

mean.se <- function(x){
  mean.x <- round(mean(x), 2)
  se.x <- round(sd(x)/sqrt(length(x)), 2)
  result <- paste(mean.x, "±", se.x , sep = " ")
}
mean.se(x)

table1 <- NULL
for (comp in c("Taxonomic", "Functional")){
  for (seasn in unique(site.info$seasn)){
    for (metric in c("bsor", "bsim", "bnes")){
      
      if (comp == "Functional" & metric == "bnes"){
        metric <- "bsne"
      }
      
      if (comp == "Taxonomic"){
        av.pw.frame <- tax.mantels[tax.mantels$metric == metric & 
                                     tax.mantels$seasn == seasn, ] 
        r2.frame <- r2.summary[r2.summary$metric == metric & 
                                 r2.summary$seasn == seasn, ]
        hp.frame <- hp.summary[hp.summary$metric == metric & 
                                 hp.summary$seasn == seasn, ]
      }
      else {
        av.pw.frame <- fun.mantels[fun.mantels$metric == metric & 
                                     fun.mantels$seasn == seasn, ]
        r2.frame <- fr2.summary[fr2.summary$metric == metric & 
                                fr2.summary$seasn == seasn, ]
        fhp.frame <- fhp.summary[fhp.summary$metric == metric & 
                                 fhp.summary$seasn == seasn, ]
      }
      
      av.pw <- mean.se(av.pw.frame$av.pairwise)
      av.R2 <- mean.se(r2.frame$r2)
      e <- mean.se(hp.frame[hp.frame$var == "e", "I"])
      e2 <- mean.se(hp.frame[hp.frame$var == "e2", "I"])
      s <- mean.se(hp.frame[hp.frame$var == "s", "I"])
      s2 <- mean.se(hp.frame[hp.frame$var == "s2", "I"])
      
      table1.temp <- data.frame(comp, metric, seasn, av.pw, av.R2, e, e2, s, s2)
      table1 <- rbind(table1, table1.temp)
      
    }  
  }
}  

tapply(tax.mantels$av.pairwise, list(tax.mantels$seasn, tax.mantels$metric), mean)
tax.mantels

table1