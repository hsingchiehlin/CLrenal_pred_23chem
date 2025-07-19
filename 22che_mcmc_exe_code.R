library(rstan)
library(tidyverse)
library(ggplot2)
library(scales)
library(patchwork)
library(plyr)
library(ggridges)
source("MCSim/function.R")
set_PATH()
makemod()


# -------------- MCMC fitting --------------
model <- "PFAS_kidneychip.model.R"
makemcsim(model)
# Prepare .in file
Transwell_conc_data <- read.csv("Data/Transwell_22che_conc_data.csv")
Transwell_conc_data$Compound[which(Transwell_conc_data$Compound == "4:2 FTS")] <- "4_2_FTS"
Transwell_conc_data$Compound[which(Transwell_conc_data$Compound == "6:2 FTS")] <- "6_2_FTS"
Transwell_conc_data$Compound[which(Transwell_conc_data$Compound == "8:2 FTS")] <- "8_2_FTS"
Transwell_conc_data$Compound[which(Transwell_conc_data$Compound == "Tenofovir TFV")] <- "Tenofovir_TFV"

Transwell_conc_data_used <- aggregate(.~ Compound + Addition.Type + Time + Compartment, 
                                      data = Transwell_conc_data[, c("Compound","Addition.Type","Time","Compartment", 
                                                                     "Input.Amount.umol", "Estimated.Concentration.uM")], mean)
#Transwell_conc_data_used$Compound[which(Transwell_conc_data_used$Compound == "4:2 FTS")] <- "4_2_FTS"
#Transwell_conc_data_used$Compound[which(Transwell_conc_data_used$Compound == "6:2 FTS")] <- "6_2_FTS"
#Transwell_conc_data_used$Compound[which(Transwell_conc_data_used$Compound == "8:2 FTS")] <- "8_2_FTS"
che.list <- unique(Transwell_conc_data_used$Compound) # i
seed.list <- c(1014004711, 1290171961, 657100342, 2068737417) # j

#  prepare mcmc.in file
for(i in 1:length(che.list)){
  
  ind.che.data <- Transwell_conc_data_used[which(Transwell_conc_data_used$Compound == che.list[i]),]
  for(j in 1:4){
    
    file.name <- paste("MCSim/",che.list[i], "_mc_", j, ".in", sep="")
    seed <- seed.list[j]
    Time_cell <- 24
    Time_conc <- "4, 24"
    Init_amount_atob <- mean(ind.che.data$Input.Amount.umol[which(ind.che.data$Addition.Type == "AtoB")])
    Conc_top_atob <- ind.che.data$Estimated.Concentration.uM[which(ind.che.data$Addition.Type == "AtoB" & ind.che.data$Compartment == "Top")]
    Conc_bottom_atob <- ind.che.data$Estimated.Concentration.uM[which(ind.che.data$Addition.Type == "AtoB" & ind.che.data$Compartment == "Bottom")]
    Conc_cell_atob <- ind.che.data$Estimated.Concentration.uM[which(ind.che.data$Addition.Type == "AtoB" & ind.che.data$Compartment == "Lysate")]
    
    Init_amount_btoa <- mean(ind.che.data$Input.Amount.umol[which(ind.che.data$Addition.Type == "BtoA")])
    Conc_top_btoa <- ind.che.data$Estimated.Concentration.uM[which(ind.che.data$Addition.Type == "BtoA" & ind.che.data$Compartment == "Top")]
    Conc_bottom_btoa <- ind.che.data$Estimated.Concentration.uM[which(ind.che.data$Addition.Type == "BtoA" & ind.che.data$Compartment == "Bottom")]
    Conc_cell_btoa <- ind.che.data$Estimated.Concentration.uM[which(ind.che.data$Addition.Type == "BtoA" & ind.che.data$Compartment == "Lysate")]
    
    cat("Integrate (Lsodes, 1e-6, 1e-6, 1);\n\n", file = file.name)
    sink(file.name, append=TRUE)
    cat('MCMC ("MCMC.default.out",      # output file
      "",                      # name of restart file
      "",                      # name of data file
      300000, 0,               # iterations, print predictions flag,
      30, 300000,             # printing frequency, iters to print\n')
    seedline <- paste("      ", seed, ");# random seed\n", sep="")
    cat(seedline)
    cat('Level { # top
  
  Distrib(p_1, Normal, -6.907755, 11.51293); #(1e-3, 1e5)
  Distrib(p_2, Normal, -6.907755, 11.51293);
  Distrib(p_3, Normal, -6.907755, 11.51293);
  Distrib(p_4, Normal, -6.907755, 11.51293);
  Distrib(GSD_error_C1_blood, LogUniform, 1.1, 10);
  Distrib(GSD_error_C1_cell, LogUniform, 1.1, 10);
  Distrib(GSD_error_C1_lumen, LogUniform, 1.1, 10);
  Distrib(GSD_error_C2_blood, LogUniform, 1.1, 10);
  Distrib(GSD_error_C2_cell, LogUniform, 1.1, 10);
  Distrib(GSD_error_C2_lumen, LogUniform, 1.1, 10);

  
  Likelihood(Data(C1_blood), LogNormal, Prediction(C1_blood), GSD_error_C1_blood);
  Likelihood(Data(C1_cell), LogNormal, Prediction(C1_cell), GSD_error_C1_cell);
  Likelihood(Data(C1_lumen), LogNormal, Prediction(C1_lumen), GSD_error_C1_lumen);
  Likelihood(Data(C2_blood), LogNormal, Prediction(C2_blood), GSD_error_C2_blood);
  Likelihood(Data(C2_cell), LogNormal, Prediction(C2_cell), GSD_error_C2_cell);
  Likelihood(Data(C2_lumen), LogNormal, Prediction(C2_lumen), GSD_error_C2_lumen);\n\n')
    cat('  Simulation { 
  area = 0.143;                # cm^2
  V_blood = 0.000235;          #  L
  V_cell = 2.529479e-08;       # (0.00065)^3*3.14*(4/3)*50000/1000
  V_lumen = 0.000075;          #  L\n')
    
    A1_lumen_init <- paste("    A1_lumen_init = Events (A1_lumen, 1, 0, Add, ", Init_amount_atob, ");\n", sep="")
    A2_blood_init <- paste("    A2_blood_init = Events (A2_blood, 1, 0, Add, ", Init_amount_btoa, ");\n", sep="")
    Print.conc <- paste("    Print(C1_blood, C1_lumen, C2_blood, C2_lumen, ", Time_conc, ");\n", sep="")
    C1_blood <- paste("    Data(C1_blood, ", paste(Conc_bottom_atob, collapse = ","), ");\n", sep="")
    C1_lumen <- paste("    Data(C1_lumen, ", paste(Conc_top_atob, collapse = ","), ");\n", sep="")
    C2_blood <- paste("    Data(C2_blood, ", paste(Conc_bottom_btoa, collapse = ","), ");\n", sep="")
    C2_lumen <- paste("    Data(C2_lumen, ", paste(Conc_top_btoa, collapse = ","), ");\n", sep="")
    Print.cell <- paste("    Print(C1_cell, C2_cell, ", Time_cell, ");\n", sep="")
    C1_cell <- paste("    Data(C1_cell, ", Conc_cell_atob, ");\n", sep="")
    C2_cell <- paste("    Data(C2_cell, ", Conc_cell_btoa, ");\n", sep="")
    
    cat(A1_lumen_init, A2_blood_init, Print.conc, C1_blood, C1_lumen, C2_blood, C2_lumen, Print.cell, C1_cell, C2_cell)
    
    cat('  }\n')
    cat('}\n')
    cat('END.')
    sink()
    
  }
}

#  prepare setpts.in file
for(i in 1:length(che.list)){
  ind.che.data <- Transwell_conc_data_used[which(Transwell_conc_data_used$Compound == che.list[i]),]
  file.name <- paste("MCSim/",che.list[i], "_setpts",".in", sep="")
  Init_amount_atob <- mean(ind.che.data$Input.Amount.umol[which(ind.che.data$Addition.Type == "AtoB")])
  Init_amount_btoa <- mean(ind.che.data$Input.Amount.umol[which(ind.che.data$Addition.Type == "BtoA")])
  
  cat("Integrate (Lsodes, 1e-6, 1e-6, 1);\n\n", file = file.name)
  sink(file.name, append=TRUE)
  out.file <- paste('"', che.list[i], '_setpts.out"', sep="")
  SetPoints <- paste('SetPoints ("",', out.file, ",0,p_1, p_2, p_3, p_4, GSD_error_C1_blood, GSD_error_C1_cell, GSD_error_C1_lumen, GSD_error_C2_blood, GSD_error_C2_cell, GSD_error_C2_lumen);\n\n", sep="")
  cat(SetPoints)
  cat('Simulation { 
  area = 0.143;                # cm^2
  V_blood = 0.000235;          #  L
  V_cell = 2.529479e-08;       # (0.00065)^3*3.14*(4/3)*50000/1000
  V_lumen = 0.000075;          #  L\n')
  
  A1_lumen_init <- paste("A1_lumen_init = Events (A1_lumen, 1, 0, Add, ", Init_amount_atob, ");\n", sep="")
  A2_blood_init <- paste("A2_blood_init = Events (A2_blood, 1, 0, Add, ", Init_amount_btoa, ");\n", sep="")
  cat(A1_lumen_init, A2_blood_init)
  
  cat('PrintStep(  
    # 1: Lumen (A) to Blood (B) study
    C1_blood,      # Quantity in blood compartment (micromoles)
    C1_cell,       # Quantity in cell lysate compartment (micromoles)
    C1_lumen,      # Quantity in lumen compartment (micromoles)
    
    # 2: Blood (B) to Lumen (A) study
    C2_blood,      # Quantity in blood compartment (micromoles)
    C2_cell,       # Quantity in cell lysate compartment (micromoles)
    C2_lumen,      # Quantity in lumen compartment (micromoles)
    
    # 3: Both add study
    #A3_blood,      # Quantity in blood compartment (micromoles)
    #A3_cell,       # Quantity in cell lysate compartment (micromoles)
    #A3_lumen,       # Quantity in lumen compartment (micromoles)
    0, 24, 1);

}

End.')
  sink()
  
}

#  MCMC modeling and generate the trace and density plots, as well as scatter plot to compare obs. and pred., for each chemical
for (i in 1:length(che.list)){
  
  file.name.1 <- paste(che.list[i], "_mc_1", ".in", sep="")
  file.name.2 <- paste(che.list[i], "_mc_2", ".in", sep="")
  file.name.3 <- paste(che.list[i], "_mc_3", ".in", sep="")
  file.name.4 <- paste(che.list[i], "_mc_4", ".in", sep="")
  out_1 <- mcsim(model, file.name.1)
  out_2 <- mcsim(model, file.name.2)
  out_3 <- mcsim(model, file.name.3)
  out_4 <- mcsim(model, file.name.4)
  
  out_all <- rbind(out_1, out_2, out_3, out_4)
  write.csv(out_all, paste("MCMC_outputs/", che.list[i], "_alliters.csv", sep = ""))
  
  R_hat.results <- monitor(mcmc_array(list(out_1, out_2, out_3, out_4)), digits = 4, warmup = 2000)
  write.csv(R_hat.results, paste("MCMC_outputs/", che.list[i], "_R_hat.results.csv", sep = ""))
  
  str <- ceiling(nrow(out_1)/5) + 1
  end <- nrow(out_1)
  j <- c(str:end) # discard burn-in
  
  X <- mcmc_array(list(out_1, out_2, out_3, out_4))[j,,] %>% matrix(nrow = 32000) 
  write.table(X, file = paste(che.list[i], "_setpts.out", sep = ""), row.names = F, sep = "\t")
  colnames(X) <- colnames(out_1)
  write.csv(X, file = paste("MCMC_outputs/", che.list[i], "_setpts.csv", sep = ""), row.names = F)
  
  setpts.file <- paste(che.list[i], "_setpts",".in", sep="")
  X_setpts <- mcsim(model, setpts.file)
  vars <- names(X_setpts)
  index <- which(vars == "C1_blood_1.1" | vars == "C2_lumen_1.25")
  simulation_result <- as.data.frame(apply(X_setpts[index[1]:index[2]], 2, quantile,  c(0.5, 0.05, 0.95)) %>% t())
  colnames(simulation_result) <- c("median", "LCL", "UCL")
  simulation_result$Time <- rep(c(0:24), 6)
  simulation_result$Compartment <- rep(rep(c("Blood", "Cell", "Lumen"), each = 25), 2)
  simulation_result$Addition.Type <- rep(c("AtoB", "BtoA"), each = 75)
  write.csv(simulation_result, file = paste("MCMC_outputs/", che.list[i], "_simulation_result.csv", sep = ""), row.names = F)
  
  
  ind.che.data <- Transwell_conc_data[which(Transwell_conc_data$Compound == che.list[i]),]
  ind.che.data$Compartment <- ind.che.data$Compartment_new
  plot.1 <- 
    ggplot(simulation_result)+
    geom_ribbon(aes(x = Time, ymin = LCL, ymax = UCL), fill="#E3ECFE") +
    geom_line(aes(x = Time, y = median), linewidth = 1.2, color = "#1E69F4")+
    geom_point(data = ind.che.data, aes(x = Time, y = Estimated.Concentration.uM), size = 2, color = "#0B275B", alpha = 0.8)+
    facet_grid(Compartment~Addition.Type, scales = "free")+
    xlab("Time") +
    ylab("Concentration (uM)")+
    ggtitle(che.list[i])+
    scale_x_continuous(lim = c(0, 24),
                       breaks = c(0, 4, 8, 12, 16, 20, 24),
                       labels = c(0, 4, 8, 12, 16, 20, 24)) +
    theme_bw() +
    theme(text=element_text(family="sans", face="plain", color="#000000", size=12, hjust=0.5, vjust=0.5), 
          plot.title = element_text(face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_line(color = "black"),
          axis.title = element_text(color = "black"),
          axis.text = element_text(color = "black"),
          strip.text = element_text(color = "black"))
  
  ind.che.data$seach <- paste(ind.che.data$Addition.Type, ind.che.data$Time, ind.che.data$Compartment, sep = "_")
  simulation_result$seach <- paste(simulation_result$Addition.Type, simulation_result$Time, simulation_result$Compartment, sep = "_")                              
  for(k in 1:length(ind.che.data$Estimated.Concentration.uM)){
    ind.che.data$Pred_conc[k] <- simulation_result$median[which(simulation_result$seach == ind.che.data$seach[k])]
  }
  ind.che.data$ratio <- log10(ind.che.data$Pred_conc/ind.che.data$Estimated.Concentration.uM)
  ind.che.data$abs.dev <- abs(ind.che.data$ratio)
  MAD <- mean(ind.che.data$abs.dev)
  spea.corr <- cor(log10(ind.che.data$Pred_conc), log10(ind.che.data$Estimated.Concentration.uM), method = "spearman")
  max <- 10^round(log10(max(c(ind.che.data$Pred_conc, ind.che.data$Estimated.Concentration.uM)))+0.5)
  min <- 10^round(log10(min(c(ind.che.data$Pred_conc, ind.che.data$Estimated.Concentration.uM)))-0.5)
  plot.2 <- 
    ggplot(ind.che.data, aes(x = Estimated.Concentration.uM, y = Pred_conc))+
    geom_point(aes(color = Compartment, shape = Addition.Type), alpha=0.8, size = 2) +
    xlab("Experimental concentration") +
    ylab("Predicted concentration")+
    annotate(geom="text", x=10^(log10(min)+1), y=max, label=paste("MAD = ", round(MAD,2), ",", " Rho = ", round(spea.corr, 2), sep = ""), color="black")+
    geom_abline(intercept = 0, slope = 1, linetype=2)+
    scale_x_log10(lim = c(min, max),
                  breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(lim = c(min, max),
                  breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                  labels = trans_format("log10", math_format(10^.x))) +
    theme_bw() +
    theme(text=element_text(family="sans", face="plain", color="#000000", size=12, hjust=0.5, vjust=0.5), 
          legend.title = element_text(color = "black", size=10),
          legend.text = element_text(color = "black", size=8),
          legend.background = element_blank(),
          legend.key.size = unit(0.4, 'cm'),
          legend.position = c(0.8, 0.3), 
          plot.title = element_text(face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_line(color = "black"),
          axis.title = element_text(color = "black"),
          axis.text = element_text(color = "black"),
          strip.text = element_text(color = "black", size=12))
  
  alliters <- stack(out_all[,2:5])
  alliters$iter <- rep(rep(c(0:10000), 4), 4)
  alliters$chain.num <- rep(rep(c(1:4), each = 10001), 4)
  alliters$chain.num <- factor(alliters$chain.num)
  plot.3 <- 
    ggplot(data = alliters)+
    geom_line(aes(x = iter, y = values, color = chain.num), size = 0.6) + 
    geom_vline(xintercept = 2000, color = "red", linetype = "dashed")+
    facet_wrap(ind ~ ., nrow = 2, scales = "free_y")+
    xlab("Iteration")+
    ylab("Value")+
    scale_color_manual(values = c("#768FDf", "#ec7b7d","#FDC168","#91d1be"))+
    theme_bw() +
    theme(text=element_text(family="sans", face="plain", color="#000000", size=12, hjust=0.5, vjust=0.5), 
          legend.position = "none",
          plot.title = element_text(face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_line(color = "black"),
          axis.title = element_text(color = "black"),
          axis.text = element_text(color = "black"),
          strip.text = element_text(color = "black"))
  
  postdis <- stack(X_setpts[,2:5])
  postdis$chain.num <- rep(rep(c(1:4), each = 8000),4)
  postdis$chain.num <- factor(postdis$chain.num)
  postdis$values <- exp(postdis$values)
  R_hat <- data.frame(ind=c("p_1", "p_2", "p_3", "p_4"),
                      R_hat = R_hat.results$Rhat[2:5])
  max.postdis <- 10^round(log10(max(postdis$values))+0.5)
  min.postdis <- 10^round(log10(min(postdis$values))-0.5)
  plot.4 <- 
    ggplot(data = postdis)+
    geom_density_ridges(aes(x = values, y = ind, color = chain.num, fill = chain.num), size = 0.8, alpha = 0.5, height = 1, scale = 1.25) + 
    geom_text(
      data    = R_hat,
      mapping = aes(x = 10^(log10(max.postdis)-1), y = ind, label = round(R_hat, 3)),
      color = "black", vjust = -0.5)+
    scale_x_log10(lim = c(min.postdis, max.postdis),
                  breaks = trans_breaks("log10", function(x) 10^x, n = 8),
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("Value (cm/hr)")+
    ylab("Parameter")+
    scale_color_manual(values = c("#EAFDFC", "#BFEAF5","#91D8E4","#82AAE3"))+
    scale_fill_manual(values = c("#EAFDFC", "#BFEAF5","#91D8E4","#82AAE3"))+
    theme_bw() +
    theme(text=element_text(family="sans", face="plain", color="#000000", size=12, hjust=0.5, vjust=0.5), 
          legend.position = "none",
          plot.title = element_text(face="bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = "gray", linetype = "dashed"),
          panel.grid.minor = element_blank(),
          axis.ticks = element_line(color = "black"),
          axis.title = element_text(color = "black"),
          axis.text = element_text(color = "black"),
          strip.text = element_text(color = "black"))
  
  plot.1 + plot.3 + plot.2 + plot.4 + plot_layout(widths = c(1/2, 1/2), heights = c(1/2, 1/2))
  ggsave(filename = paste("MCMC_plots/", che.list[i], "_mcmc_all_plot.pdf", sep = ""), width = 9.5, height = 8.2)
  
  print(che.list[i])
}

#  Overall scatter plot to compare obs. and pred.
compare_obs_sim_dat <- c()
for (i in 1:length(che.list)){
  ind.che.data <- Transwell_conc_data[which(Transwell_conc_data$Compound == che.list[i]),]
  ind.che.data$Compartment <- ind.che.data$Compartment_new
  ind.che.data$seach <- paste(ind.che.data$Addition.Type, ind.che.data$Time, ind.che.data$Compartment, sep = "_")
  simulation_result <- read.csv(paste("MCMC_outputs/", che.list[i], "_simulation_result.csv", sep = ""))
  simulation_result$seach <- paste(simulation_result$Addition.Type, simulation_result$Time, simulation_result$Compartment, sep = "_")                              
  for(k in 1:length(ind.che.data$Estimated.Concentration.uM)){
    ind.che.data$Pred_conc[k] <- simulation_result$median[which(simulation_result$seach == ind.che.data$seach[k])]
  }
  ind.che.data$ratio <- log10(ind.che.data$Pred_conc/ind.che.data$Estimated.Concentration.uM)
  ind.che.data$abs.dev <- abs(ind.che.data$ratio)
  compare_obs_sim_dat <- rbind(compare_obs_sim_dat, ind.che.data)
}
write.csv(compare_obs_sim_dat, file = "MCMC_outputs/compare_obs_sim_dat.csv")

Stat_dat <- ddply(compare_obs_sim_dat, .(Compound), summarize, 
                  MAD = mean(abs.dev),
                  ME = mean(ratio),
                  Pear.Corr = cor(log10(Estimated.Concentration.uM), log10(Pred_conc)),
                  Spea.Corr = cor(log10(Estimated.Concentration.uM), log10(Pred_conc), method = "spearman"))
max <- 10^round(log10(max(c(compare_obs_sim_dat$Pred_conc, compare_obs_sim_dat$Estimated.Concentration.uM)))+0.5)
min <- 10^round(log10(min(c(compare_obs_sim_dat$Pred_conc, compare_obs_sim_dat$Estimated.Concentration.uM)))-0.5)
ggplot(compare_obs_sim_dat, aes(x = Estimated.Concentration.uM, y = Pred_conc))+
  geom_point(aes(color = Compartment, shape = Addition.Type), alpha=0.8, size = 2) +
  geom_text(
    data    = Stat_dat,
    mapping = aes(x=10^(log10(min)+2), y=10^(log10(max)-0.5),
                  label = paste("r = ", round(Pear.Corr, 2), sep = "")),
    color = "black",
    size = 3.5
  )+
  geom_text(
    data    = Stat_dat,
    mapping = aes(x=10^(log10(min)+2), y=10^(log10(max)-1.5),
                  label = paste("rho = ", round(Spea.Corr, 2), sep = "")),
    color = "black",
    size = 3.5
  )+
  geom_text(
    data    = Stat_dat,
    mapping = aes(x=10^(log10(max)-2), y=10^(log10(min)+1.5),
                  label = paste("MAD = ", round(MAD, 2), sep = "")),
    color = "black",
    size = 3.5
  )+
  geom_text(
    data    = Stat_dat,
    mapping = aes(x=10^(log10(max)-2), y=10^(log10(min)+0.5),
                  label = paste("ME = ", round(ME, 2), sep = "")),
    color = "black",
    size = 3.5
  )+
  xlab("Experimental concentration") +
  ylab("Predicted concentration")+
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  geom_abline(intercept = 1, slope = 1, linetype=2, color="lightgray")+
  geom_abline(intercept = -1, slope = 1, linetype=2, color="lightgray")+
  scale_x_log10(lim = c(min, max),
                breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(lim = c(min, max),
                breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", math_format(10^.x))) +
  facet_wrap(~Compound, nrow = 3)+
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=12, hjust=0.5, vjust=0.5), 
        legend.title = element_text(color = "black", size=10),
        legend.text = element_text(color = "black", size=8),
        legend.background = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        #legend.position = c(0.8, 0.3), 
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=10))
ggsave("MCMC_plots/Comp_allobs_pred.scatter.jpeg", width = 12.8, height = 5.5)
ggsave("MCMC_plots/Comp_allobs_pred.scatter.pdf", width = 12.8, height = 5.5)


compare_obs_sim_dat_mean <- ddply(compare_obs_sim_dat, .(Compound,Addition.Type,Compartment,Time), summarize, 
                                  Pred_conc = mean(Pred_conc),
                                  Estimated.Concentration.uM = mean(Estimated.Concentration.uM))
compare_obs_sim_dat_mean$ratio <- log10(compare_obs_sim_dat_mean$Pred_conc/compare_obs_sim_dat_mean$Estimated.Concentration.uM)
compare_obs_sim_dat_mean$abs.dev <- abs(compare_obs_sim_dat_mean$ratio)

Stat_dat <- ddply(compare_obs_sim_dat_mean, .(Compound), summarize, 
                  MAD = mean(abs.dev),
                  ME = mean(ratio),
                  Pear.Corr = cor(log10(Estimated.Concentration.uM), log10(Pred_conc)),
                  Spea.Corr = cor(log10(Estimated.Concentration.uM), log10(Pred_conc), method = "spearman"))
ggplot(compare_obs_sim_dat_mean, aes(x = Estimated.Concentration.uM, y = Pred_conc))+
  geom_point(aes(color = Compartment, shape = Addition.Type), alpha=0.8, size = 2) +
  geom_text(
    data    = Stat_dat,
    mapping = aes(x=10^(log10(min)+2), y=10^(log10(max)-0.5),
                  label = paste("r = ", round(Pear.Corr, 2), sep = "")),
    color = "black",
    size = 3.5
  )+
  geom_text(
    data    = Stat_dat,
    mapping = aes(x=10^(log10(min)+2), y=10^(log10(max)-1.5),
                  label = paste("rho = ", round(Spea.Corr, 2), sep = "")),
    color = "black",
    size = 3.5
  )+
  geom_text(
    data    = Stat_dat,
    mapping = aes(x=10^(log10(max)-2), y=10^(log10(min)+1.5),
                  label = paste("MAD = ", round(MAD, 2), sep = "")),
    color = "black",
    size = 3.5
  )+
  geom_text(
    data    = Stat_dat,
    mapping = aes(x=10^(log10(max)-2), y=10^(log10(min)+0.5),
                  label = paste("ME = ", round(ME, 2), sep = "")),
    color = "black",
    size = 3.5
  )+
  xlab("Experimental concentration") +
  ylab("Predicted concentration")+
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  geom_abline(intercept = 1, slope = 1, linetype=2, color="lightgray")+
  geom_abline(intercept = -1, slope = 1, linetype=2, color="lightgray")+
  scale_x_log10(lim = c(min, max),
                breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(lim = c(min, max),
                breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", math_format(10^.x))) +
  facet_wrap(~Compound, nrow = 3)+
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=12, hjust=0.5, vjust=0.5), 
        legend.title = element_text(color = "black", size=10),
        legend.text = element_text(color = "black", size=8),
        legend.background = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        #legend.position = c(0.8, 0.3), 
        plot.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black", size=10))
ggsave("MCMC_plots/Comp_meanobs_pred.scatter.jpeg", width = 12.8, height = 5.5)
ggsave("MCMC_plots/Comp_meanobs_pred.scatter.pdf", width = 12.8, height = 5.5)


                        

