library(tidyverse)
library(ggplot2)
library(scales)
library(patchwork)
library(plyr)
library(ggridges)
library(ggh4x)
library(ggpattern)
library(gridpattern)
library(xlsx)

# --------- Data prepare ----------
# 2D culture (include 22 chemical study, Lucie's non target analysis, and Lin et al data)
Che22_Kidney_data <- read.csv("Data/Kidney_CL_2D_alldata.csv")
Che22_Kidney_data <- ddply(Che22_Kidney_data, .(Study, Invitro.type, Compound, Compound.type, Cell.type), summarize, 
                           ratio = mean(ratio), CL_SD = sd(log10(CL), na.rm = T), CL = mean(CL), CL_obs = mean(CL_obs))
Che22_Kidney_data$CL_UB <- 10^(log10(Che22_Kidney_data$CL) + 1.645*(Che22_Kidney_data$CL_SD))
Che22_Kidney_data$CL_LB <- 10^(log10(Che22_Kidney_data$CL) - 1.645*(Che22_Kidney_data$CL_SD))
Che22_Kidney_data <- Che22_Kidney_data[,-6]
Che22_Kidney_data$Study[which(Che22_Kidney_data$Study == "Lin et al., 2024")] <- "Lin et al."
#Che22_Kidney_data$Invitro.type <- "2D plate"
Che22_Kidney_data$Compound[which(Che22_Kidney_data$Compound == "Tenofovir (TFV)")] <- "Tenofovir_TFV"
Che22_Kidney_data$Compound[which(Che22_Kidney_data$Compound == "Gentamycin")] <- "Gentamicin"

# Transwell for 22 chemical study
Transwell_conc_data <- read.csv("Transwell_22che_conc_data.csv")
Transwell_conc_data$Compound[which(Transwell_conc_data$Compound == "4:2 FTS")] <- "4_2_FTS"
Transwell_conc_data$Compound[which(Transwell_conc_data$Compound == "6:2 FTS")] <- "6_2_FTS"
Transwell_conc_data$Compound[which(Transwell_conc_data$Compound == "8:2 FTS")] <- "8_2_FTS"
Transwell_conc_data$Compound[which(Transwell_conc_data$Compound == "Tenofovir TFV")] <- "Tenofovir_TFV"
che.list <- unique(Transwell_conc_data$Compound) # i
fub <- read.csv("Data/fup_22che_dat.csv")
fub$Che.type <- c(rep("PFAS", 14), rep("Non-PFAS", 8))
Qk <- 29.64556*1000*24
Qu <- 6.7732*1000*24*0.3

P_ratio_CL_22chedat <- c()
for (i in 1:length(che.list)){
  
  setpts <- read.csv(paste("MCMC_outputs/", che.list[i], "_setpts.csv", sep = ""))
  colnames(setpts)[2:5] <- c("p1", "p2", "p3", "p4")
  setpts[, 2:5] <- exp(setpts[, 2:5])
  #setpts$chain.num <- rep(c(1:4), each = 8000)
  setpts$P_ratio <- (setpts$p2*setpts$p4)/(setpts$p1*setpts$p3)
  setpts$CL <- fub$fub[which(fub$Compound == che.list[i])]*(Qk*(Qu*setpts$P_ratio)/(Qu*setpts$P_ratio+(Qk-Qu)))/70
  setpts$Compound <- che.list[i]
  setpts$Compound.type <- fub$Che.type[which(fub$Compound == che.list[i])]
  setpts$Cell.type <- "OAT1"
  setpts$CL_obs <- fub$CL_obs[which(fub$Compound == che.list[i])]
  P_ratio_CL_22chedat <- rbind(P_ratio_CL_22chedat, setpts)
  
}
Transwell_Summ_22chedat <- ddply(P_ratio_CL_22chedat, .(Compound,Compound.type,Cell.type,CL_obs), summarize,
                               P_ratio.median = median(P_ratio, na.rm = T),
                               P_ratio.Q5 = quantile(P_ratio, probs = c(0.05)),
                               P_ratio.Q25 = quantile(P_ratio, probs = c(0.25)),
                               P_ratio.Q75 = quantile(P_ratio, probs = c(0.75)),
                               P_ratio.Q95 = quantile(P_ratio, probs = c(0.95)),
                               CL.median = median(CL, na.rm = T),
                               CL.Q5 = quantile(CL, probs = c(0.05)),
                               CL.Q25 = quantile(CL, probs = c(0.25)),
                               CL.Q75 = quantile(CL, probs = c(0.75)),
                               CL.Q95 = quantile(CL, probs = c(0.95)))

Transwell_Summ_22chedat_boxplot <- data.frame(Study = "22 che study",
                                              Compound = rep(Transwell_Summ_22chedat$Compound, 5),
                                              Compound.type = rep(Transwell_Summ_22chedat$Compound.type, 5),
                                              Cell.type = rep(Transwell_Summ_22chedat$Cell.type, 5),
                                              ratio = stack(Transwell_Summ_22chedat[,5:9])$values,
                                              CL = stack(Transwell_Summ_22chedat[,10:14])$values,
                                              stat = rep(c("median", "Q5", "Q25", "Q75", "Q95"), each = 21),
                                              CL_obs = rep(Transwell_Summ_22chedat$CL_obs, 5))
#Transwell_Summ_22chedat_boxplot$Compound[which(Transwell_Summ_22chedat_boxplot$Compound == "Gentamicin")] <- "Gentamycin"
#Transwell_Summ_22chedat_boxplot$Compound[which(Transwell_Summ_22chedat_boxplot$Compound == "Tenofovir TFV")] <- "Tenofovir (TFV)"


# Transwell for Lin et al. data
che.list.3pfas <- c("pfos", "pfbs", "pfhxa")
che.lable.3pfas <- c("PFOS", "PFBS", "PFHxA")
fub.3pfas <- c(0.0049, 0.0128, 0.0068) #using fraction ubound in plasma -> Smeltz et al: UC assay + UPLC-MS/MS 

P_ratio_CL_3chedat <- c()
for (i in 1:length(che.list.3pfas)){
  
  setpts <- read.csv(paste("MCMC_outputs/", che.list.3pfas[i], "_Transwell_single_setpts.csv", sep = ""))
  colnames(setpts)[2:5] <- c("p1", "p2", "p3", "p4")
  setpts[, 2:5] <- exp(setpts[, 2:5])
  #setpts$chain.num <- rep(c(1:4), each = 8000)
  setpts$P_ratio <- (setpts$p2*setpts$p4)/(setpts$p1*setpts$p3)
  setpts$CL <- fub.3pfas[i]*(Qk*(Qu*setpts$P_ratio)/(Qu*setpts$P_ratio+(Qk-Qu)))/70
  setpts$Compound <- che.lable.3pfas[i]
  setpts$Compound.type <- "PFAS"
  setpts$Cell.type <- "OAT1"
  setpts$CL_obs <- fub$CL_obs[which(fub$Compound == che.lable.3pfas[i])]
  P_ratio_CL_3chedat <- rbind(P_ratio_CL_3chedat, setpts)
  
}
Transwell_Summ_3chedat <- ddply(P_ratio_CL_3chedat, .(Compound,Compound.type,Cell.type,CL_obs), summarize,
                                 P_ratio.median = median(P_ratio, na.rm = T),
                                 P_ratio.Q5 = quantile(P_ratio, probs = c(0.05)),
                                 P_ratio.Q25 = quantile(P_ratio, probs = c(0.25)),
                                 P_ratio.Q75 = quantile(P_ratio, probs = c(0.75)),
                                 P_ratio.Q95 = quantile(P_ratio, probs = c(0.95)),
                                 CL.median = median(CL, na.rm = T),
                                 CL.Q5 = quantile(CL, probs = c(0.05)),
                                 CL.Q25 = quantile(CL, probs = c(0.25)),
                                 CL.Q75 = quantile(CL, probs = c(0.75)),
                                 CL.Q95 = quantile(CL, probs = c(0.95)))

Transwell_Summ_3chedat_boxplot <- data.frame(Study = "Lin et al.",
                                              Compound = rep(Transwell_Summ_3chedat$Compound, 5),
                                              Compound.type = rep(Transwell_Summ_3chedat$Compound.type, 5),
                                              Cell.type = rep(Transwell_Summ_3chedat$Cell.type, 5),
                                              ratio = stack(Transwell_Summ_3chedat[,5:9])$values,
                                              CL = stack(Transwell_Summ_3chedat[,10:14])$values,
                                              stat = rep(c("median", "Q5", "Q25", "Q75", "Q95"), each = 3),
                                              CL_obs = rep(Transwell_Summ_3chedat$CL_obs, 5))

Transwell_ALL_data <- rbind(Transwell_Summ_22chedat_boxplot, Transwell_Summ_3chedat_boxplot)
Transwell_ALL_data$Invitro.type <- "Transwell"


# --------- Plotting: P ratio ----------
quantiles <- function(x) {
  r <- sort(x)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

Che22_Kidney_2D_ratiodata <- read.csv("Data/Kidney_CL_2D_alldata.csv")
Che22_Kidney_2D_ratiodata <- ddply(Che22_Kidney_2D_ratiodata, .(Study, Invitro.type, Compound, Compound.type, Cell.type), summarize, 
                                   ratio.mean = mean(ratio), ratio_SD = sd(log10(ratio), na.rm = T))
Che22_Kidney_2D_ratiodata$ratio_UB <- 10^(log10(Che22_Kidney_2D_ratiodata$ratio.mean) + 1.645*(Che22_Kidney_2D_ratiodata$ratio_SD))
Che22_Kidney_2D_ratiodata$ratio_LB <- 10^(log10(Che22_Kidney_2D_ratiodata$ratio.mean) - 1.645*(Che22_Kidney_2D_ratiodata$ratio_SD))


Che22_Kidney_2D_ratiodata$Compound[which(Che22_Kidney_2D_ratiodata$Compound == "Tenofovir (TFV)")] <- "Tenofovir_TFV"
Che22_Kidney_2D_ratiodata$Compound[which(Che22_Kidney_2D_ratiodata$Compound == "Gentamycin")] <- "Gentamicin"
Che22_Kidney_2D_ratiodata <- Che22_Kidney_2D_ratiodata[which(Che22_Kidney_2D_ratiodata$Study == "22 che study"),]

ggplot(Transwell_ALL_data[which(Transwell_ALL_data$Study == "22 che study"),], aes(y=Compound, x=log10(ratio),fill = Cell.type)) + #reorder(Compound, values)
  stat_summary(fun.data = quantiles, geom="errorbar", position = position_dodge(0.8), size=0.35, width=0.3)+
  stat_summary(fun.data = quantiles, geom="boxplot", position = position_dodge(0.8), size=0.35, width=0.7)+
  geom_point(data = Che22_Kidney_2D_ratiodata, aes(x = log10(ratio.mean), y = Compound, color = Cell.type),size = 2,position = position_dodge(width=0.75)) +
  geom_errorbarh(data = Che22_Kidney_2D_ratiodata, aes(y = Compound, xmin = log10(ratio_LB), xmax = log10(ratio_UB), color = Cell.type), 
                 inherit.aes = F, position = position_dodge(width=0.75)) +
  geom_vline(xintercept = 0, linetype = 2)+
  xlab("log10(P ratio)") +
  ylab("") +
  facet_nested(Compound.type~Invitro.type,nest_line = element_line(linetype = 1), scales = "free_y", space = "free", render_empty = F)+
  #scale_shape_manual(values = c(16, 21))+
  #ggtitle("A. TERT1 RPTEC - Parent")+
  #scale_shape_manual(values = c(21, 24))+
  #scale_color_manual(values = c("#167CF2","#6FCC33","#F28E16"))+
  #scale_fill_manual(values = c("#167CF2","#6FCC33","#F28E16"))+
  #coord_cartesian(xlim = c(1E-4, 2E+3), ylim = c(1E-4, 2E+3))+
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=12),
        title = element_text(color = "black", size=12),
        plot.title = element_text(color = "black", face = "bold", size=14),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 0.6),
        #legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size=10),
        legend.margin = margin(-0.5,0,0,0, unit="cm"),
        #legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        #legend.key.size = unit(0.4, 'cm'),
        strip.text.x = element_text(color = "black", size=12),
        strip.text.y = element_text(color = "black", size=12, angle = 0),
        strip.background = element_blank(),
        axis.ticks = element_line(color = "black", size = 0.3),
        axis.title = element_text(color = "black", size=12),
        axis.text = element_text(color = "black", size=9))
ggsave("Pratio_22che.jpeg",width = 9, height = 5.5)
write.csv(Che22_Kidney_2D_ratiodata, "Results/Che22_Kidney_2D_ratiodata_07102025.csv")

Che22_Kidney_Transwell_ratiodata <- P_ratio_CL_22chedat[, c(17,15)]
Che22_Kidney_Transwell_ratiodata$iter <- rep(c(1:32000), 21)
Che22_Kidney_Transwell_ratiodata <-  pivot_wider(data = Che22_Kidney_Transwell_ratiodata,
                                                 names_from = Compound,
                                                 values_from = P_ratio)
write.csv(Che22_Kidney_Transwell_ratiodata, "Results/Che22_Kidney_Transwell_ratiodata.csv")
write.csv(Transwell_ALL_data, "Results/Che22_Kidney_Transwell_ratiodata_quantile.csv")





# --------- Plotting ----------
Transwell_CL_data <- Transwell_ALL_data[which(Transwell_ALL_data$stat == "median"), -7]
Transwell_CL_data$CL_UB <- Transwell_ALL_data$CL[which(Transwell_ALL_data$stat == "Q95")]
Transwell_CL_data$CL_LB <- Transwell_ALL_data$CL[which(Transwell_ALL_data$stat == "Q5")]

# PLOT 1 ---- point plot for CL (only compare 22chemical study and Lin et al.)
CL_Che22_dat <- rbind(Che22_Kidney_data[,c(1:5, 7:10)], Transwell_CL_data[,c(1:4, 6:10)])
CL_Che22_dat <- CL_Che22_dat[-which(CL_Che22_dat$Study == "Lucie - single study"),]
CL_Che22_dat$Study[which(CL_Che22_dat$Study == "22 che study")] <- "This study"
CL_Che22_dat$Study <- factor(CL_Che22_dat$Study, levels = c("This study", "Lin et al.")) 
CL_Che22_dat$Cell.type <- factor(CL_Che22_dat$Cell.type, levels = c("Parent", "OAT1"))
CL_Che22_dat$label <- NA
CL_Che22_dat$order <- NA
for (i in 1:length(CL_Che22_dat$Study)) {
  CL_Che22_dat$label[i] <- fub$label[which(fub$Compound == CL_Che22_dat$Compound[i])]
  CL_Che22_dat$order[i] <- fub$Order[which(fub$Compound == CL_Che22_dat$Compound[i])]
}
ggplot(CL_Che22_dat, aes(x = CL, y = reorder(label, -order)))+
  geom_point(aes(x = CL_obs, y = reorder(label, -order)), shape = 17, size = 2.5) +
  geom_errorbar(aes(xmin = CL_LB, xmax = CL_UB, color = Cell.type), width = 0, size = 1) +
  geom_point(aes(color = Cell.type), fill = "white", size = 2.5) +
  xlab("Renal clearance (mL/Kg per day)") +
  ylab("") +
  facet_nested(Study~Invitro.type, nest_line = element_line(linetype = 1), render_empty = F, scales = "free_y", space = "free_y")+
  scale_shape_manual(values = c(21, 16))+
  #scale_x_log10() + 
  scale_x_log10(lim = c(1E-5, 2E+4),
                breaks = trans_breaks("log10", function(x) 10^x, n = 6),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=12),
        title = element_text(color = "black", size=12),
        plot.title = element_text(color = "black", face = "bold", size=14),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(size = 0.6),
        #legend.position = "none",
        #legend.title = element_blank(),
        legend.text = element_text(color = "black", size=10),
        #legend.margin = margin(-0.5,0,0,0, unit="cm"),
        #legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        #legend.key.size = unit(0.4, 'cm'),
        strip.text.x = element_text(color = "black", size=12),
        strip.text.y = element_text(color = "black", size=12, angle = 0),
        strip.background = element_blank(),
        axis.ticks = element_line(color = "black", size = 0.3),
        axis.title = element_text(color = "black", size=12),
        axis.text = element_text(color = "black", size=9))
ggsave("CL_pointplot_Comparemultiplestudy.jpeg",width = 9.5, height = 6.2)
write.xlsx(CL_Che22_dat, "Results/Chem22_CL_predictions.xlsx")


# PLOT 2 ---- scatter plot for CL (only 22chemical study: compare between 2D and Transwell assay)
#CL_Che22_dat <- rbind(Che22_Kidney_data[,c(1:5, 7:10)], Transwell_CL_data[,c(1:4, 6:10)])
CL_Che22_dat$Compound[which(CL_Che22_dat$Compound == "Gentamicin")] <- "Gentamycin"
CL_Che22_dat.rmna <- CL_Che22_dat[-which(is.na(CL_Che22_dat$CL_obs)),]
CL_Che22_dat.rmna <- CL_Che22_dat.rmna[-which(CL_Che22_dat.rmna$Study == "Lin et al."),]
CL_Che22_dat.rmna$diff <- log10(CL_Che22_dat.rmna$CL) - log10(CL_Che22_dat.rmna$CL_obs)
CL_Che22_dat.rmna$abs.diff <- abs(log10(CL_Che22_dat.rmna$CL) - log10(CL_Che22_dat.rmna$CL_obs))

stat <- ddply(CL_Che22_dat.rmna, .(Cell.type,Study,Invitro.type,Compound.type), summarize, 
              rho = cor(log10(CL_obs), log10(CL), use = "na.or.complete", method = "spearman"),
              r = cor(log10(CL_obs), log10(CL), use = "na.or.complete", method = "pearson"),
              MAE = mean(abs.diff), ME = mean(diff))
#write.xlsx(stat, "Results/Chem22_CL_stat.xlsx")

cor.dat <- ddply(CL_Che22_dat.rmna, .(Cell.type,Study,Invitro.type,Compound.type), summarize, cor = cor(log10(CL_obs), log10(CL), use = "na.or.complete", method = "spearman"))
cor.dat.p <- ddply(CL_Che22_dat.rmna, .(Cell.type,Study,Invitro.type,Compound.type), summarize, cor = cor(log10(CL_obs), log10(CL), use = "na.or.complete", method = "pearson"))
data2 <- data.frame(x = c(10^(seq(-7, 4, 1)), 1E+5), y = c(10^(seq(-7, 4, 1)), 1E+5))
data2$min <- data2$y/10
data2$max <- data2$y*10

ggplot(CL_Che22_dat.rmna, aes(x = CL_obs, y = CL))+
  geom_abline(intercept = 0, slope = 1, linetype=3, size = 0.35)+
  geom_ribbon(data=data2, aes(x=x,ymin=min,ymax=max),inherit.aes = FALSE, fill = "lightgray", alpha = 0.2) +
  #geom_abline(intercept = 1, slope = 1, linetype=1, size = 0.35, color = "lightgray")+
  #geom_abline(intercept = -1, slope = 1, linetype=1, size = 0.35, color = "lightgray")+
  geom_errorbar(aes(ymin = CL_LB, ymax = CL_UB, color = Compound), width = 0, size = 1) +
  geom_smooth(aes(linetype=Compound.type),method=lm, formula = y~x, se=FALSE, size = 0.8, color="black")+
  geom_point(aes(color = Compound, shape = Compound.type), fill = "white", size = 2.5) +
  geom_text(
    data    = cor.dat.p[which(cor.dat.p$Compound.type == "PFAS"),],
    mapping = aes(x=1E-3, y=5000, label = paste0("r[PFAS] ==", round(cor, digits = 2))), parse = TRUE,
    color = "black"
  )+
  geom_text(
    data    = cor.dat[which(cor.dat$Compound.type == "PFAS"),],
    mapping = aes(x=1E-3, y=500, label = paste0("rho[PFAS] ==", round(cor, digits = 2))), parse = TRUE,
    color = "black"
  )+
  
  geom_text(
    data    = cor.dat.p[which(cor.dat.p$Compound.type == "Non-PFAS"),],
    mapping = aes(x=1E+2, y=1E-3, label = paste0("r[Non-PFAS] ==", round(cor, digits = 2))), parse = TRUE,
    color = "black"
  )+
  geom_text(
    data    = cor.dat[which(cor.dat$Compound.type == "Non-PFAS"),],
    mapping = aes(x=1E+2, y=1E-4, label = paste0("rho[Non-PFAS] ==", round(cor, digits = 2))), parse = TRUE,
    color = "black"
  )+
  xlab("Observed renal clearance (mL/Kg per day)") +
  ylab("Predicted renal clearance (mL/Kg per day)") +
  facet_nested(Cell.type~Invitro.type, nest_line = element_line(linetype = 1), render_empty = F)+
  scale_shape_manual(values = c(24, 16))+
  coord_cartesian(xlim = c(1E-6, 9E+4), ylim = c(1E-6, 9E+4))+
  scale_x_log10(expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10^x, n = 6),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(expand = c(0,0),
                breaks = trans_breaks("log10", function(x) 10^x, n = 6),
                labels = trans_format("log10", math_format(10^.x))) +

  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=12),
        title = element_text(color = "black", size=12),
        plot.title = element_text(color = "black", face = "bold", size=14),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 0.6),
        #legend.position = "none",
        #legend.title = element_blank(),
        legend.text = element_text(color = "black", size=10),
        #legend.margin = margin(-0.5,0,0,0, unit="cm"),
        #legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        #legend.key.size = unit(0.4, 'cm'),
        strip.text.x = element_text(color = "black", size=12),
        strip.text.y = element_text(color = "black", size=12, angle = 0),
        strip.background = element_blank(),
        axis.ticks = element_line(color = "black", size = 0.3),
        axis.title = element_text(color = "black", size=12),
        axis.text = element_text(color = "black", size=9))
ggsave("CL_scatterplot_updated.jpeg",width = 13, height = 6.3)
ggsave("CL_scatterplot_updated.pdf",width = 13, height = 6.3)
write.xlsx(CL_Che22_dat.rmna, "Results/Chem22_CL_scatterplot_data_updated.xlsx")


# PLOT 3 ---- scatter plot for CL (compare among 22 chemicals study, Lucie's non-target analysis and Lin et al. (2024) study)
CL_Che22_dat <- rbind(Che22_Kidney_data[,c(1:5, 7:10)], Transwell_CL_data[,c(1:4, 6:10)])
cor.dat <- ddply(CL_Che22_dat, .(Cell.type,Study,Invitro.type,Compound.type), summarize, cor = cor(log10(CL_obs), log10(CL), use = "na.or.complete", method = "spearman"))
cor.dat.p <- ddply(CL_Che22_dat, .(Cell.type,Study,Invitro.type,Compound.type), summarize, cor = cor(log10(CL_obs), log10(CL), use = "na.or.complete", method = "pearson"))
ggplot(CL_Che22_dat, aes(x = CL_obs, y = CL))+
  geom_abline(intercept = 0, slope = 1, linetype=2, size = 0.35)+
  geom_abline(intercept = 1, slope = 1, linetype=2, size = 0.35, color = "lightgray")+
  geom_abline(intercept = -1, slope = 1, linetype=2, size = 0.35, color = "lightgray")+
  geom_errorbar(aes(ymin = CL_LB, ymax = CL_UB, color = Compound), width = 0, size = 1) +
  geom_point(aes(color = Compound, shape = Compound.type), fill = "white", size = 2.5) +
  geom_text(
    data    = cor.dat.p,
    mapping = aes(x=1E-6, y=1E+3, label = paste0("r ==", round(cor, digits = 2))), parse = TRUE,
    color = "black"
  )+
  geom_text(
    data    = cor.dat,
    mapping = aes(x=1E-6, y=1E+2, label = paste0("rho ==", round(cor, digits = 2))), parse = TRUE,
    color = "black"
  )+
  xlab("Observed renal clearance (mL/Kg per day)") +
  ylab("Predicted renal clearance (mL/Kg per day)") +
  facet_nested(Cell.type~Study + Compound.type + Invitro.type, nest_line = element_line(linetype = 1))+
  scale_shape_manual(values = c(21, 16))+
  scale_x_log10(lim = c(1E-8, 1E+4),
                breaks = trans_breaks("log10", function(x) 10^x, n = 6),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(lim = c(1E-8, 1E+4),
                breaks = trans_breaks("log10", function(x) 10^x, n = 6),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000", size=12),
        title = element_text(color = "black", size=12),
        plot.title = element_text(color = "black", face = "bold", size=14),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 0.6),
        #legend.position = "none",
        #legend.title = element_blank(),
        legend.text = element_text(color = "black", size=10),
        #legend.margin = margin(-0.5,0,0,0, unit="cm"),
        #legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        #legend.key.size = unit(0.4, 'cm'),
        strip.text.x = element_text(color = "black", size=12),
        strip.text.y = element_text(color = "black", size=12, angle = 0),
        strip.background = element_blank(),
        axis.ticks = element_line(color = "black", size = 0.3),
        axis.title = element_text(color = "black", size=12),
        axis.text = element_text(color = "black", size=9))
ggsave("CL_scatter_Comparemultiplestudy.jpeg",width = 16, height = 5.24)
ggsave("CL_scatter_Comparemultiplestudy.pdf",width = 16, height = 5.24)
