# R code for replicating the dataset and the analyses in the paper #
# Psychedelics as an intervention for psychological, existential distress in terminally ill patients: a systematic review and network meta-analysis #
# R code by Mattia Marchi (mattiamarchimd@gmail.com) 
# March 29, 2024

###----------------------------------------------------------------------------------------
###----------Pairwise Meta Analysis - Psychedelics in terminally ill patients--------------
###----------------------------------------------------------------------------------------
#Load required packages
library(meta)
library(tidyverse)
library(metafor)
library(netmeta)

###------------------------------1. Depression--------------------------------------------
#Import data
dep <- structure(list(ID = c("Fan et al, 2017", "Gasser et al, 2014", "Griffiths et al, 2016", "Grob et al, 2011", "Holze et al, 2023", 
                             "Liu et al, 2020", "Ross et al, 2016", "Wolfson et al, 2020", "Xu et al, 2017"),
                      mean1_depre = c(25.09, 7.5, 6.64, 10, 12.4, 8.21, 6.5, 9, 13.45), sd1_depre = c(7.07, 3.3, 5.2, 9.35, 12.2, 3.13, 7.03, 9, 5.21),
                      ncont1_depre = c(20L, 8L, 25L, 12L, 42L, 203L, 14L, 13L, 25L), mean2_depre = c(32.03, 8.7, 14.8, 16.1, 19, 11, 14.24, 12.2, 17.36),
                      sd2_depre = c(7.21, 2.9, 7.25, 12.47, 8, 3.8, 7.13, 5.3, 6.25), ncont2_depre = c(17L, 3L, 25L, 12L, 42L, 100L, 15L, 5L, 25L),
                      t1 = c("Ketamine", "LSD", "Psilocybin", "Psilocybin", "LSD", "Ketamine", "Psilocybin", "MDMA", "Ketamine"), t2 = c("Midazolam", "PBO", "PBO", "PBO", "PBO", "PBO", "PBO", "PBO", "PBO"),
                      contrast = c("Ketamine - Midazolam", "LSD - PBO", "Psilocybin - PBO", "Psilocybin - PBO", "LSD - PBO", "Ketamine - PBO", "Psilocybin - PBO", "MDMA - PBO", "Ketamine - PBO"),
                      Study.design = c("RCT", "RCT", "Crossover RCT", "Crossover RCT", "Crossover RCT", "RCT", "Crossover RCT", "RCT", "RCT"),
                      Country = c("China", "Switzerland", "USA", "USA", "Switzerland", "China", "USA", "USA", "China"),
                      Diagnosis = c("Cancer", "Life-threatening illness (most cancer)", "Life-threatening cancer", "Life-threatening cancer (advanced stage)", "Life-threatening illness (most cancer)", "Breast cancer", 
                                    "Life-threatening cancer", "Life-threatening illness (most cancer)", "Breast cancer"),
                      Setting = c("Inpatient clinic", "Psychotherapy", "Calm room", "Calm room", "Psychotherapy", "Intra-operative administration", "Psychotherapy", "Psychotherapy", "Intra-operative administration"),
                      Follow.up = c("1 and 3 days", "8 weeks", "5 and 26 weeks", "2 weeks", "2, 8, and 16 weeks", "3 days; 1, 4, and 12 weeks", "6 and 26 weeks", "4 weeks", "1 week"),
                      X..Females = c(25L, 36L, 49L, 92L, 48L, 100L, 62L, 78L, 100L), Mean.age = c(45.8, 51.7, 56.3, 47, 45, 47.4, 56.3, 54.9, 42.8), Std.mean.baseline.depression = c(3.63, 2.32, 4.88, 1.29, 2, 7.35, 1.91, 2.6, 6.25),
                      Comparator = c("Active", "Active", "Active", "Placebo", "Placebo", "Placebo", "Placebo", "Placebo", "Placebo")), row.names = c(NA, -9L), class = "data.frame")
#Random-effects meta-analysis of depression level
d <- metacont(n.e = ncont1_depre, mean.e = mean1_depre, sd.e = sd1_depre,
              n.c = ncont2_depre, mean.c = mean2_depre, sd.c = sd2_depre,
              data = dep, studlab = ID, sm = "SMD", method.tau = "DL")
d
forest(d, layout = "RevMan5", digits.sd = 2, random = T, fixed = F,
       label.e = "Psychedelics", label.c = "Controls",
       label.left = "Favours psychedelics", label.right = "Favours controls", allstudies = F)
#Publication Bias
#Funnel Plot
funnel(d, xlab = "Hedges' g")
#Egger's test
d_re_data <- escalc(measure = "SMD",
                    m1i = mean1_depre, sd1i = sd1_depre, n1i = ncont1_depre,
                    m2i = mean2_depre, sd2i = sd2_depre, n2i = ncont2_depre,
                    data = dep)
#Pooling ES
d_re <- rma(yi = d_re_data$yi, vi = d_re_data$vi)
d_re
regtest(d_re)
#Carry out trim-and-fill analysis
d_taf <- trimfill(d_re)
d_taf
#Leave-one-out
dep_loo <- leave1out(d_re)
#Meta-regression
dep$Study.design <- as.factor(dep$Study.design)
dep$Country <- as.factor(dep$Country)
dep$Setting <- as.factor(dep$Setting)
dep$Comparator <- as.factor(dep$Comparator)
#Age
mreg_d_age <- metareg(d, ~ Mean.age)
mreg_d_age
#% Female
mreg_d_fem <- metareg(d, ~ X..Females)
mreg_d_fem
#Depre baseline
mreg_d_depbaseline <- metareg(d, ~ Std.mean.baseline.depression)
mreg_d_depbaseline
#Setting
mreg_d_setting <- metareg(d, ~ Setting)
mreg_d_setting
#Study design
mreg_d_design <- metareg(d, ~ Study.design)
mreg_d_design
#Country
mreg_d_country <- metareg(d, ~ Country)
mreg_d_country

###----------------------------------2. Anxiety-----------------------------------------------
#Import data
anx <- structure(list(ID = c("Gasser et al, 2014", "Griffiths et al, 2016", "Holze et al, 2023", "Ross et al, 2016", "Wolfson et al, 2020"),
                      Timepoints = c("2 months", "5 weeks (reported) and 6 months", "2, 8, 16 weeks", "6 weeks", "1 month"),
                      mean1_anx = c(41.5, 34.64, 89.3, 33.5, 38.9), sd1_anx = c(3.2, 9.2, 37.7, 10.59, 10.6), ncont1_anx = c(8L, 25L, 42L, 14L, 13L),
                      mean2_anx = c(51.7, 40.48, 111, 47.05, 48.6), sd2_anx = c(5.3, 10.55, 20, 10.84, 12.6), ncont2_anx = c(3L, 25L, 42L, 15L, 5L),
                      t1 = c("LSD", "Psilocybin", "LSD", "Psilocybin", "MDMA"), t2 = c("PBO", "PBO", "PBO", "PBO", "PBO"),
                      contrast = c("LSD - LSD placebo", "Psilocybin - Psilocibyn placebo", "LSD - Placebo", "Psilocybin - PBO", "MDMA - Placebo"),
                      Study.design = c("RCT", "Crossover RCT", "Crossover RCT", "Crossover RCT", "RCT"), Country = c("Switzerland", "USA", "Switzerland", "USA", "USA"),
                      Diagnosis = c("Life-threatening illness (most cancer)", "Life-threatening cancer", "Life-threatening illness (most cancer)", "Life-threatening cancer", "Life-threatening illness (most cancer)"),
                      Setting = c("Psychotherapy", "Calm room", "Psychotherapy", "Psychotherapy", "Psychotherapy"), X..Females = c(36L, 49L, 48L, 62L, 78L),
                      Mean.age = c(51.7, 56.3, 45, 56.3, 54.9), Std.mean.baseline.anxiety = c(8.9, 5.39, 5.55, 4, 5.27), Comparator = c("Active", "Active", "Placebo", "Placebo", "Placebo")), row.names = c(NA, -5L), class = "data.frame")
#Random-effects meta-analysis of anxiety level
a <- metacont(n.e = ncont1_anx, mean.e = mean1_anx, sd.e = sd1_anx,
              n.c = ncont2_anx, mean.c = mean2_anx, sd.c = sd2_anx,
              data = anx, studlab = ID, sm = "SMD", method.tau = "DL")
a
forest(a, layout = "RevMan5", digits.sd = 2, random = T, fixed = F,
       label.e = "Psychedelics", label.c = "Controls",
       label.left = "Favours psychedelics", label.right = "Favours controls", allstudies = F)
#Publication Bias
#Funnel Plot
funnel(a, xlab = "Hedges' g")
#Egger's test
a_re_data <- escalc(measure = "SMD",
                    m1i = mean1_anx, sd1i = sd1_anx, n1i = ncont1_anx,
                    m2i = mean2_anx, sd2i = sd2_anx, n2i = ncont2_anx,
                    data = anx)
#Pooling ES
a_re <- rma(yi = a_re_data$yi, vi = a_re_data$vi)
d_re
regtest(a_re)
#Carry out trim-and-fill analysis
a_taf <- trimfill(a_re)
a_taf
#Leave-one-out
anx_loo <- leave1out(a_re)
#Meta-regression
anx$Study.design <- as.factor(anx$Study.design)
anx$Country <- as.factor(anx$Country)
anx$Setting <- as.factor(anx$Setting)
anx$Comparator <- as.factor(anx$Comparator)
#Age
mreg_a_age <- metareg(a, ~ Mean.age)
mreg_a_age
#% Female
mreg_a_fem <- metareg(a, ~ X..Females)
mreg_a_fem
#Anxiety baseline
mreg_a_anxbaseline <- metareg(a, ~ Std.mean.baseline.anxiety)
mreg_a_anxbaseline
#Setting
mreg_a_setting <- metareg(a, ~ Setting)
mreg_a_setting
#Study design
mreg_a_design <- metareg(a, ~ Study.design)
mreg_a_design
#Country
mreg_a_country <- metareg(a, ~ Country)
mreg_a_country


###----------------------------------------------------------------------------------------
###-----------Network Meta Analysis - Psychedelics in terminally ill patients--------------
###----------------------------------------------------------------------------------------

###----------------------------------1. Depression-----------------------------------------
#Organize data/calculate pairwise comparisons
pw1 <- pairwise(treat = list(t1, t2), n = list(ncont1_depre, ncont2_depre),
                mean = list(mean1_depre, mean2_depre), sd = list(sd1_depre, sd2_depre),
                studlab = ID, data = dep, sm = "SMD")
#Perform standard NMA
net1 <- netmeta(TE, seTE, treat1, treat2, studlab, data = pw1,
                common = FALSE, ref = "PBO")
dep_net <- netgraph(net1, seq = "optimal", plastic = FALSE, multiarm = TRUE,
                    cex.points = 6, number.of.studies = TRUE, cex.number = 1,
                    pos.number.of.studies = 0.3, iterate = F)
print(summary(net1))
forest(net1, label.left = "Favours active drug", label.right = "Favours placebo", sortvar = TE,
       smlab = "Std. Mean Difference
       Random, 95% CI")
plot(net1)
#To make all the contrasts
plot(net1, ref = c("K", "L", "MD", "Mid", "Psi", "PB"))
#Ranking treatments
netrank(net1)
set.seed(1909)
ran1 <- rankogram(net1, small.values = "good")
plot(ran1)#, cumulative.rankprob = T)
set.seed(1909)
sucra1 <- netrank(ran1)
netleague(net1, seq = netrank(net1), ci = T)
#Decompose heterogenity
decomp.design(net1)
netsplit(net1)

###------------------------------------2. Anxiety--------------------------------------------
#Organize data/calculate pairwise comparisons
pw2 <- pairwise(treat = list(t1, t2), n = list(ncont1_anx, ncont2_anx),
                mean = list(mean1_anx, mean2_anx), sd = list(sd1_anx, sd2_anx),
                studlab = ID, data = anx, sm = "SMD")
#Perform standard NMA
net2 <- netmeta(TE, seTE, treat1, treat2, studlab, data = pw2,
                common = FALSE, ref = "PBO")
anx_net <- netgraph(net2, seq = "optimal", plastic = FALSE, multiarm = TRUE,
                    cex.points = 6, number.of.studies = TRUE, cex.number = 1,
                    pos.number.of.studies = 0.3, iterate = F)
print(summary(net2))
forest(net2, label.left = "Favours active drug", label.right = "Favours placebo", sortvar = TE,
       smlab = "Std. Mean Difference
       Random, 95% CI")
plot(net2)
#To make all the contrasts
plot(net2, ref = c("L", "M", "Psi", "PB"))
#Ranking treatments
netrank(net2)
set.seed(1909)
ran2 <- rankogram(net2, small.values = "good")
plot(ran2)
set.seed(1909)
sucra2 <- netrank(ran2)
netleague(net2, seq = netrank(net2), ci = F)
#Decompose heterogenity
decomp.design(net2)
netsplit(net2)

###------------------------------------3. Death Acceptance----------------------------------------
#Import data
da <- structure(list(ID = c("Griffiths et al, 2016", "Ross et al, 2016", "Wolfson et al, 2020"),
                     mean1_deathacc = c(36.17, 13.64, 3.5), sd1_deathacc = c(7.95, 22.30, 1.6), ncont1_deathacc = c(25L, 14L, 13L),
                     mean2_deathacc = c(29.14, 11.27, 3), sd2_deathacc = c(11.25, 21.96, 0.7), ncont2_deathacc = c(25L, 15L, 5L),
                     t1 = c("Psilocybin", "Psilocybin", "MDMA"), t2 = c("PBO", "PBO", "PBO"), contrast = c("Psilocybin-PBO", "Psilocybin-PBO", "MDMA-PBO"),
                     note = c("death acceptance", "death transcendence (also death anxiety available)", "death acceptance (also fear of death available)")), row.names = c(NA, -3L), class = "data.frame")
#Organize data/calculate pairwise comparisons
pw6 <- pairwise(treat = list(t1, t2), n = list(ncont1_deathacc, ncont2_deathacc),
                mean = list(mean1_deathacc, mean2_deathacc), sd = list(sd1_deathacc, sd2_deathacc),
                studlab = ID, data = da, sm = "SMD")
#Perform standard NMA
net6 <- netmeta(TE, seTE, treat1, treat2, studlab, data = pw6,
                common = FALSE, ref = "PBO")
da_net <- netgraph(net6, seq = "optimal", plastic = FALSE, multiarm = TRUE,
                   cex.points = 6, number.of.studies = TRUE, cex.number = 1,
                   pos.number.of.studies = 0.3, iterate = F)
print(summary(net6))
forest(net6, label.left = "Favours placebo", label.right = "Favours active drug", sortvar = TE,
       smlab = "Std. Mean Difference
       Random, 95% CI")
plot(net6)
#To make all the contrasts
plot(net6, ref = c("M", "Psi", "PB"))
#Ranking treatments
netrank(net6)
set.seed(1909)
ran6 <- rankogram(net6, small.values = "bad")
plot(ran6)
set.seed(1909)
sucra6 <- netrank(ran6)
netleague(net6, seq = netrank(net6), ci = F)
#Decompose heterogenity
decomp.design(net6)
netsplit(net6)

###--------------------------------4. Quality of life------------------------------------
#Import data
qol <- structure(list(ID = c("Gasser et al, 2014", "Griffiths et al, 2016", "Ross et al, 2016", "Wolfson et al, 2020"),
                      mean1_QoL = c(50, 7.14, 15.43, 15), sd1_QoL = c(14.9, 1.45, 2.47, 3.9), ncont1_QoL = c(8L, 25L, 14L, 13L),
                      mean2_QoL = c(36, 6.17, 12.28, 16.3), sd2_QoL = c(12.8, 1.6, 2.56, 6.7), ncont2_QoL = c(3L, 25L, 15L, 5L),
                      t1 = c("LSD", "Psilocybin", "Psilocybin", "MDMA"), t2 = c("PBO", "PBO", "PBO", "PBO"),
                      contrast = c("LSD-PBO", "Psilocybin-PBO", "Psilocybin-PBO", "MDMA-PBO")), row.names = c(NA, -4L), class = "data.frame")
#Organize data/calculate pairwise comparisons
pw7 <- pairwise(treat = list(t1, t2), n = list(ncont1_QoL, ncont2_QoL),
                mean = list(mean1_QoL, mean2_QoL), sd = list(sd1_QoL, sd2_QoL),
                studlab = ID, data = qol, sm = "SMD")
#Perform standard NMA
net7 <- netmeta(TE, seTE, treat1, treat2, studlab, data = pw7,
                common = FALSE, ref = "PBO")
qol_net <- netgraph(net7, seq = "optimal", plastic = FALSE, multiarm = TRUE,
                    cex.points = 6, number.of.studies = TRUE, cex.number = 1,
                    pos.number.of.studies = 0.3, iterate = F)
print(summary(net7))
forest(net7, label.left = "Favours placebo", label.right = "Favours active drug", sortvar = TE,
       smlab = "Std. Mean Difference
       Random, 95% CI")
plot(net7)
#To make all the contrasts
plot(net7, ref = c("L", "MD", "Psi", "PB"))
#Ranking treatments
netrank(net7)
set.seed(1909)
ran7 <- rankogram(net7, small.values = "bad")
plot(ran7)
set.seed(1909)
sucra7 <- netrank(ran7)
netleague(net7, seq = netrank(net7), ci = F)
#Decompose heterogenity
decomp.design(net7)
netsplit(net7)

###---------------------------5. Safety and Tolerability-------------------------------
#Import data
tol <- structure(list(Author..year = c("Fan et al, 2017", "Gasser et al, 2014", "Griffiths et al, 2016", "Grob et al, 2011", "Holze et al, 2023", "Liu et al, 2020", "Ross et al, 2016", "Wolfson et al, 2020", "Xu et al, 2017"),
                      n1 = c(20L, 8L, 56L, 12L, 42L, 203L, 29L, 13L, 25L), n2 = c(17L, 3L, 56L, 12L, 42L, 10L, 29L, 5L, 25L),
                      t1 = c("Ketamine", "LSD", "Psilocybin", "Psilocybin", "LSD", "Ketamine", "Psilocybin", "MDMA", "Ketamine"), t2 = c("Midazolam", "PBO", "PBO", "PBO", "PBO", "PBO", "PBO", "PBO", "PBO"),
                      contrast = c("Ketamine-Midazolam", "LSD-LSD PBO", "Psilocybin-Psilocybin PBO", "Psilocybin-PBO", "LSD-PBO", "Ketamine-PBO", "Psilocybin-PBO", "MDMA-PBO", "Ketamine-PBO"),
                      N_Drop.out_any.cause_t = c(0L, 0L, 5L, 0L, 2L, 0L, 3L, 1L, 0L), N_Drop.out_any.cause_c = c(0L, 0L, 5L, 0L, 4L, 0L, 4L, 0L, 0L), N_Drop.out_severe.adverse.reactions_t = c(0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L),
                      N_Drop.out_severe.adverse.reactions_c = c(0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L), N_Death_t = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L), N_Death_c = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                      N_Adverse.effects_t = c(0L, 6L, 10L, NA, 6L, 72L, 1L, 12L, 14L), N_Adverse.effects_c = c(0L, 0L, 23L, NA, 3L, 40L, 0L, 2L, 17L)), class = "data.frame", row.names = c(NA, -9L))

#----A. Dropout any cause
pw3 <- pairwise(treat = list(t1, t2), n = list(n1, n2),
                event = list(N_Drop.out_any.cause_t, N_Drop.out_any.cause_c),
                studlab = Author..year, data = tol, sm = "OR")
net3 <- netmeta(TE, seTE, treat1, treat2, studlab, data = pw3,
                common = FALSE, ref = "PBO")
anycause_net <- netgraph(net3, seq = "optimal", plastic = FALSE, multiarm = TRUE,
                         cex.points = 6, number.of.studies = TRUE, cex.number = 1,
                         pos.number.of.studies = 0.3, iterate = F)
print(summary(net3))
forest(net3, label.left = "Favours active drug", label.right = "Favours placebo", sortvar = TE,
       smlab = "Odds Ratio
       Random, 95% CI")
plot(net3)
#To make all the contrasts
plot(net3, ref = c("L", "M", "Psi", "PB"))
#Ranking treatments
netrank(net3)
set.seed(1909)
ran3 <- rankogram(net3, small.values = "good")
plot(ran3)
set.seed(1909)
sucra3 <- netrank(ran3)
netleague(net3, seq = netrank(net3), ci = F)
#Decompose heterogenity
decomp.design(net3)
netsplit(net3)

#----B. N adverse effects
pw5 <- pairwise(treat = list(t1, t2), n = list(n1, n2),
                event = list(N_Adverse.effects_t, N_Adverse.effects_c),
                studlab = Author..year, data = tol, sm = "OR")
net5 <- netmeta(TE, seTE, treat1, treat2, studlab, data = pw5,
                common = FALSE, ref = "PBO")
ae_net <- netgraph(net5, seq = "optimal", plastic = FALSE, multiarm = TRUE,
                   cex.points = 6, number.of.studies = TRUE, cex.number = 1,
                   pos.number.of.studies = 0.3, iterate = F)
print(summary(net5))
forest(net5, label.left = "Favours active drug", label.right = "Favours placebo", sortvar = TE,
       smlab = "Odds Ratio
       Random, 95% CI")
plot(net5)
#To make all the contrasts
plot(net5, ref = c("K" ,"L", "M", "Psi", "PB"))
#Ranking treatments
netrank(net5, small.values = "good")
set.seed(1909)
ran5 <- rankogram(net5, small.values = "good")
plot(ran5)
set.seed(1909)
sucra5 <- netrank(ran5)
netleague(net5, seq = netrank(net5), ci = F)
#Decompose heterogenity
decomp.design(net5)
netsplit(net5)