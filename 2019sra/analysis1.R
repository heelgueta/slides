##SETUP
#setup
options(scipen=999) #disables scientific notation
options(max.print=1000000) #enable long outputs

#load data
df402 <- read.csv("egic3.csv")
#subset data for study egic402
df402 <- df402[which(df402$source==402),]
#data imputation using mice
df402 <- mice::mice(df402, m=1, maxit=5, method = 'pmm', seed= 500); summary(df402)
df402 <- mice::complete(df402,1)
#setup for some factor variables
df402$expsrc <- as.factor(df402$expsrc)
df402$expval <- as.factor(df402$expval)
df402$expsrc <- factor(df402$expsrc, levels=c(-1,0,1),labels=c("Anti-Establishment","None","Establishment"))
df402$expval <- factor(df402$expval, levels=c(-1,0,1),labels=c("Anti-hydro","Neutral","Pro-hydro"))

#internal consistency analyses for some scales
round(psych::alpha(df402[,02:05])$total$std.alpha,3) #sich, alpha = .718
round(psych::alpha(df402[,10:13])$total$std.alpha,3) #simi, alpha = .749
round(psych::alpha(df402[,c(14:19,25:30)])$total$std.alpha,3) #trus, alpha = .874
round(psych::alpha(df402[,33:35])$total$std.alpha,3) #ttex, alpha = .936
round(psych::alpha(df402[,36:43])$total$std.alpha,3) #evhy, alpha = .899

#new variables: factor scores
cfamod <- 'sich =~ sich1d + sich2d + sich3d + sich4r';fitmod <- lavaan::cfa(cfamod, data=df402,estimator="MLR");df402$sich <- c(lavaan::lavPredict(fitmod))
cfamod <- 'simi =~ simi1d + simi2d + simi3r + simi4d';fitmod <- lavaan::cfa(cfamod, data=df402,estimator="MLR");df402$simi <- c(lavaan::lavPredict(fitmod))
cfamod <- 'trus =~ trugo1 + trugo2 + trugo3 + truco1';fitmod <- lavaan::cfa(cfamod, data=df402,estimator="MLR");df402$trus <- c(lavaan::lavPredict(fitmod))
cfamod <- 'ttex =~ ttext1 + ttext2 + ttext3';fitmod <- lavaan::cfa(cfamod, data=df402,estimator="MLR");df402$ttex <- c(lavaan::lavPredict(fitmod))
cfamod <- 'evhy =~ ehydr5 + ehydr6 + ehydr7 + ehydr8';fitmod <- lavaan::cfa(cfamod, data=df402,estimator="MLR");df402$evhy <- c(lavaan::lavPredict(fitmod))

#checking distribution of new variables (factor scores)
hist(df402$sich)
hist(df402$simi)
hist(df402$trus)
hist(df402$ttex)
hist(df402$evhy)

#variable labels
df402 = expss::apply_labels(df402,
                            partid = "Participant ID",
                            sich1d = "",
                            sich2d = "",
                            sich3d = "",
                            sich4r = "",
                            nati1d = "",
                            nati2r = "",
                            nati3d = "",
                            nati4d = "",
                            simi1d = "",
                            simi2d = "",
                            simi3r = "",
                            simi4d = "",
                            trutm1 = "",
                            trutm2 = "",
                            trutm3 = "",
                            trunm1 = "",
                            trunm2 = "",
                            trunm3 = "",
                            atnucl = "Pre-test attitudes Nucl power",
                            atther = "Pre-test attitudes Ther power",
                            athydr = "Pre-test attitudes Hydr power",
                            atgeot = "Pre-test attitudes Geot power",
                            atwind = "Pre-test attitudes Wind power",
                            trugo1 = "",
                            trugo2 = "",
                            trugo3 = "",
                            truco1 = "",
                            truco2 = "",
                            truco3 = "",
                            expval = "Exp. condition: Valence",
                            expsrc = "Exp. condition: Source",
                            ttext1 = "",
                            ttext2 = "",
                            ttext3 = "",
                            ehydr1 = "",
                            ehydr2 = "",
                            ehydr3 = "",
                            ehydr4 = "",
                            ehydr5 = "",
                            ehydr6 = "",
                            ehydr7 = "",
                            ehydr8 = "",
                            acnucl = "Post-test acceptance Nucl power",
                            acther = "Post-test acceptance Ther power",
                            achydr = "Post-test acceptance Hydr power",
                            acgeot = "Post-test acceptance Geot power",
                            acwind = "Post-test acceptance Wind power",
                            gender = "",
                            ethnic = "",
                            ageyea = "",
                            edupar = "",
                            polori = "",
                            source = "",
                            sich = "Degree of identification with the Majority",
                            simi = "Degree of identification with the Minorities",
                            trus = "General trust",
                            ttex = "Perceived trustworthiness",
                            evhy = "Post-test evaluation of Hydr power"
)

##ANALYSES
##cfa & sems using 'lavaan' package
cfamod <- '
sich =~ sich1d + sich3d + sich4r # + sich2d
simi =~ simi1d + simi2d + simi3r # + simi4d
trus =~ trugo1 + trugo2 + trugo3 # + truco1
evhy =~ ehydr6 + ehydr7 + ehydr8 #ehydr5 + 
#achydr ~  evhy + trus + sich + simi
evhy =~ trus + sich + simi
trus ~  sich + simi
sich ~~ simi
'

cfafit <- lavaan::cfa(model = cfamod, data = df402, estimator="MLR")
lavaan::fitmeasures(cfafit, fit.measures = c("chisq.scaled","df.scaled","pvalue.scaled","cfi.scaled","rmsea.scaled","srmr"))
#lavaan::standardizedsolution(cfafit)
#lavaan::modificationindices(cfafit, sort.=TRUE, minimum.value = 4, maximum.number = 20, standardized = TRUE)

#cfa diagram using 'semPlot' package
#semPlot::semPaths(cfafit)

##ancova model for analyses of experiment using 'car' package
aovmod <- 'ttex ~ sich + simi + trus + athydr +
                  expval + expsrc + 
                  expval * expsrc +
                  expval * athydr +
                  expsrc * simi +
                  expsrc * sich '
anofit <- lm(aovmod, data=df402)
knitr::kable(sjstats::anova_stats(car::Anova(anofit, type="III"))[2:12,c(1,4:5,7,6)])

#simple slopes analysis using 'emmeans' package
knitr::kable(emmeans::emtrends(anofit, pairwise ~ expsrc, var="sich"))

#plotting interactions using 'sjPlot' package
sjPlot::plot_model(anofit,type="pred",terms=c("trus"))
sjPlot::plot_model(anofit,type="pred",terms=c("athydr","expval"),line.size=3)
sjPlot::plot_model(anofit,type="pred",terms=c("sich","expsrc"),
                   colors = "Set1",line.size = 3,
                   title = "Interaction analysis",
                   legend.title = "Newspaper",
                   axis.title = c("Degree of identification with Majority","Perceived trustworthiness"))
