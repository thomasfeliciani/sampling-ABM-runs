

################################################################################
# Cleaning environment, loading resources.

rm(list = ls())

#library("lhs")
#library("ParamHelpers") ##### supposedly works natively with factors?
library("clhs") # for Latin hypercube sampling
library("AlgDesign") # for generating optimal designs
library("skpr")
library("lubridate") # for handling time strings
library("fpc")#("philentropy") # for computing Bhattacharyya distances
library("reshape2")
library("ggplot2")
library("openxlsx") # for exporting results tables
library("compiler") # for speeding up Monte Carlo experiments.

source("./util.r")


# Importing factorial design from NetLogo output
d1 <- read.csv(
  file = 
    "./input/SegregationExtended_TF sampling approaches - scenario 1-table.csv",
  skip = 6, header = TRUE, sep = ","
)
d2 <- read.csv(
  file = 
    "./input/SegregationExtended_TF sampling approaches - scenario 2-table.csv",
  skip = 6, header = TRUE, sep = ","
)
d1$scenario <- 1
d2$scenario <- 2

d <- rbind(d1, d2)
rm(d1,d2)

################################################################################
# Cleaning / formatting data.


# Removing redundant / unnecessary variables
d$BlackWhiteVis <- NULL
d$X.step. <- NULL
d$maxTicks <- NULL
d$X.run.number. <- NULL
d$noise <- NULL
d$stepsBetweenUpdates <- NULL

# Relabeling
varnames <- c(
  "density",
  "simWanted",
  "nbRadius", ###
  "fourGroups",
  "seed",
  "pop",
  "popRed",
  "popGreen",
  "popBlue",
  "popYellow",
  "percRed",
  "percGreen",
  "percBlue",
  "percYellow",
  "percSimilarT0",
  "percSimilarRedT0",
  "percSimilarGreenT0",
  "percSimilarBlueT0",
  "percSimilarYellowT0",
  "percUnhappyT0",
  "percUnhappyRedT0",
  "percUnhappyGreenT0",
  "percUnhappyBlueT0",
  "percUnhappyYellowT0",
  "percSimilar",
  "percSimilarRed",
  "percSimilarGreen",
  "percSimilarBlue",
  "percSimilarYellow",
  "percUnhappy",
  "percUnhappyRed",
  "percUnhappyGreen",
  "percUnhappyBlue",
  "percUnhappyYellow",
  "percClustRed",
  "percClustGreen",
  "percClustBlue",
  "percClustYellow",
  "time",
  "CPUtime",
  "scenario"
)
# View(cbind(names(d), varnames))
names(d) <- varnames; rm(varnames)

# Some variables need to be recoded:
d$nGroups <- dplyr::recode(
  .x = d$fourGroups,
  "true" = 4,
  "false" = 2
)
d$nGroups <- factor(
  x = d$nGroups,
  levels = c(2, 4),
  labels = c("two-group condition", "four-group condition")
)
# table(d$nGroups, d$fourGroups)
d$fourGroups <- NULL

# reverse coding %unhappy into contentedness:
d$contentednessT0 <- 100 - d$percUnhappyT0
d$contentedness <- 100 - d$percUnhappy

# recoding % same color neighbors into a proportion:
d$propSameT0 <- d$percSimilarT0 / 100
d$propSame <- d$percSimilar / 100

# Lastly, dealing with the convergence condition:
d$convergence <- factor(
  x = d$time == 500,
  levels = c(TRUE, FALSE),
  labels = c("not converged by t=500", "converged by t=500")
)

# Ensuring that each variable has the proper class:
for (v in 1:ncol(d)){
  print("__________________________________________")
  print(paste0(
    "[", v, "/", ncol(d), "] ", names(d)[v]
  ))
  print(paste("class:", class(d[,v])))
  nlevels <- length(unique(d[,v]))
  ifelse(nlevels <= 5, print(table(d[,v])), print(paste("levels:", nlevels)))
  print(summary(d[,v]))
}
# All seems alright.




################################################################################
# Drawing samples

# Taking care of reproducibility:
set.seed(100)

parameters <- c(# Selecting our indep vars (here: all initialization parameters)
  "density",
  "simWanted", 
  "nbRadius", 
  "nGroups"
)
measuresT0 <- c( # Configuration of the initialized environment (time = 0)
  "pop",
  #"percRed", "percGreen", "percBlue", "percYellow",
  "propSameT0"#,
  #"contentednessT0" ########### 
)
measuresTfinal <- c( # Model outcomes (configuration at time = final)
  "propSame",
  "contentedness",
  "time",
  "convergence"
)

# Identifying unique parameter configurations for the two scenarios:
paramSpace1 <- base::unique(d[d$scenario == 1, parameters])
paramSpace2 <- base::unique(d[d$scenario == 2, parameters])
paramSpace1$paramConfigID <- paramSpace1$designConfigID <- 1:nrow(paramSpace1)
paramSpace2$paramConfigID <- paramSpace2$designConfigID <- 1:nrow(paramSpace2)

# Determining which unique run belongs to which of the unique parameter configs:
d$paramConfigID <- NA
d1 <- d[d$scenario == 1,]
d2 <- d[d$scenario == 2,]

for (i in 1:nrow(paramSpace1)) { # Scenario 1
  belongs <- d1$density == paramSpace1$density[i] &
    d1$simWanted == paramSpace1$simWanted[i] &
    d1$nbRadius == paramSpace1$nbRadius[i] &
    d1$nGroups == paramSpace1$nGroups[i]
  
  d1$paramConfigID[belongs] <- i
}
d1$designConfigID <- d1$paramConfigID

for (i in 1:nrow(paramSpace2)) { # Scenario 2
  belongs <- d2$density == paramSpace2$density[i] &
    d2$simWanted == paramSpace2$simWanted[i] &
    d2$nbRadius == paramSpace2$nbRadius[i] &
    d2$nGroups == paramSpace2$nGroups[i]
  
  d2$paramConfigID[belongs] <- i
}
d2$designConfigID <- d2$paramConfigID

################################################################################
# Step 1: Selecting parameter configurations


#d <- drawSamples(d = d, designSpace = paramSpace)









################################################################################
################################################################################
################################################################################
# Design evaluation
#
# Here is a function that takes a design (e.g. "d[d_a,]"), a model to estimate 
# and outputs the performance metrics we need.


evalDesign <- function(
    data,
    designSpace = designSpace,
    paramSpace = designSpace,
    sizeExperiment = 15,
    design = paste0("design_", letters[1:9]),
    includeMultiRuns = FALSE,
    resample = 10, # how many times each sampling design is to be tested
    replacement = TRUE, # For OD designs; ignored in all other cases
    model#,
    #returnSampledDataset = TRUE
) {
  
  # First we establish our ground truth: the regression model estimated using
  # all available simulation runs (in jargon: the ground truth is defined as 
  # the multi-run full factorial design).
  gtm <- lm(data = data, formula = model)
  
  
  # Now let's draw the samples and calculate their performance relative to the
  # ground truth:
  results <- list()
  
  pb <- txtProgressBar(min = 0, max = resample, style = 3)
  for(t in 1:resample) {
    data <- drawSamples(
      d = data,
      designSpace = designSpace,
      paramSpace = paramSpace,
      designs = design,
      includeMultiRuns = includeMultiRuns,
      replacement = replacement,
      sizeExperiment = sizeExperiment
    )
    
    results[[t]] <- list(
      data = data,
      resultsTable = data.frame(
        design = character(),#deparse(substitute(design)), #design name
        size = numeric(), # how many runs
        sizeProp = numeric(), # proportion of factorial
        CPUtimeS = numeric(), # CPU time in seconds
        CPUtime = character(),# CPU time (as a formatted, human-readable string)
        R2 = numeric(),
        bhattacharyya = numeric(),
        model = character()
      ),
      model = list()
    )
    
    for (i in 1:length(design)) {
      dd <- data[data[,design[i]],] #selecting runs from the specified design
      
      # Calculating performance metrics:
      #   - CPU time:
      CPUtime <- dd$CPUtime |>
        sum() |> round() |> lubridate::seconds_to_period() |> as.character()
      
      #   - R squared: 
      m <- lm(data = dd, formula = model)
      R2 <- summary(m)$r.squared
      
      #   - Bhattacharyya distance from the multi-run full-factorial (design_am)
      bhattacharyya <- fpc::bhattacharyya.dist(
        mu1 = gtm$coefficients,
        mu2 = m$coefficients,
        Sigma1 = vcov(gtm),
        Sigma2 = vcov(m)
      )
      
      results[[t]]$model[[design[i]]] <- m
      results[[t]]$resultsTable[i,] <- data.frame(
        design = design[i],#deparse(substitute(design)), #design name
        size = nrow(dd), # how many runs
        sizeProp = sizeExperiment,
        CPUtimeS = sum(dd$CPUtime), # CPU time in seconds
        CPUtime = CPUtime,# CPU time (as a formatted, human-readable string)
        R2 = R2,
        bhattacharyya = bhattacharyya,
        model = deparse(model)
      )
    }
    
    setTxtProgressBar(pb, t)
  }
  
  for (i in 1:resample) ifelse(
    i == 1,
    results$summary <- results[[i]]$resultsTable,
    results$summary <- rbind(results$summary, results[[i]]$resultsTable)
  )
  
  close(pb)
  return(results)
}


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# First step - Sampling parameter configurations



# Ground truth _________________________________________________________________
png(
  filename = "./output/Figure 2 - sampling A linear.png",
  height = 1150, width = 2000, res = 300
)
ggplot(
  data = d1,
  aes(
    x = nbRadius,#as.factor(nbRadius),
    y = propSame,
    fill = nGroups,
    color = nGroups
  )
) +
  #geom_jitter(aes(shape = convergence), width = 0.18, height = 0, alpha = 0.2)+
  geom_violin(
    data = d1[d1$nGroups == "four-group condition",],
    aes(group = nbRadius, color = nGroups), 
    position = position_nudge(x = 0.08),
    scale = "width", alpha = 0.5
    #color = "blue", fill = "blue"
  ) +
  #geom_jitter(
  #  data = d1[d1$nGroups == "four-group condition",],
  #  aes(x = nbRadius + 0.08), 
  #  width = 0.07, height = 0, size = 0.1, alpha = 0.15) +
  geom_boxplot(
    data = d1[d1$nGroups == "four-group condition",],
    aes(group = nbRadius, color = nGroups),
    position = position_nudge(x = 0.08),
    alpha = 0, color = "gray10", width = 0.08
  ) +
  
  geom_violin(
    data = d1[d1$nGroups == "two-group condition",],
    aes(group = nbRadius, color = nGroups), 
    position = position_nudge(x = -0.08),
    scale = "width", alpha = 0.5
  ) +
  #geom_jitter(
  #  data = d1[d1$nGroups == "two-group condition",],
  #  aes(x = nbRadius - 0.08), 
  #  width = 0.07, height = 0, size = 0.1, alpha = 0.15) +
  geom_boxplot(
    data = d1[d1$nGroups == "two-group condition",],
    aes(group = nbRadius, color = nGroups),
    position = position_nudge(x = -0.08),
    alpha = 0, color = "gray10", width = 0.08
  ) +
  
  scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.8, direction = -1) +
  scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8, direction = -1) +
  #scale_shape_manual(
  #  values = c("converged by t=500" = 1, "not converged by t=500" = 4)
  #) +
  geom_smooth(method = "lm", level = 0.99, linewidth = 0.5, fullrange = TRUE) +
  labs(
    x = "neighborhood radius",
    #y = "average proportion of same-color neighbors\n(t=final)",
    y = "segregation (t=final)",
    color = "color group composition",
    fill = "color group composition",
    #shape = ""
  ) +
  guides(
    color = guide_legend(reverse = TRUE),
    fill = guide_legend(reverse = TRUE)
  ) +
  scale_y_continuous(
    #limits = c(0, 1),
    expand = c(0,0), 
    breaks = seq(from = 0, to = 1, by = 0.05)
  ) +
  scale_x_continuous(limits = c(1, 8.5), breaks = unique(d$nbRadius)) +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid.major = element_line(color = "#F5F5F5", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA, color = NA)
  )
dev.off()



# Samples ______________________________________________________________________
#
# Taking the samples and evaluating them:
#designsToTest <- paste0("design_", letters[1:9]) # designs A-I
designsToTest <- paste0("design_", letters[1:5]) # designs A-E
sizeParamSpace <- nrow(paramSpace1)

compiler::enableJIT(1)
r20 <- evalDesign( # one fifth
  data = d1,
  designSpace = paramSpace1,
  paramSpace = paramSpace1,
  sizeExperiment = 16,#round(sizeParamSpace * 0.15),
  design = designsToTest,
  resample = 100,
  includeMultiRuns = FALSE,
  replacement = FALSE,
  model = propSame ~ nbRadius * nGroups
)
r25 <- evalDesign( # one fourth
  data = d1,
  designSpace = paramSpace1,
  paramSpace = paramSpace1,
  sizeExperiment = 20,#round(sizeParamSpace * 0.15),
  design = designsToTest,
  resample = 100,
  includeMultiRuns = FALSE,
  replacement = FALSE,
  model = propSame ~ nbRadius * nGroups
)
r50 <- evalDesign( # one half
  data = d1,
  designSpace = paramSpace1,
  paramSpace = paramSpace1,
  sizeExperiment = 40,#round(sizeParamSpace * 0.15),
  design = designsToTest,
  resample = 100,
  includeMultiRuns = FALSE,
  replacement = FALSE,
  model = propSame ~ nbRadius * nGroups
)
r100 <- evalDesign( # same size as factorial
  data = d1,
  designSpace = paramSpace1,
  paramSpace = paramSpace1,
  sizeExperiment = 80,#round(sizeParamSpace * 0.15),
  design = designsToTest,
  resample = 100,
  includeMultiRuns = FALSE,
  replacement = FALSE,
  model = propSame ~ nbRadius * nGroups
)
compiler::enableJIT(0)


summary <- rbind(
  r20$summary,
  r25$summary[r25$summary$design != "design_a",],
  r50$summary[r50$summary$design != "design_a",],
  r100$summary[r100$summary$design != "design_a",]
)


# Keeping tally of how many times each run was chosen: #########################

for (rsize in c(20, 25, 50, 100)) { # For every sample size...
  print(paste("Tallying sampling events for designs of size", rsize))
  for (design in letters[1:5]) { # ...and or every deisgn...
    
    # ... initialize the tally by setting each run to 0 sampling events...
    d1[,paste0("choice_design_", design, "_r", rsize)] <- 0 
    
    # ... and then add 1 to it each time the run was chosen:
    for (element in 1:nrow(paramSpace1)) {
      d1[,paste0("choice_design_", design, "_r", rsize)] <- 
        d1[,paste0("choice_design_", design, "_r", rsize)] + 
        get(paste0("r", rsize))[[element]]$data[,paste0("design_", design)]
      #r25[[element]]$data$design_a
    }
    
    # Also saving the same info in a more convenient data format:
    temp <- data.frame(
      seed = d1$seed,
      scenario = 1,
      size = rsize,
      design = paste0("design_", design),
      tally = d1[, paste0("choice_design_", design, "_r", rsize)]
    )
    ifelse(
      exists("tallies"),
      tallies <- rbind(tallies, temp),
      tallies <- temp
    )
    rm(temp)
    
  }
}

# Creating a saving point ______________________________________________________
# Explorting output. Writing to file may tike some time.
save(tallies, summary, r20, r25, r50, r100, file = "./output/scenario1.RData")

# May resume here:
#load("./output/scenario1.RData")

# Plotting tallies.
# First I recode the variables:
tallies$design <- factor(
  x = tallies$design,
  levels = designsToTest,
  labels = c(
    "a) factorial", 
    "b) uniform sampling", 
    "c) LHS", 
    "d) D-OD\n(linear + unif.sampl.)", 
    "e) D-OD\n(quadratic)"
  )
)
tallies$sizeLabel <- factor(
  x = tallies$size, 
  levels = c(20, 25, 50, 100),#c(1, 0.5, 0.25, 0.15),
  labels = c("N=16 (20%)", "N=20 (25%)", "N=40 (50%)", "N=80 (100%)")
)

png(
  filename = "./output/Figure tbd - sampling tallies in scenario 1.png",
  height = 1500, width = 2500, res = 300
)
ggplot(
  data = tallies[!(tallies$design == "a) factorial" & tallies$size < 100),], 
  aes(x = tally, fill = sizeLabel)
) +
  geom_bar() +
  facet_grid(rows = vars(sizeLabel), cols = vars(design)) +
  scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.8) +
  scale_x_continuous(breaks = 0:max(tallies$tally)) +
  scale_y_continuous(breaks = 0) +
  labs(y = "frequency", x = "times sampled") +
  theme(
    plot.background = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.line.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    #axis.text.x = element_text(size = 7),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    #axis.title.x = element_blank(),
    strip.text = element_text(size = 10),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "#F5F5F5", linewidth = 0.5),
    #panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size = 10),
    legend.position = "NA"
  )
dev.off()



summary$design <- factor(
  x = summary$design,
  levels = designsToTest,
  labels = c(
    "a) factorial", 
    "b) uniform sampling", 
    "c) LHS", 
    "d) D-OD (linear + unif.sampl.)", 
    "e) D-OD (quadratic)"#,
    #"f) A-OD (linear + unif.sampl.)",
    #"g) A-OD (quadratic)",
    #"h) I-OD (linear + unif.sampl.)",
    #"i) I-OD quadratic"
  )
)
summary$CPUtimeS <- summary$CPUtimeS / 60 # converting to minutes

#summary$sizeLabel[summary$design == "a) factorial"] <- 1
summary$sizeLabel <- factor(
  x = summary$size, 
  levels = c(80, 40, 20, 16),#c(1, 0.5, 0.25, 0.15),
  labels = c("N=80 (100%)", "N=40 (50%)", "N=20 (25%)", "N=16 (20%)")
)

# For readability:
#summary$bhattacharyya_log <- log(summary$bhattacharyya, base = 2)
summary$bhattacharyya_log <- log(summary$bhattacharyya, base = 10)



png(
  filename = "./output/Figure 3 - performance with model 1.png",
  height = 1500, width = 2500, res = 300
)
ggplot(
  data = reshape2::melt(
    summary, measure.vars = c("CPUtimeS", "bhattacharyya_log")), 
  aes(x = design, y = value, color = sizeLabel)
) + 
  geom_violin(
    aes(color = sizeLabel, fill = sizeLabel),
    position = position_dodge(0.5), alpha = 0.5, scale = "width"
  ) +
  #geom_point(
  #  position = position_jitterdodge(dodge.width = 0.5),
  #  shape = 16, alpha = 0.6, size = 0.2
  #) +
  geom_boxplot(
    aes(fill = sizeLabel),
    position = position_dodge(0.5),
    width = 0.3, color = "black", alpha = 0, outlier.shape = NA
  ) +
  scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +
  scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.8) +
  facet_wrap(
    ~ variable, 
    scales = "free_y", 
    labeller = labeller(variable = c(
      "CPUtimeS" = "CPU time in minutes",
      #"bhattacharyya" = "log\u2082 (Bhattacharyya distance)"
      "bhattacharyya_log" = "log10(Bhattacharyya distance)"
    ))
  ) +
  labs(
    color = "Number of trials\n(% of factorial design)", 
    fill = "Number of trials\n(% of factorial design)"
  ) +
  theme(
    plot.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
    axis.title = element_blank(),
    strip.text = element_text(size = 11),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "#F5F5F5", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size = 10),
    legend.position = "top"
  )
dev.off()      




# Showing which runs were sampled ______________________________________________
# preparing data structure for visualizing all samples side-by-side

run <- 100 # which run do we look at

for (s in 1:length(designsToTest)) {
  temp <- r20[[run]]$data[
    r20[[run]]$data[,designsToTest[s]],  # rows we want
    c(parameters, measuresT0, measuresTfinal) # variables we want
  ]
  temp$design <- designsToTest[s]
  ifelse(s == 1, rs <- temp, rs <- rbind(rs, temp))
}
rs$design <- factor(
  x = rs$design,
  levels = designsToTest,
  labels = c(
    "a) factorial\nN=80", 
    "b) uniform sampling\nN=16", 
    "c) LHS\nN=16", 
    "d) D-OD (linear + unif.sampl.)\nN=16", 
    "e) D-OD (quadratic)\nN=16"#,
    #"f) A-OD (linear + unif.sampl.)\nN=12",
    #"g) A-OD (quadratic)\nN=12",
    #"h) I-OD (linear + unif.sampl.)\nN=12",
    #"i) I-OD quadratic\nN=12"
  )
)


png(
  filename = "./output/Figure 4 - comparing samples.png",
  #height = 1850, 
  height = 1400, 
  width = 2200, res = 300
)
ggplot(
  data = rs[order(rs$nGroups),], # ordering ensures that 4groups appear on top
  aes(x = nbRadius, y = density, color = nGroups)
) +
  #geom_jitter(width = 0.1, height = 0.1, size = 2.5, alpha = 0.6) +
  geom_point(aes(shape = nGroups), size = 3, alpha = 1) +
  scale_shape_manual(values = c(
    "two-group condition" = 2,
    "four-group condition" = 0
  )) +
  facet_wrap(~design) +
  scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +
  labs(
    x = "neighborhood radius",
    y = "density", fill = "", shape = "", color = ""
  ) +
  scale_x_continuous(breaks = unique(rs$nbRadius)) +
  scale_y_continuous(breaks = unique(d$density), labels = addSuffix) +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "#F5F5F5", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA, color = NA),
    legend.position = c(0.85, 0.15)
  )
dev.off()


# Showing fitted models ________________________________________________________
png(
  filename = "./output/Figure 5 - fitted models.png",
  #height = 1850, 
  height = 1400, 
  width = 2200, res = 300
)
ggplot(
  data = rs,
  #data = d,
  aes(
    x = nbRadius,
    y = propSame,
    fill = nGroups,
    color = nGroups
  )
) +
  #geom_jitter(aes(shape = convergence), width = 0.2, height = 0) +
  #geom_point(aes(shape = convergence), size = 2.5) +
  geom_point(shape = 1, size = 2.5) +
  facet_wrap(~design) + 
  geom_smooth(method = "lm", level = 0.95, linewidth = 0.5, fullrange = TRUE) +
  scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.8) +
  scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +
  #scale_shape_manual(
  #  values = c("converged by t=500" = 1, "not converged by t=500" = 4)
  #) +
  labs(
    x = "neighborhood radius",
    #y = "average proportion of same-color neighbors\n(t=final)",
    y = "segregation (t=final)",
    color = "",
    fill = "",
    shape = ""
  ) +
  scale_y_continuous(
    #limits = c(0, 1),
    expand = c(0,0)#, 
    #breaks = seq(from = 0, to = 1, by = 0.1)
  ) +
  scale_x_continuous(breaks = unique(d$nbRadius)) +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "#F5F5F5", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA),
    #legend.position = "bottom"
    legend.position = c(0.85, 0.15)
  )
dev.off()

# Exporting regression tables __________________________________________________
#
# Here I retrieve, format and export the regression table from the
# experiment results shown in the previous figures.
benchmark <- as.data.frame(summary(lm(
  data = d1,
  formula = r20$summary$model[1]#propSame ~ nbRadius * nGroups
))$coefficients)
benchmark <- data.frame(
  coefficient = rownames(benchmark),
  estimate = round(benchmark$Estimate, digits = 3),
  pvalue = ptostars(x = benchmark$"Pr(>|t|)"), # significance stars
  stderr = round(benchmark$"Std. Error", digits = 3)
)
rtable <- data.frame(Coefficient = benchmark$coefficient)
rtable[,2] <- paste0(
  as.character(benchmark$estimate), benchmark$pvalue,
  " (", as.character(benchmark$stderr), ")"
)
names(rtable)[2] <- "ground truth"

for (s in 1:length(designsToTest)) {
  temp <- summary(r20[[run]]$model[[designsToTest[s]]])$coefficients
  temp <- as.data.frame(temp)
  temp <- data.frame(
    coefficient = rownames(temp),
    estimate = round(temp$Estimate, digits = 3),
    pvalue = ptostars(x = temp$"Pr(>|t|)"), # significance stars
    stderr = round(temp$"Std. Error", digits = 3)
  )
  
  #if (s == 1) rtable <- data.frame(Coefficient = temp$coefficient)
  rtable[,2 + s] <- paste0(
    as.character(temp$estimate), temp$pvalue,
    " (", as.character(temp$stderr), ")"
  )
  names(rtable)[2 + s] <- designsToTest[s]
}

openxlsx::write.xlsx(x = rtable, file = "./output/Table 5.4.xlsx") 






# Appendix heatmap _____________________________________________________________
# The heatmap shows which areas of the design space are preferred by each
# samplign approach. For this visualization, too, we assume a design size of 15%

# For each point in the design space we tally how many times it was picked by
# each samplig approach:
rs <- r20[[1]]$data

ds <- expand.grid(
  density = unique(rs$density),
  nbRadius = unique(rs$nbRadius)#,
  #nGroups = unique(rs$nGroups) ##
)
ds[,designsToTest] <- 0

rs[,designsToTest] <- 0
for (s in 1:length(designsToTest)) for (run in 1:10) {#foreach design and sample
  rs[,designsToTest[s]] <- 
    rs[,designsToTest[s]] + r20[[run]]$data[,designsToTest[s]]
}

for (i in 1:nrow(ds)) {
  rows <- subset(
    x = rs, 
    rs$density == ds$density[i] & rs$nbRadius == ds$nbRadius[i]
  )
  for (s in 1:length(designsToTest)) ds[i, designsToTest[s]] <- sum(
    rows[,designsToTest[s]]
  )
}

ds <- reshape2::melt(data = ds, id.vars = c("density", "nbRadius"))#,"nGroups"))
ds$design <- factor(
  x = ds$variable,
  levels = designsToTest,
  labels = c(
    "a) factorial\nN=80 per experiment", 
    "b) uniform sampling\nN=16 per experiment", 
    "c) LHS\nN=16 per experiment", 
    "d) D-OD (linear + unif.sampl.)\nN=16 per experiment", 
    "e) D-OD (quadratic)\nN=16 per experiment"#,
    #"f) A-OD (linear + unif.sampl.)\nN=12 per experiment",
    #"g) A-OD (quadratic)\nN=12 per experiment",
    #"h) I-OD (linear + unif.sampl.)\nN=12 per experiment",
    #"i) I-OD quadratic\nN=12 per experiment"
  )
)
ds$variable <- NULL
ds$nbRadius <- as.factor(ds$nbRadius)

ds$value <- ds$value / max(ds$value)

png(
  filename = "./output/Figure A1 - heatmaps first scenario.png",
  #height = 1850, 
  height = 1400, 
  width = 2200, res = 300
)
ggplot(data = ds, aes(x = nbRadius, y = density, fill = value)) +
  geom_tile() +
  facet_wrap(~design) +
  scale_fill_viridis_c(
    option = "C", begin = 0.1, end = 0.8,
    breaks = c(0, 0.5, 1),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      #ticks = FALSE,
      barheight = 0.7
      #nbin = length(unique(ds$value))
    )
  ) +
  labs(
    x = "neighborhood radius",
    y = "density", 
    fill = "relative frequency"
  ) +
  scale_x_discrete(expand = c(0,0), breaks = unique(rs$nbRadius)) +
  scale_y_continuous(
    expand = c(0,0),
    breaks = unique(ds$density),
    labels = addSuffix
  ) +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "#F5F5F5", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA),
    #legend.direction = "horizontal",
    legend.position = c(0.85, 0.15) #"bottom"
  )
dev.off()












################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Second scenario - Sampling initialized environments


# Ground truth _________________________________________________________________

png(
  filename = "./output/Figure 6 - sampling B non-linear.png",
  #height = 1300, width = 1300, res = 300
  height = 1150, width = 2000, res = 300
)
ggplot(data = d2, aes(x = simWanted, y = propSame)) +
  geom_violin(
    aes(group = as.factor(simWanted)),
    scale = "width", color = "gray83", fill = "gray83"
  ) +
  geom_boxplot(
    aes(group = as.factor(simWanted)),
    alpha = 0, color = "gray35", width = 1
  ) +
  geom_smooth(
    aes(x = simWanted, y = propSame),
    formula = y ~ x + poly(x, 2),
    level = 0.99,
    linewidth = 0.5,
    color = "black",
    fill = "gray70",
    fullrange = TRUE
  ) +
  labs(
    x = "tolerance",
    y = "segregation (t=final)"
  ) +
  scale_y_continuous(
    breaks = seq(from = 0, to = 1, by = 0.1)
  ) +
  scale_x_continuous(
    breaks = seq(from = 0, to = 100, by = 5),
    labels = addSuffix
  ) +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "#F5F5F5", linewidth = 0.5),
    panel.grid.minor = element_blank(),
  )
dev.off()



# Samples ______________________________________________________________________
#
#
# First of all, we need to update our deisgn space. W
designsToTest <- paste0("design_", letters[1:5]) # designs A-E

#designSpace <- unique.data.frame(d2[,c(parameters, "contentednessT0")])
designSpace <- unique.data.frame(d2[,c(parameters)])

designSpace$designConfigID <- 1:nrow(designSpace)
for (i in 1:nrow(designSpace)) {
  belongs <- d2$density == designSpace$density[i] &
    d2$simWanted == designSpace$simWanted[i] &
    d2$nbRadius == designSpace$nbRadius[i] &
    d2$nGroups == designSpace$nGroups[i] #&
    #d2$contentednessT0 == designSpace$contentednessT0[i]
  
  d2$designConfigID[belongs] <- i
}
sizeParamSpace <- nrow(paramSpace2)



#sizeParamSpace <- length(table(d2$density)) * length(table(d2$simWanted))
#sizeParamSpace <- nrow(paramSpace[
#  paramSpace$nGroups == "two-group condition" & paramSpace$nbRadius == 1.5,
#])

# Taking the samples and evaluating them:

#designsToTest <- paste0("design_", letters[1:9]) # designs A-I

d2$QUADsimWanted <- d2$simWanted ^ 2
#d2$QUADcontentednessT0 <- d2$contentednessT0 ^2


# Taking the samples and evaluating them:
compiler::enableJIT(1)
r20 <- evalDesign(
  data = d2,
  designSpace = designSpace,
  paramSpace = paramSpace2,
  sizeExperiment = 20,#round(sizeParamSpace * 0.15),
  design = designsToTest,
  resample = 100,
  includeMultiRuns = FALSE,
  replacement = FALSE,     #replacement = TRUE
  model = propSame ~ simWanted + QUADsimWanted
)
r25 <- evalDesign(
  data = d2,
  designSpace = designSpace,
  paramSpace = paramSpace2,
  sizeExperiment = 25,
  design = designsToTest,
  resample = 100,
  includeMultiRuns = FALSE,
  replacement = FALSE,
  model = propSame ~ simWanted + QUADsimWanted
)
r50 <- evalDesign(
  data = d2,
  designSpace = designSpace,
  paramSpace = paramSpace2,
  sizeExperiment = 50,
  design = designsToTest,
  resample = 100,
  includeMultiRuns = FALSE,
  replacement = FALSE,
  model = propSame ~ simWanted + QUADsimWanted
)
r100 <- evalDesign(
  data = d2,
  designSpace = designSpace,
  paramSpace = paramSpace2,
  sizeExperiment = 100,
  design = designsToTest,
  resample = 100,
  includeMultiRuns = FALSE,
  replacement = FALSE,
  model = propSame ~ simWanted + QUADsimWanted
)
compiler::enableJIT(0)


summary <- rbind(
  r20$summary,
  r25$summary[r25$summary$design != "design_a",],
  r50$summary[r50$summary$design != "design_a",],
  r100$summary[r100$summary$design != "design_a",]
)



###################################
# Taking tallies
rm(tallies)
for (rsize in c(20, 25, 50, 100)) { # For every sample size...
  print(paste("Tallying sampling events for designs of size", rsize))
  for (design in letters[1:5]) { # ...and or every deisgn...
    
    # ... initialize the tally by setting each run to 0 sampling events...
    d2[,paste0("choice_design_", design, "_r", rsize)] <- 0 
    
    # ... and then add 1 to it each time the run was chosen:
    for (element in 1:nrow(paramSpace1)) {
      d2[,paste0("choice_design_", design, "_r", rsize)] <- 
        d2[,paste0("choice_design_", design, "_r", rsize)] + 
        get(paste0("r", rsize))[[element]]$data[,paste0("design_", design)]
      #r25[[element]]$data$design_a
    }
    
    # Also saving the same info in a more convenient data format:
    temp <- data.frame(
      seed = d2$seed,
      scenario = 2,
      size = rsize,
      design = paste0("design_", design),
      tally = d2[, paste0("choice_design_", design, "_r", rsize)]
    )
    ifelse(
      exists("tallies"),
      tallies <- rbind(tallies, temp),
      tallies <- temp
    )
    rm(temp)
    
  }
}


# Creating a saving point ______________________________________________________
save(tallies, summary, r20, r25, r50, r100, file = "./output/scenario2.RData")
#load("./output/scenario2.RData")




# Plotting tallies.
# First I recode the variables:
tallies$design <- factor(
  x = tallies$design,
  levels = designsToTest,
  labels = c(
    "a) factorial", 
    "b) uniform sampling", 
    "c) LHS", 
    "d) D-OD\n(linear + unif.sampl.)", 
    "e) D-OD\n(quadratic)"
  )
)
tallies$sizeLabel <- factor(
  x = tallies$size, 
  levels = c(20, 25, 50, 100),#c(1, 0.5, 0.25, 0.15),
  labels = c("N=20 (20%)", "N=25 (25%)", "N=50 (50%)", "N=100 (100%)")
)


png(
  filename = "./output/Figure tbd - sampling tallies in scenario 2.png",
  height = 1500, width = 2500, res = 300
)
ggplot(
  data = tallies[!(tallies$design == "a) factorial" & tallies$size < 100),], 
  aes(x = tally, fill = sizeLabel)
) +
  geom_bar() +
  facet_grid(rows = vars(sizeLabel), cols = vars(design)) +
  scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.8) +
  #scale_x_continuous(breaks = 0:max(tallies$tally)) +
  scale_x_continuous(breaks = seq(from = 0, to = max(tallies$tally), by = 2)) +
  scale_y_continuous(breaks = 0) +
  labs(y = "frequency", x = "times sampled") +
  theme(
    plot.background = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.line.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    #axis.text.x = element_text(size = 7),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    #axis.title.x = element_blank(),
    strip.text = element_text(size = 10),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "#F5F5F5", linewidth = 0.5),
    #panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size = 10),
    legend.position = "NA"
  )
dev.off()




summary$design <- factor(
  x = summary$design,
  levels = designsToTest,
  labels = c(
    "a) factorial", 
    "b) uniform sampling", 
    "c) LHS", 
    "d) D-OD (linear + unif.sampl.)", 
    "e) D-OD (quadratic)"#,
    #"f) A-OD (linear + unif.sampl.)",
    #"g) A-OD (quadratic)",
    #"h) I-OD (linear + unif.sampl.)",
    #"i) I-OD quadratic"
  )
)
summary$CPUtimeS <- summary$CPUtimeS / 60 # converting to minutes

#summary$sizeLabel[summary$design == "a) factorial"] <- 1
summary$sizeLabel <- factor(
  x = summary$size, 
  levels = c(100, 50, 25, 20),
  labels = c("100%, N=100", "50%, N=50", "25%, N=25", "20%, N=20")
)

# For readability:
summary$bhattacharyya_log <- log(summary$bhattacharyya, base = 10)


png(
  filename = "./output/Figure 7 - performance with model 2.png",
  height = 1500, width = 2500, res = 300
)
ggplot(
  data = reshape2::melt(
    summary, measure.vars = c("CPUtimeS", "bhattacharyya_log")), 
  aes(x = design, y = value, color = sizeLabel)
) + 
  geom_violin(
    aes(color = sizeLabel, fill = sizeLabel),
    position = position_dodge(0.5), alpha = 0.5, scale = "width"
  ) +
  #geom_point(
  #  position = position_jitterdodge(dodge.width = 0.5),
  #  shape = 16, alpha = 0.6, size = 0.2
  #) +
  geom_boxplot(
    aes(fill = sizeLabel),
    position = position_dodge(0.5),
    width = 0.3, color = "black", alpha = 0, outlier.shape = NA
  ) +
  scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +
  scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.8) +
  facet_wrap(
    ~ variable, 
    scales = "free_y", 
    labeller = labeller(variable = c(
      "CPUtimeS" = "CPU time in minutes",
      #"bhattacharyya" = "log\u2082 (Bhattacharyya distance)"
      "bhattacharyya_log" = "log10(Bhattacharyya distance)"
    ))
  ) +
  labs(
    color = "Number of trials\n(% of factorial design)", 
    fill = "Number of trials\n(% of factorial design)"
  ) +
  theme(
    plot.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
    axis.title = element_blank(),
    strip.text = element_text(size = 11),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "#F5F5F5", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size = 10),
    legend.position = "top"
  )
dev.off()      



# Showing which runs were sampled ______________________________________________
# preparing data structure for visualizing all samples side-by-side
for (s in 1:length(designsToTest)) {
  temp <- r20[[run]]$data[
    r20[[run]]$data[,designsToTest[s]],  # rows we want
    c(parameters, measuresT0, measuresTfinal) # variables we want
  ]
  temp$design <- designsToTest[s]
  ifelse(s == 1, rs <- temp, rs <- rbind(rs, temp))
}
rs$design <- factor(
  x = rs$design,
  levels = designsToTest,
  labels = c(
    "a) factorial\nN=100", 
    "b) uniform sampling\nN=20", 
    "c) LHS\nN=20", 
    "d) D-OD (linear + unif.sampl.)\nN=20", 
    "e) D-OD (quadratic)\nN=20"#,
    #"f) A-OD (linear + unif.sampl.)\nN=15",
    #"g) A-OD (quadratic)\nN=15",
    #"h) I-OD (linear + unif.sampl.)\nN=15",
    #"i) I-OD quadratic\nN=15"
  )
)

png(
  filename = "./output/Figure 8 - comparing samples.png",
  #height = 1850,
  height = 1400,
  width = 2200, res = 300
)
#ggplot(data = rs, aes(x = contentednessT0, y = density, color = simWanted)) + 
ggplot(data = rs, aes(x = simWanted, y = density)) + 
  #geom_jitter(width = 0, height = 0.15, size = 2.5, alpha = 0.6, shape = 1) +
  geom_point(size = 2.5, alpha = 0.6, shape = 1) +
  facet_wrap(~design) + 
  #scale_color_viridis_c(option = "C", begin = 0.2, end = 0.8) +
  labs(x = "tolerance", y = "density", color = "tolerance (%)") +
  scale_x_continuous(
    breaks = seq(from = 0, to = 100, by = 10),
    labels = addSuffix
  ) +
  scale_y_continuous(breaks = unique(d$density), labels = addSuffix) +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "#F5F5F5", linewidth = 0.5),
    panel.grid.minor = element_blank()#,
    #legend.background = element_blank(),
    #legend.key = element_rect(fill = NA),
    #legend.position = c(0.85, 0.15)
  )
dev.off()


# Showing fitted models ________________________________________________________
png(
  filename = "./output/Figure 9 - fitted models.png",
  #height = 1850,
  height = 1400,
  width = 2200, res = 300
)
ggplot(
  data = rs,
  aes(x = simWanted, y = propSame) #, fill = nGroups, color = nGroups
) +
  facet_wrap(~design) + 
  geom_smooth(
    method = "lm", 
    formula = y ~ x + poly(x, 2), 
    level = 0.99, 
    linewidth = 0.5,
    color = "gray40",
    fill = "gray80",
    fullrange = TRUE,
  ) +
  #geom_point(aes(shape = convergence), alpha = 0.6, size = 2.5) +
  geom_point(shape = 1, size = 2.5) +
  #scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.8) +
  #scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +
  #scale_shape_manual(
  #  values = c("converged by t=500" = 1, "not converged by t=500" = 4)
  #) +
  labs(
    x = "tolerance",
    #y = "average proportion of same-color neighbors\n(t=final)",
    y = "segregation (t=final)",
    #color = "",
    #fill = "",
    shape = ""
  ) +
  scale_y_continuous(
    #limits = c(0, 1),
    #expand = c(0,0), 
    breaks = seq(from = 0, to = 5, by = 0.1)
  ) +
  scale_x_continuous(
    breaks = seq(from = 0, to = 100, by = 10),
    labels = addSuffix
  ) +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "#F5F5F5", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA),
    #legend.position = "bottom"
    legend.position = c(0.85, 0.15)
  )
dev.off()

# Showing regression tables ____________________________________________________
#
# Here I retrieve, format and export the regression table from the
# experiment results shown in the previous figures.

benchmark <- as.data.frame(summary(lm(
  data = d2,
  formula = r20$summary$model[1]#propSame ~ nbRadius * nGroups
))$coefficients)
benchmark <- data.frame(
  coefficient = rownames(benchmark),
  estimate = round(benchmark$Estimate, digits = 5),
  pvalue = ptostars(x = benchmark$"Pr(>|t|)"), # significance stars
  stderr = round(benchmark$"Std. Error", digits = 5)
)
rtable <- data.frame(Coefficient = benchmark$coefficient)
rtable[,2] <- paste0(
  as.character(benchmark$estimate), benchmark$pvalue,
  " (", as.character(benchmark$stderr), ")"
)
names(rtable)[2] <- "ground truth"

for (s in 1:length(designsToTest)) {
  temp <- summary(r20[[run]]$model[[designsToTest[s]]])$coefficients
  temp <- as.data.frame(temp)
  temp <- data.frame(
    coefficient = rownames(temp),
    estimate = round(temp$Estimate, digits = 5),
    pvalue = ptostars(x = temp$"Pr(>|t|)"), # significance stars
    stderr = round(temp$"Std. Error", digits = 5)
  )
  
  #if (s == 1) rtable <- data.frame(Coefficient = temp$coefficient)
  rtable[,2 + s] <- paste0(
    as.character(temp$estimate),
    temp$pvalue,
    "(",
    as.character(temp$stderr),
    ")"
  )
  names(rtable)[2 + s] <- designsToTest[s]
}

openxlsx::write.xlsx(x = rtable, file = "./output/Table 5.6.xlsx")



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

# Appendix 1 ___________________________________________________________________
# The heatmap shows which areas of the design space are preferred by each
# sampling approach. For this visualization, too, we assume a design size of 20%

# For each point in the design space we tally how many times it was picked by
# each samplig approach:


## Scenario 2 __________________________________________________________________
load("./output/scenario2.RData")

rs <- r20[[1]]$data

# Discretizing the continuous variable "contentednessT0"
rs$contentednessT0 <- plyr::round_any(rs$contentednessT0, 5, f = round)
rs20 <- r20
for (run in 1:10) {
  rs20[[run]]$data$contentednessT0 <- plyr::round_any(
    rs20[[run]]$data$contentednessT0, 5, f = round
  )
  
}


ds <- expand.grid(
  density = unique(rs$density),
  simWanted = 1:20 * 5#,
  #contentednessT0 = 0:20 * 5 #,
  #nGroups = unique(rs$nGroups) ##
)
ds[,designsToTest] <- 0

rs[,designsToTest] <- 0
for (s in 1:length(designsToTest)) for (run in 1:100) {#foreach design and sample
  rs[,designsToTest[s]] <- 
    rs[,designsToTest[s]] + rs20[[run]]$data[,designsToTest[s]]
}

for (i in 1:nrow(ds)) {
  rows <- subset(
    x = rs, 
    rs$density == ds$density[i] & rs$simWanted == ds$simWanted[i]
    #rs$density == ds$density[i] & rs$contentednessT0 == ds$contentednessT0[i]
  )
  for (s in 1:length(designsToTest)) ds[i, designsToTest[s]] <- sum(
    rows[,designsToTest[s]]
  )
}

#ds <- reshape2::melt(data = ds, id.vars = c("density", "contentednessT0"))
ds <- reshape2::melt(data = ds, id.vars = c("density", "simWanted"))
ds$design <- factor(
  x = ds$variable,
  levels = designsToTest,
  labels = c(
    "a) factorial\nN=100 per experiment", 
    "b) uniform sampling\nN=20 per experiment", 
    "c) LHS\nN=20 per experiment", 
    "d) D-OD (linear + unif.sampl.)\nN=20 per experiment", 
    "e) D-OD (quadratic)\nN=20 per experiment"#,
    #"f) A-OD (linear + unif.sampl.)\nN=15 per experiment",
    #"g) A-OD (quadratic)\nN=15 per experiment",
    #"h) I-OD (linear + unif.sampl.)\nN=15 per experiment",
    #"i) I-OD quadratic\nN=15 per experiment"
  )
)
ds$variable <- NULL



png(
  filename = "./output/Figure A2 - heatmaps second scenario.png",
  #height = 1850, 
  height = 1400, 
  width = 2200, res = 300
)
#ggplot(data = ds, aes(x = contentednessT0, y = density, fill = value)) +
ggplot(data = ds, aes(x = simWanted, y = density, fill = value)) +
  geom_tile()+#geom_bin_2d() +#geom_tile() +#geom_bin_2d() +
  facet_wrap(~design) +
  scale_fill_viridis_c(
    option = "C", begin = 0.05, end = 0.82,
    #breaks = c(0, 0.5, 1),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barheight = 0.7
    )
  ) +
  labs(
    x = "tolerance",#"contentedness (t=0)",
    y = "density", 
    fill = "frequency"
  ) +
  scale_x_continuous(
    expand = c(0,0),
    breaks = seq(from = 0, to = 100, by=10),
    labels = addSuffix
  ) +
  scale_y_continuous(
    expand = c(0,0),
    breaks = unique(ds$density),
    labels = addSuffix
  ) +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "#F5F5F5", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA),
    #legend.direction = "horizontal",
    legend.position = c(0.85, 0.15) #"bottom"
  )
dev.off()


