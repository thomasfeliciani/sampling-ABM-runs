#  Utility functions

# Extract samples from a dataset according to the various sampling approaches:

drawSamples <- function(
    d, 
    designSpace,
    paramSpace = designSpace,
    designs = paste0("design_", letters[1:9]),# designs A-I where A is factorial
    includeMultiRuns = FALSE,
    sizeExperiment = 15 # proportion of a factorial design.
) {
  
  ##############################################################################
  # Step 1: Selecting parameter configurations
  
  # a) full factorial __________________________________________________________
  # Starting with one run per condition:
  if ("design_a" %in% designs) {
    d_a <- numeric()
    for (i in 1:nrow(paramSpace)) d_a[i] <- sample(
      which(d$paramConfigID == paramSpace$paramConfigID[i]),
      size = 1
    )
    d$design_a <- FALSE
    d$design_a[d_a] <- TRUE
    
    # Then multiple runs per condition:
    if (includeMultiRuns) {
      d_am <- 1:nrow(d)
      d$design_am <- FALSE
      d$design_am[d_am] <- TRUE
    }
  }
  
  
  
  # b) uniform sampling ________________________________________________________
  if ("design_b" %in% designs) {
    # Selecting unique parameter configurations at random uniform:
    s <- sample(
      1:nrow(designSpace),
      size = sizeExperiment,
      replace = FALSE
    )
    
    # Next, we identify corresponding unique simulation runs from d, either
    # assuming one run per condition (i.e. one run per sample in s), or multiple
    # 
    # One run per condition:
    d_b <- numeric()
    for (i in 1:length(s)) {
      IDs <- which(d$designConfigID == designSpace$designConfigID[s[i]])
      ifelse(
        length(IDs) == 1,
        d_b[i] <- IDs,
        d_b[i] <- sample(IDs, size = 1)
      )
    }
    d$design_b <- FALSE
    d$design_b[d_b] <- TRUE
    
    # Multiple runs per condition:
    if (includeMultiRuns) {
      d_bm <- which(d$designConfigID %in% designSpace$designConfigID[s])
      d$design_bm <- FALSE
      d$design_bm[d_bm] <- TRUE
    }
  }
  
  
  
  # c) LHS _____________________________________________________________________
  # Latin hypercube sampling.
  # For this I use the library "clhs": it allows for factors / ordinal variables
  # as we have; plus, if needed, it allows us to associate a cost for each
  # experiment run.
  if("design_c" %in% designs) {
    s <- clhs::clhs(x = designSpace, size = sizeExperiment)
    
    # One run per condition
    d_c <- numeric()
    for (i in 1:length(s)) {
      IDs <- which(d$designConfigID == designSpace$designConfigID[s[i]])
      ifelse(
        length(IDs) == 1,
        d_c[i] <- IDs,
        d_c[i] <- sample(IDs, size = 1)
      )
    }
    d$design_c <- FALSE
    d$design_c[d_c] <- TRUE
    
    # Multiple runs per condition:
    if (includeMultiRuns) {
      d_cm <- which(d$designConfigID %in% designSpace$designConfigID[s])
      d$design_cm <- FALSE
      d$design_cm[d_cm] <- TRUE
    }
  }
  
  
  
  # Optimal designs ____________________________________________________________
  # Here we prepare the materials we need to generate optimal designs.
  #
  # We start by identifying the parameters on which we need to sample.
  # The first thing we need to do is to remove ID variables from the list
  # of parameters:
  toremove <- which(names(designSpace) %in% c("designConfigID","paramConfigID"))
  params <- names(designSpace)[- toremove]
  
  # Second, we ensure that we have been given variation in all parameters
  # (necessary to avoid so-called "singular" designs):
  #params <- c("density", "simWanted", "nbRadius", "nGroups")
  toremove <- c()
  for(i in 1:length(params)) if (length(unique(designSpace[,params[i]])) < 2) {
    toremove <- c(toremove, i)
  }
  if(!is.null(toremove)) params <- params[- toremove]
  
  # Let's define which parameters can be used as quadratic (they need at least
  # three levels):
  quadraticTerms <- c()
  for (par in params) if (length(unique(designSpace[,par])) > 2) {
    quadraticTerms <- c(quadraticTerms, par)
  }
  
  # Second: we can define a formula with all first-level interactions and all
  # valid quadratic terms
  quadraticFormula <- as.formula(paste0(
    "~ (.) ^ 2 + ",
    paste0("I(", quadraticTerms, " ^ 2)", collapse = " + ")
  ))
  # Now everything's ready for generating optimal designs:
  
  
  # d) D-OD linear + uniform sampling __________________________________________
  if ("design_d" %in% designs) {
    # Here we split the sample into two subsamples: one from an opimal D design
    # and one from a uniform sample, here called "urs".
    # The function AlgDesign::optFederov sort of handles this out of the box.
    # 
    s <- redoIfError(
      attempts = 5,
      expression = {
        urs <- sample( # 
          1:nrow(designSpace),
          size = floor(sizeExperiment / 2),
          replace = FALSE
        )
        
        # And then we draw an optimal D design, including the uniform sample:
        AlgDesign::optFederov(
          frml = ~ (.) ^ 2, # Linear effects, with all two-way interactions
          center = TRUE,
          nTrials = sizeExperiment,
          augment = TRUE,
          rows = urs, # this is where the uniform random sample is fed in
          data = designSpace[,params],
          criterion = "D" #####################
        )$rows
      }
    )
    
    # Now, our resulting samples for (e) D-OD + uniform sampling are:
    #
    # One run per condition:
    d_d <- numeric()
    for (i in 1:length(s)) {
      IDs <- which(d$designConfigID == designSpace$designConfigID[s[i]])
      ifelse(
        length(IDs) == 1,
        d_d[i] <- IDs,
        d_d[i] <- sample(IDs, size = 1)
      )
    }
    d$design_d <- FALSE
    d$design_d[d_d] <- TRUE
    
    # Multiple runs per condition:
    if (includeMultiRuns) {
      d_dm <- which(d$designConfigID %in% designSpace$designConfigID[s])
      d$design_dm <- FALSE
      d$design_dm[d_dm] <- TRUE
    }
  }
  
  
  
  # e) D-OD quadratic __________________________________________________________
  # An optimal D design, this time specified on a quadratic model.
  if("design_e" %in% designs){
    s <- AlgDesign::optFederov(
      frml = quadraticFormula,
      #frml = ~ (.) ^ 2,
      center = TRUE,
      nTrials = sizeExperiment,
      data = designSpace[,params],
      criterion = "D"
    )$rows
    
    # One run per condition
    d_e <- numeric()
    for (i in 1:length(s)) {
      IDs <- which(d$designConfigID == designSpace$designConfigID[s[i]])
      ifelse(
        length(IDs) == 1,
        d_e[i] <- IDs,
        d_e[i] <- sample(IDs, size = 1)
      )
    }
    d$design_e <- FALSE
    d$design_e[d_e] <- TRUE
    
    # Multiple runs per condition:
    if (includeMultiRuns) {
      d_em <- which(d$designConfigID %in% designSpace$designConfigID[s])
      d$design_em <- FALSE
      d$design_em[d_em] <- TRUE
    }
  }
  
  
  
  # f) A-OD linear + uniform sampling __________________________________________
  if ("design_f" %in% designs) {
    # Here we split the sample into two subsamples: one from an opimal D design
    # and one from a uniform sample, here called "urs".
    # The function AlgDesign::optFederov sort of handles this out of the box.
    # 
    s <- redoIfError(
      attempts = 5,
      expression = {
        urs <- sample( # 
          1:nrow(designSpace),
          size = floor(sizeExperiment / 2),
          replace = FALSE
        )
        
        # And then we draw an optimal D design, including the uniform sample:
        AlgDesign::optFederov(
          frml = ~ (.) ^ 2, # Linear effects, with all two-way interactions
          center = TRUE,
          nTrials = sizeExperiment,
          augment = TRUE,
          rows = urs, # this is where the uniform random sample is fed in
          data = designSpace[,params],
          criterion = "A" #####################
        )$rows
      }
    )
    
    # Now, our resulting samples for (e) D-OD + uniform sampling are:
    #
    # One run per condition:
    d_f <- numeric()
    for (i in 1:length(s)) {
      IDs <- which(d$designConfigID == designSpace$designConfigID[s[i]])
      ifelse(
        length(IDs) == 1,
        d_f[i] <- IDs,
        d_f[i] <- sample(IDs, size = 1)
      )
    }
    d$design_f <- FALSE
    d$design_f[d_f] <- TRUE
    
    # Multiple runs per condition:
    if (includeMultiRuns) {
      d_fm <- which(d$designConfigID %in% designSpace$designConfigID[s])
      d$design_fm <- FALSE
      d$design_fm[d_fm] <- TRUE
    }
  }
  
  
  
  # g) A-OD quadratic __________________________________________________________
  # An optimal D design, this time specified on a quadratic model.
  if("design_g" %in% designs){
    s <- AlgDesign::optFederov(
      frml = quadraticFormula, #################
      center = TRUE,
      nTrials = sizeExperiment,
      data = designSpace[,params],
      criterion = "A" ####################
    )$rows
    
    # One run per condition
    d_g <- numeric()
    for (i in 1:length(s)) {
      IDs <- which(d$designConfigID == designSpace$designConfigID[s[i]])
      ifelse(
        length(IDs) == 1,
        d_g[i] <- IDs,
        d_g[i] <- sample(IDs, size = 1)
      )
    }
    d$design_g <- FALSE
    d$design_g[d_g] <- TRUE
    
    # Multiple runs per condition:
    if (includeMultiRuns) {
      d_gm <- which(d$designConfigID %in% designSpace$designConfigID[s])
      d$design_gm <- FALSE
      d$design_gm[d_gm] <- TRUE
    }
  }
  
  
  
  # h) I-OD linear + uniform sampling __________________________________________
  if ("design_h" %in% designs) {
    # Here we split the sample into two subsamples: one from an opimal D design
    # and one from a uniform sample, here called "urs".
    # The function AlgDesign::optFederov sort of handles this out of the box.
    # 
    s <- redoIfError(
      attempts = 5,
      expression = {
        urs <- sample( # 
          1:nrow(designSpace),
          size = floor(sizeExperiment / 2),
          replace = FALSE
        )
        
        # And then we draw an optimal D design, including the uniform sample:
        AlgDesign::optFederov(
          frml = ~ (.) ^ 2, # Linear effects, with all two-way interactions
          center = TRUE,
          nTrials = sizeExperiment,
          augment = TRUE,
          rows = urs, # this is where the uniform random sample is fed in
          data = designSpace[,params],
          criterion = "I" #####################
        )$rows
      }
    )
    
    # Now, our resulting samples for (e) D-OD + uniform sampling are:
    #
    # One run per condition:
    d_h <- numeric()
    for (i in 1:length(s)) {
      IDs <- which(d$designConfigID == designSpace$designConfigID[s[i]])
      ifelse(
        length(IDs) == 1,
        d_h[i] <- IDs,
        d_h[i] <- sample(IDs, size = 1)
      )
    }
    d$design_h <- FALSE
    d$design_h[d_h] <- TRUE
    
    # Multiple runs per condition:
    if (includeMultiRuns) {
      d_hm <- which(d$designConfigID %in% designSpace$designConfigID[s])
      d$design_hm <- FALSE
      d$design_hm[d_hm] <- TRUE
    }
  }
  
  
  
  # i) I-OD quadratic __________________________________________________________
  # An optimal D design, this time specified on a quadratic model.
  if("design_i" %in% designs){
    s <- AlgDesign::optFederov(
      frml = quadraticFormula, #################
      center = TRUE,
      nTrials = sizeExperiment,
      data = designSpace[,params],
      criterion = "I" ####################
    )$rows
    
    # One run per condition
    d_i <- numeric()
    for (i in 1:length(s)) {
      IDs <- which(d$designConfigID == designSpace$designConfigID[s[i]])
      ifelse(
        length(IDs) == 1,
        d_i[i] <- IDs,
        d_i[i] <- sample(IDs, size = 1)
      )
    }
    d$design_i <- FALSE
    d$design_i[d_i] <- TRUE
    
    # Multiple runs per condition:
    if (includeMultiRuns) {
      d_im <- which(d$designConfigID %in% designSpace$designConfigID[s])
      d$design_im <- FALSE
      d$design_im[d_im] <- TRUE
    }
  }
  
  
  
  return(d)
}





# Miscellanea___________________________________________________________________
#

# Transforms vecor x to range in [0,1]
normalize <- function(x) {
  min <- min(x, na.rm = TRUE)
  max <- max(x, na.rm = TRUE)
  ifelse(
    min - max == 0,
    x[!is.na(x)] <- 1, # if no variability in x, all values are assumed to be 1.
    x <- (x - min) / (max - min)
  )
  return(x)
}

standardize <- function(x) (x - mean(x)) / sd(x)


# Clips vector x by replacing all values above 1 with 1, and all values below 0
# with 0
truncate <- function(x, min = 0, max = 1) {
  sapply(
    x, \(x) ifelse(x > max, return(1), ifelse(x < min, return(0), return(x)))
  )
}


# A function that tries to re-evaluate an expression a bunch of times if it
# fails with an error:
redoIfError <- function(attempts, expression) {
  count <- 1
  s <- NULL
  while (count <= attempts & is.null(s)) {
    s <- tryCatch(
      expr = expression,
      error = \(msg) {
        count <<- count + 1
        return(NULL)
      }
    )
  }
  return(s)
}


# A function that helps adding a suffix to the axis text in a ggplot. It can be
# used e.g. to add units to the plot.
addSuffix <- \(x, suffix = "%", ...) format(paste0(x, suffix), ...)
