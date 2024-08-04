#  Utility functions



# For drawing optimal designs from a candidate set. It returns row indices of
# the candidate set.
drawOD <- \(
  candidateSet, # a data.frame with all conditions to choose from
  nTrials, # experiment size
  formula,
  augment = FALSE,
  augmentRows = NA, # if augmented, specify the vector of row indexes.
  replace = TRUE, # FALSE uses "algDesign"; TRUE uses "skpr"
  optimality = "D" # optimality criterion. "D", "A", "I" are generally supported
  #designSpace
) {
  
  if (replace) {
    ifelse(
      augment,
      OD <- skpr::gen_design(
        candidateset = candidateSet, 
        model = formula,#~ (.) ^ 2,
        trials = nTrials,
        augmentdesign = candidateSet[augmentRows,],
        optimality = optimality
      )[,names(candidateSet)],
      OD <- skpr::gen_design(
        candidateset = candidateSet, 
        model = formula,#~ (.) ^ 2,
        trials = nTrials,
        #augmentdesign = augmentdesign,
        optimality = optimality
      )[,names(candidateSet)]
    )
    
    x <- c() # Find the chosen row indices. Very slow implementation. ##########
    for (i in 1:nTrials) {
      for (j in 1:nrow(candidateSet)) {
        if (all(candidateSet[j,] == OD[i,])) {
          x[i] <- j
          next
        }
      }
    }
    OD <- x
    
    
  } else { # without replacement. We use algDesign::optFederov
    redoIfError(
      attempts = 5,
      expression = {
        #urs <- sample( # 
        #  1:nrow(designSpace),
        #  size = floor(sizeExperiment / 2),
        #  replace = FALSE
        #)
        ifelse(
          augment,
          OD <- AlgDesign::optFederov( ## If augment __________________________
            frml = formula, # Linear effects, with all two-way interactions
            center = TRUE,
            nTrials = nTrials,
            augment = TRUE,
            rows = augmentRows, # this is where the uniform random sample is fed in
            data = candidateSet,#designSpace[,params],#####################
            criterion = "D" 
          )$rows,
          OD <- AlgDesign::optFederov( ## If no augment _______________________
            frml = formula, # Linear effects, with all two-way interactions
            center = TRUE,
            nTrials = nTrials,
            augment = FALSE, ####
            #rows = rows, # this is where the uniform random sample is fed in
            data = candidateSet,#designSpace[,params],#####################
            criterion = "D"
          )$rows
        )
        
        ## And then we draw an optimal D design, including the uniform sample:
        #OD <- AlgDesign::optFederov(
        #  frml = formula, # Linear effects, with all two-way interactions
        #  center = TRUE,
        #  nTrials = nTrials,
        #  augment = TRUE,
        #  rows = rows, # this is where the uniform random sample is fed in
        #  data = designSpace[,params],
        #  criterion = "D" #####################
        #)$rows
      }
    )
  }
  
  return(OD)
}


# Extract samples from a dataset according to the various sampling approaches:

drawSamples <- function(
    d, 
    designSpace,
    paramSpace = designSpace,
    designs = paste0("design_", letters[1:9]),# designs A-I where A is factorial
    includeMultiRuns = FALSE, # FALSE or an integer number.
    replacement = TRUE,# for optimal designs, FALSE uses algDesign; TRUE "skpr"
    sizeExperiment = 15 # proportion of a factorial design.
) {
  
  ##############################################################################
  # Step 1: Selecting parameter configurations
  
  # a) full factorial __________________________________________________________
  # Starting with one run per condition:
  if ("design_a" %in% designs) {
    d_a <- d_am <- numeric()
    
    for (i in 1:nrow(paramSpace)) {
      #IDs <- which(d$designConfigID == paramSpace$designConfigID[i])
      IDs <- which(d$paramConfigID == paramSpace$designConfigID[i])
      IDs <- IDs[!IDs %in% d_a] 
      ifelse(
        length(IDs) == 1,
        d_a[i] <- IDs,
        d_a[i] <- sample(IDs, size = 1)
      )
      # Multiple runs per condition
      if (includeMultiRuns != FALSE) {
        if (length(IDs) > 1) { # If there are many runs to choose from...
          d_am <- c(
            d_am,
            sample( # ... choose as many as possible ...
              IDs,
              size = min(c(length(IDs),includeMultiRuns)),
              replace = FALSE
            )
          )
        } else d_am <- c(d_am, IDs) # ... but if there's only one pick only that
      }
      
    }
    #for (i in 1:nrow(paramSpace)) d_a[i] <- sample(
    #  which(d$paramConfigID == paramSpace$paramConfigID[i]),
    #  size = 1
    #)
    d$design_a <- FALSE
    d$design_a[d_a] <- TRUE
    
    if (includeMultiRuns != FALSE) {
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
    d_b <- d_bm <- numeric()
    for (i in 1:length(s)) {
      IDs <- which(d$designConfigID == designSpace$designConfigID[s[i]])
      ifelse(
        length(IDs) == 1,
        d_b[i] <- IDs,
        d_b[i] <- sample(IDs, size = 1)
      )
      # Multiple runs per condition
      if (includeMultiRuns != FALSE) {
        if (length(IDs) > 1) { # If there are many runs to choose from...
          d_bm <- c(
            d_bm,
            sample( # ... choose as many as possible ...
              IDs,
              size = min(c(length(IDs),includeMultiRuns)),
              replace = FALSE
            )
          )
        } else d_bm <- c(d_bm, IDs) # ... but if there's only one pick only that
      }
      
    }
    d$design_b <- FALSE
    d$design_b[d_b] <- TRUE
    
    if (includeMultiRuns != FALSE) {
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
    ifelse(
      sizeExperiment == nrow(designSpace),
      s <- 1:sizeExperiment, # return a factorial
      s <- clhs::clhs(x = designSpace, size = sizeExperiment) # return a LHS
    )
    
    
    # One run per condition
    d_c <- d_cm <- numeric()
    for (i in 1:length(s)) {
      IDs <- which(d$designConfigID == designSpace$designConfigID[s[i]])
      ifelse(
        length(IDs) == 1,
        d_c[i] <- IDs,
        d_c[i] <- sample(IDs, size = 1)
      )
      # Multiple runs per condition
      if (includeMultiRuns != FALSE) {
        if (length(IDs) > 1) { # If there are many runs to choose from...
          d_cm <- c(
            d_cm,
            sample( # ... choose as many as possible ...
              IDs,
              size = min(c(length(IDs),includeMultiRuns)),
              replace = FALSE
            )
          )
        } else d_cm <- c(d_cm, IDs) # ... but if there's only one pick only that
      }
    }
    d$design_c <- FALSE
    d$design_c[d_c] <- TRUE
    
    if (includeMultiRuns != FALSE) {
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
  
  # ... as well as a linear one, with all interactions:
  formula <- as.formula("~ (.) ^ 2")
  
  # Now everything's ready for generating optimal designs:
  
  
  # d) D-OD linear + uniform sampling __________________________________________
  if ("design_d" %in% designs) {
    # Here we split the sample into two subsamples: one from an opimal D design
    # and one from a uniform sample, here called "urs".
    # The function AlgDesign::optFederov handles this out of the box.
    # 
    urs <- sample(
      1:nrow(designSpace),
      size = floor(sizeExperiment / 2),
      replace = FALSE
    )
    s <- drawOD(
      candidateSet = designSpace[,params],
      nTrials = sizeExperiment,
      formula = formula,
      augment = TRUE,
      augmentRows = urs,
      replace = replacement,
      optimality = "D"
    )
    
    # Now, our resulting samples for (e) D-OD + uniform sampling are:
    #
    # One run per condition:
    d_d <- d_dm <- numeric()
    for (i in 1:length(s)) {
      IDs <- which(d$designConfigID == designSpace$designConfigID[s[i]])
      #print(IDs)
      IDs <- IDs[!IDs %in% d_d]# removing runs already chosen
      ifelse(
        length(IDs) == 1,
        d_d[i] <- IDs,
        d_d[i] <- sample(IDs, size = 1)
      )
      # Multiple runs per condition
      if (includeMultiRuns != FALSE) {
        if (length(IDs) > 1) { # If there are many runs to choose from...
          d_dm <- c(
            d_dm,
            sample( # ... choose as many as possible ...
              IDs,
              size = min(c(length(IDs),includeMultiRuns)),
              replace = FALSE
            )
          )
        } else d_dm <- c(d_dm, IDs) # ... but if there's only one pick only that
      }
    }
    d$design_d <- FALSE
    d$design_d[d_d] <- TRUE
    
    if (includeMultiRuns != FALSE) {
      d$design_dm <- FALSE
      d$design_dm[d_dm] <- TRUE
    }
  }
  
  
  
  # e) D-OD quadratic __________________________________________________________
  # An optimal D design, this time specified on a quadratic model.
  if("design_e" %in% designs){
    
    s <- drawOD(
      candidateSet = designSpace[,params],
      nTrials = sizeExperiment,
      formula = quadraticFormula,
      augment = FALSE,
      replace = replacement,
      optimality = "D"
    )
    #s <- AlgDesign::optFederov(
    #  frml = quadraticFormula,
    #  #frml = ~ (.) ^ 2,
    #  center = TRUE,
    #  nTrials = sizeExperiment,
    #  data = designSpace[,params],
    #  criterion = "D"
    #)$rows
    
    # One run per condition
    d_e <- d_em <- numeric()
    for (i in 1:length(s)) {
      IDs <- which(d$designConfigID == designSpace$designConfigID[s[i]])
      IDs <- IDs[!IDs %in% d_e]# removing runs already chosen
      ifelse(
        length(IDs) == 1,
        d_e[i] <- IDs,
        d_e[i] <- sample(IDs, size = 1)
      )
      # Multiple runs per condition
      if (includeMultiRuns != FALSE) {
        if (length(IDs) > 1) { # If there are many runs to choose from...
          d_em <- c(
            d_em,
            sample( # ... choose as many as possible ...
              IDs,
              size = min(c(length(IDs),includeMultiRuns)),
              replace = FALSE
            )
          )
        } else d_em <- c(d_em, IDs) # ... but if there's only one pick only that
      }
    }
    d$design_e <- FALSE
    d$design_e[d_e] <- TRUE
    
    if (includeMultiRuns != FALSE) {
      d$design_em <- FALSE
      d$design_em[d_em] <- TRUE
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

# translate p-values into significance stars
ptostars <- \( 
  x, # p-value(s)
  labels = c("***", "**", "*", ""),
  levels = c(0, 0.001, 0.01, 0.05, 1)
) {
  labels[findInterval(x, levels)]
}

# A function that tries to re-evaluate an expression a bunch of times if it
# fails with an error:
redoIfError <- \(attempts, expression) {
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
