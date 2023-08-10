# ==============================================================================
# Formula and mass utility functions
# ==============================================================================
# load required libraries ------------------------------------------------------
library(MetaboCoreUtils)

# Retention time indexing for LC-MS --------------------------------------------
# data frame with RI markers
naps <- data.frame(rtime = c(0.456, 0.497, 0.595, 0.874, 2.500, 3.736,
                             4.414,	4.957, 5.432,	5.875),
                   rindex = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000))

# data frame with example compounds
metabolites <- data.frame(name = c("O-acetylcarnitine", "2-methylbutyrylcarnitine",
                                   "tauroursodeoxycholic acid"),
                          rtime = c(0.46, 2.759, 5.222))

# perform indexing using linear interpolation
rindex <- indexRtime(metabolites$rtime, naps)
rindex

# perform secondary correction
ref <- data.frame(rindex = c(109.7561, 855.7895),
                  refindex = c(110, 850))

correctRindex(rindex, ref)
