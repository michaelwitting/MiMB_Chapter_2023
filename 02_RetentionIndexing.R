# ==============================================================================
# Retention time indexing and effective mobility transformation
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

# Effective mobility transformatin for CE-MS -----------------------------------
# data frame with EOF marker
eof <- data.frame(rtime = c(840.796),
                  mobility = c(0))

# data frame with example compounds
metabolites <- data.frame(name = c("procaine"),
                          rtime = c(450.239))

# perform effective mobility transformation
eff <- convertMtime(metabolites$rtime / 60,
                    rtime = eof$rtime / 60,
                    mobility = eof$mobility,
                    tR = 3 / 60,
                    U = 30,
                    L = 800)
eff