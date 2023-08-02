# ==============================================================================
# Formula and mass utility functions
# ==============================================================================
# load required libraries ------------------------------------------------------
library(MetaboCoreUtils)
library(Rdisop)
library(MetaboAnnotation)

# Parsing of chemical formulas -------------------------------------------------
# create vector with chemical formulas
chem_formula <- c("C6H12O6", "C6H13NO2")

# count elements
chem_formula_count <- countElements(chem_formula)

# paste elements
chem_formula_pasted <- pasteElements(chem_formula_count)

# create isotopic formulae
iso_formula <- c("[13C6]H12O6", "[13C3]C3H12O6")
iso_formula_count <- countElements(iso_formula)
iso_formula_pasted <- pasteElements(iso_formula_count)

# standardize formula
chem_formula_nonstd <- c("H12O6C6")
chem_formula_std <- standardizeFormula(chem_formula_nonstd)

# Calculation of exact masses, Kendrick masses or Kendrick mass defect ---------
# calculate masses
calculateMass(chem_formula)
calculateMass(chem_formula_pasted)
calculateMass(iso_formula)
calculateMass(iso_formula_pasted)

# Kendrick mass defect

# Working with adducts ---------------------------------------------------------
# calculate adduct m/z
exact_mass <- calculateMass(chem_formula)
mass2mz(exact_mass, c("[M+H]+", "[M+Na]+"))

adductNames()
adductNames(polarity = "positive")
adductNames(polarity = "negative")

mz2mass(c(181.0707, 203.0526,
          132.1019, 154.0838),
        c("[M+H]+", "[M+Na]+"))

# Working with isotopic pattern ------------------------------------------------
# calculate mass (Rdisop)
getMass(getMolecule("C6H12O6"))
getMass(getMolecule("C6H13O6", z = 1))

# compare mass calculated by Rdisop and MetaboCoreUtils
getMass(getMolecule("C6H12O6")) - calculateMass("C6H12O6")

# calculate isotopic pattern
getIsotope(getMolecule("C6H13O6", z = 1))
getIsotope(getMolecule("C6H12O6Na", z = 1))

# isotopic pattern of Tryptophan
iso_pattern <- data.frame(masses = c(205.09715, 206.100189, 207.102574),
                          intensities = c(100.000, 12.854, 1.170))

plot(iso_pattern$masses, iso_pattern$intensities, type = "h")

# calculate formulas from mass
decomposeMass(iso_pattern$masses[1], z = 1, ppm = 5.0)

# calculate formulas from isotope pattern
decomposeIsotopes(iso_pattern$masses,
                  iso_pattern$intensities,
                  z = 1,
                  ppm = 5.0)

# matching of chemical formulas ------------------------------------------------
# create vectors with chemical formulas
queries <- c("C6H12O6", "C11H12N2O2")
targets <- c("C6H12O6", "C6H13NO2")

# match chemical formulas
matchFormula(queries, targets)
matchedData(matchFormula(queries, targets))

# create vectors with non-standard formulas
queries <- c("H12O6C6", "C11H12N2O2")
targets <- c("C6H12O6", "C6H13NO2")

# match chemical formulas with non-standard formulas
matchFormula(queries, targets)
matchedData(matchFormula(queries, targets))
