######################################
# 04.23.2017
# Emsemble plots example
# BISC 577
######################################

# Initialization
library(DNAshapeR)

# Extract sample sequences
fn <- "/Users/lester/BISC577/CTCF_bound_500.fa"

# Predict DNA shapes
pred <- getShape(fn)

# Generate ensemble plots
plotShape(pred$MGW)
heatShape(pred$ProT, 20)
