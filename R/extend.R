# Make intervals for IMPUTE2
# The logic is
# 1. Find min max SNP positions from gfeatures per chromosome
# 2. Create intervals (see recoup function)
# 3. Export gwe per chromosome
# 4. Convert to PED... PLINK needed...
# 5. Convert to IMPUTE gen... GTOOL needed...
# 6. Run per chromosome
# We are going to need some kind of automation for reference panels
