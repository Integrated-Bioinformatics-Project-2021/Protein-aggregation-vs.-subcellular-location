# The structural analysis

# ------------------------------------------------------------------------------
# Steps:
# 1) Download all human AlphaFold structures.
# 2) For each domain in a protein, generate a "sub-structure" based on the Alphafold structures.
# 3) Calculate the contact order of each domain (you can use the perl script from David Bakers lab).
# The R object "domains.RData" contains all domains in human proteins, determined by Gene3D software.
# ------------------------------------------------------------------------------

# 1) Download all human AlphaFold structures -----------------------------------
# Go to https://www.alphafold.ebi.ac.uk/download and download the entire human dataset.
# Save it to the /Date folder and untar it.

# 2) Generate a "sub-structure" ------------------------------------------------
# Check if all proteins in data are present in domains:
List1 <- data$Protein[data$Subcellular_location != ""]
length(unique(List1))
List2 <- unique(domains$Protein)
length(List2)
length(Reduce(intersect,list(List1,List2)))

# TODO: calculate the structure of a single domain --> pdb file

# 3) Calculate the contact order of each domain --------------------------------
# https://depts.washington.edu/bakerpg/contact_order/ (Download Contact Order Program (Perl version))

# contactOrder.pl <options> <pdbfile>
# with <options> = -r: returns relative CO
# with pdbfile the pdb file with the AlphaFold structure of a single domain
