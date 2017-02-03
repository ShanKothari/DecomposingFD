# DecomposingFD
Code to calculate functional trait diversity as defined by Scheiner et al. (2016), and its decomposition into richness, evenness, and dispersion.

## Status of package
This is eventually to become a package, though currently it is just a couple poorly annotated scripts. Stay tuned for expansions!

## How to use
Currently, the best way to use this script is to download the script (or close the repository) and source the functions into R. AlphaFD.R contains two functions -- FTD, which calculates functional trait dispersion for a single community based on its species and their traits, and FTD.comm, which wraps FTD across a set of communities. Each function also calculates how FTD is decomposed into species richness, functional evenness, and functional dispersion.

## Literature reference

Scheiner, S. M., E. Koman, S. J. Presley, and M. R. Willig (2016). Decomposing functional diversity. Methods in Ecology and Evolution (early view). DOI: 10.1111/2041-210X.12696
