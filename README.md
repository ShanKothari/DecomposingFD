# DecomposingFD
Code to calculate functional trait dispersion as defined by Scheiner et al. (2016), and its decomposition into richness, evenness, and dispersion.

## What is functional trait dispersion?

Scheiner et al. (2016) defined a new Hill number-based metric of functional diversity called functional trait dispersion that can be decomposed into three components: Species richness, functional evenness, and (mean) functional dispersion. The result is a metric of functional diversity in terms of effective number of functionally distinct species. It can be considered at different scales (alpha, beta, gamma) though it is not multiplicative. Because it can be decomposed into interpretable components, one can see exactly what aspect(s) of functional diversity are changing, e.g. across a gradient.

## Status of package
This is eventually to become a package, though currently it is just a couple poorly annotated scripts. Stay tuned for expansions! Example data (for examples soon to come) come from Etienne Laliberte's FD package, released by Dr. Laliberte under a GPL-2 license.

## How to use
Currently, the best way to use these scripts are to download them (or clone the repository) and source the functions into R.

AlphaFD.R contains two functions -- FTD, which calculates functional trait dispersion for a single community based on its species and their traits, and FTD.comm, which wraps FTD across a set of communities. Each function also calculates how FTD is decomposed into species richness, functional evenness, and functional dispersion. FTD is measured as the effective number of species that are as functionally distinct as the maximally distinct species in the data set. It can also be weighted by abundance. Importantly, to compare FTD values for two communities, they must be treated as arising from the same species pool, since the identity of the maximally distinct pair of species varies based on the species pool.

BetaFD.R contains two functions -- comm.disp, which calculates the mean distance between two communities, and FTD.beta, which uses comm.disp to calculate the beta FTD of a set of communities, as well as its decomposition into number of communities, evenness of community dispersion, and mean dispersion among communities. FTD.beta is measured as the effective number of communities that are as functionally distinct as the maximally distinct communities possible.

GammaFD.R contains code in progress to calculate unstructured functional gamma diversity.

## Literature reference

Scheiner, S. M., E. Koman, S. J. Presley, and M. R. Willig (2016). Decomposing functional diversity. Methods in Ecology and Evolution (early view). DOI: 10.1111/2041-210X.12696
