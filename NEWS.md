# Development version
* Compatibility with ggplot 3.4.0 (#5 @cregouby)

# plotluck 1.1.1
* Rerun roxygen2 for documentation format compliance
* Move example image in readme file to avoid warning
* No functional change

# plotluck 1.1.0
* Compatibility with ggplot 2.2.0
* Swap order of data and formula arguments, as requested

# plotluck 1.0.0 (2016-09-22)
* First CRAN release
* Compatibility with ggplot 2.1.0
* Changed _plotluck()_ arguments to use a formula instead of arguments _x,y,z_
* Removed _plotluck.multi()_ from being public - functionality can now be accessed using dot symbol _'.'_ in plotluck formula
* Added function sample.plotluck() to generate random plots from a data set
* Changed heat map visualization to vary rectangle size with count/weight
* Renamed and removed some options
* Added new option _verbose_ that outputs information about plot types, log scaling, etc
* Added new override option _geom_ to directly set plot type
* Eliminate _exclude.factor_ arguments; instead, use _na.rm_ to remove missing values at the start if requested, else keep throughout
* Estimate size of legend to decide wheter vertical or horizontal placement takes less space
* Estimate axis label overlap to write at an angle
* Select palette or colorbar based on variable type; use brewer palettes
* Use entropy order for lattice (2-variable) layout as well
* Improve factor ordering algorithm
* Several bug fixes

# plotluck 0.1.1
* Bug fix for scatter plots

# plotluck 0.1.0
* Initial Release
