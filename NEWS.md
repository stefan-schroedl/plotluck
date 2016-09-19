# plotluck 1.0.0

* Compatibility with ggplot 2.1.0
* Changed plotluck() arguments to use a formula instead of arguments x,y,z
* Removed plotluck.multi from being public - functionality can now be accessed by using '.' in plotluck formula
* Changed heat map visualization to change size with count/weight
* Added function 'sample.plotluck' to generate random plots from a data set
* Renamed and cleaned up some option names
* Added new option 'verbose' that outputs information about plot types, log scaling, etc
* Added new override option 'geom' to directly set plot type
* Eliminate exclude.factor arguments; instead, remove NA values at the start if requested, else keep throughout
* Estimate size of legend to decide vertical/horizontal placement
* Estimate axis label overlap to write at an angle
* Select palette or colorbar based on variable type
* Use entropy order for lattice (2-variable) layout as well
* Improve factor ordering algorithm
* Several bug fixes

# plotluck 0.1.1

* Bug fix for scatter plots

# plotluck 0.1.0

* Inital Release