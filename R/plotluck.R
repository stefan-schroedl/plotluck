#'@import ggplot2
#'@importFrom plyr ddply summarise .
#'@importFrom stats as.formula median quantile runif xtabs

safe_log <- function(x,...) {
   ifelse(x<=0,0,log(x,...))
}

group.center <- function(data, x, group.vars, w='NULL', method='median') {

   if (method == 'median') {
      fun <- median
   } else if (method == 'mean') {
      fun <- mean
   } else {
      stop(sprintf('group.center unknown function: %s', fun))
   }

   if (w == 'NULL') {
      grp.center <- ddply(data, group.vars, function(d) fun(d[[x]], na.rm=TRUE))
   } else {
      if (method == 'median') {
         # workaround for problems in wtd.quantile:
         # - if all values are NA, wtd.mean throws error: 'zero non-NA points'
         # - if weights are integer, sum can lead to NA due to overflow (conversion to double in plotluck)
         f <- function(d) if (all(is.na(d[[x]]))) { NA } else { Hmisc::wtd.quantile(d[[x]], weights=d[[w]], probs=0.5, na.rm=TRUE, normwt=TRUE) }
         grp.center <- ddply(data, group.vars, f)
      } else {
         grp.center <- ddply(data, group.vars, function(d) Hmisc::wtd.mean(d[[x]], d[[w]], na.rm=TRUE))
      }
   }
   names(grp.center)[length(grp.center)] <- '.center.'
   return(grp.center)
}

# calculate axis scale transform based on value distribution
# either:
# - identity
# - log-10
# - log-10-modulus: sign(x) * log(1+|x|)
# John & Draper 1980; see http://blogs.sas.com/content/iml/2014/07/14/log-transformation-of-pos-neg/

get.trans.fun <- function(x, expansion.thresh=2) {

   if (!is.numeric(x) ||
          length(unique(x)) <= 2) {
      return(scales::identity_trans())
   }

   q <- quantile(x, c(0, 0.25, 0.5, 0.75, 1), na.rm=TRUE)

   # heuristic:
   # - compute the ratio of the plotting range occupied by the
   #   'core' part of the distribution between the lower and upper
   #   quartile
   # - if range is positive only, candidate transform is log
   # - if range includes negative values, candidate transform is log-modulus
   # - compute this ratio for the case that transform is applied
   # - decide for transform is this ratio is greater than expansion.thresh

   range.core     <- c(q[4], q[2])
   range.complete <- c(q[5], q[1])
   ratio.before   <- diff(range.core) / diff(range.complete)

   sgn <- ifelse(q[5] <= 0,
                 -1,
                 ifelse(q[1] > 0,
                        1,
                        0))

   # offset == 1 is interpreted as using log-modulus, instead of log
   if (sgn == 1) {
      offset <- 0
   } else {
      offset <- 1
   }

   fun.trans <- log_mod_trans(base=10, offset=offset)

   range.core.after     <- fun.trans$trans(range.core)
   range.complete.after <- fun.trans$trans(range.complete)
   ratio.after          <- diff(range.core.after) / diff(range.complete.after)
   if (ratio.after > expansion.thresh * ratio.before) {
      return(fun.trans)
   } else {
      return(scales::identity_trans())
   }
}


# calculate the floor or ceiling of the exponent of a number
# note: for log-modulus (offset=1), sign of returned value signifies sign of
#   input value x; therefore, different exponents cannot be represented for |x| < 1.
# note: for log-modulus, in order to distinguish between -10^0, 0, and 10^0,
#   we encode -10^x=-x, -10=-1, -1=0, 0=1, 1=2, 10^x=x+2, ..

encode.int.power <- function(x, base=10, floor=TRUE, offset=0) {

   # if floor=FALSE, compute ceiling
   # -1 * floor(-1*x) = ceiling(x)
   floor.factor <- ifelse(floor, 1, -1)

   if (offset == 0) {
      # regular log transform
      y <- floor.factor * floor(floor.factor * safe_log(x, base))
   } else {
      # log-modulus transform
      # no fractional values!
      x <- floor.factor * floor(floor.factor * x)

      y <- floor.factor * floor(floor.factor * sign(x) * safe_log(abs(x), base))

      # encoding as described in note above
      y <- y + ifelse(x >= 1, 2, ifelse(x >= 0, 1, 0))
   }
   return(y)
}


decode.int.power <- function(x, base=10, offset=0) {
   if (offset == 0) {
      # regular log transform
      return(base^x)
   } else {
      # log-modulus transform
      return (ifelse(x <= 0, - base^abs(x),
                     ifelse(x == 1, 0,
                            base^abs(x-2))))
   }
}


# calculate breaks for log/log-modulus scale
log_mod_breaks <- function(n=5, base=10, offset=0)
{
   function(x) {

      rng <- range(x, na.rm=TRUE)

      if (any(is.infinite(rng))) {
         # Note: known issue with log scaling
         # https://github.com/hadley/ggplot2/issues/930
         # seems to be happening mostly with scatter plot, density/hex less affected
         warning('Skipping tick mark computation due to existing ggplot issue.
                 Try using hex or density instead of scatter plot or histogram.')
         return(1)
      }

      min.pos <- encode.int.power(rng[1], base=base, floor=TRUE,  offset=offset)
      max.pos <- encode.int.power(rng[2], base=base, floor=FALSE, offset=offset)

      # the range of powers that the values span
      total  <- abs(max.pos - min.pos)

      if (total + 1 <= n/2) {
         # less than half of requested ticks:
         # introduce subdivisions by lowering base
         sub <- floor(n/(total + 1))
         base <- base^(1/sub)
         min.pos <- encode.int.power(rng[1], base=base, floor=TRUE,  offset=offset)
         max.pos <- encode.int.power(rng[2], base=base, floor=FALSE, offset=offset)
         by <- 1
      } else {
         # enough or more dynamic range than breaks requested
         by <- max(1, floor(total/(n-1)))
      }

      breaks <- sapply(seq(min.pos, max.pos, by),
                       function(x) decode.int.power(x,base=base, offset=offset))

      return(breaks)
   }
}


# log (offset=0) or log-modulus (offset=1) scale
# John & Draper 1980; see e.g. http://blogs.sas.com/content/iml/2014/07/14/log-transformation-of-pos-neg/
log_mod_trans <- function(base = exp(1), offset=0)
{
   eps<-1e-100
   if (offset == 0) {
      # regular log transform
         trans <- function(x) log(x, base)
         inv   <- function(x) base^x
         domain = c(eps, Inf)
   } else {
      # log-modulus transform
      trans <- function(x) sign(x) * log(abs(x) + offset, base)
      inv   <- function(x) sign(x) * (base^abs(x) - offset)
      domain = c(-Inf, Inf)
   }
   scales::trans_new(paste0("log-mod-", format(base)),
             trans,
             inv,
             log_mod_breaks(base=base, offset=offset),
             domain=domain)
}


format.label.trans <- function(x, trans.fun) {
   if (trans.fun$name == 'identity') {
      return(x)
   }
   return(sprintf('%s (log scale)', x))
}

# change the levels of a factor such as to indicate the column name
# for the first level
format.facets <- function(data, x) {
   if (is.numeric(data[[x]])) {
      data[[x]]
   } else {
      levels <- levels(data[[x]])
      levels[1] <- sprintf('%s = %s', x, levels[1])
      levels(data[[x]]) <- levels
      return(data[[x]])
   }
}

# return levels of (unordered) factor, sorted by (weighted) frequency
order.factor.freq <- function(data, x, w='NULL', decreasing=TRUE,
                              if.ordered=FALSE, exclude.factor=NA) {
   x.data <- data[[x]]
   if ((!(is.factor(x.data) ||
             is.character(x.data))) ||
          # do not reorder ordered factors, unless explicitly overriden with 'if.ordered'
          ((!if.ordered) && is.ordered(x.data))) {
      # no change
      return(levels(x.data))
   }

   if (w=='NULL') {
      tab <- table(x.data, exclude=exclude.factor)
   } else {
      tab <- xtabs(as.formula(sprintf('%s ~ %s', w, x)), data, exclude=exclude.factor)
   }
   ord <- order(tab, decreasing=decreasing)
   return(levels(x.data)[ord])
}


# return levels of (unordered) factor x, sorted by (weighted) mean of dependent variable y
order.factor.value <- function(data, x, y, w='NULL', decreasing=TRUE, exclude.factor=NA) {

   if (y == 'NULL' || y == x) {
      return(order.factor.freq(data,x, w, decreasing, exclude.factor=exclude.factor))
   }

   x.data <- data[[x]]
   y.data <- data[[y]]
   if ((!(is.factor(x.data) || is.character(x.data)))
       || is.ordered(x.data)) {
      return(levels(x.data))
   }

   if (!(is.numeric(y.data) || is.logical(y.data))) {
      if (!is.ordered(y.data) || is.logical(y.data)) {
         return(levels(x.data))
      } else {
         # for an ordered factor, consider integer levels
         data[[y]] <- as.numeric(data[[y]])
      }
   }

   means <- group.center(data, y, x, w)

   # if medians are equal for two or more classes, order by mean
   if (length(unique(means$.center.)) != nrow(means)) {
      means <- group.center(data, y, x, w, method='mean')
   }
   ord <- order(means$.center., decreasing=decreasing)
   return(means[ord, x])
}


# return a possibly refactored vector:
# - convert character vectors into factors
# - if there are more than max.factor.levels levels, retain only the most frequent ones,
#   merge the remaining ones into '.other.'
# - convert a numeric variable with less than few.unique.as.factor distinct values
#   into a factor, if length would have allowed for more distinct values
# - exclude-factor - flag to apply in factor(...) calls, to preserve NA values
#   if required
preprocess.factors <- function(data, x, w='NULL',
                              max.factor.levels=100,
                              few.unique.as.factor=10,
                              exclude.factor=NA) {
   x.data <- data[[x]]
   non.na <- !is.na(x.data)
   u <- length(unique(x.data[non.na]))
   if (u == 1) {
      # note: case of all NA already caught in function plotluck
      stop(sprintf('Variable %s has only single level ("%s"), giving up', x, unique(x.data[non.na])))
   }

   if (is.numeric(x.data)) {
      if (u <= few.unique.as.factor &&
             nrow(data) > u * u) { # heuristic similar to histogram binning
         return(ordered(x.data, exclude=exclude.factor))
      } else {
         return(x.data)
      }
   }

   if (!is.factor(x.data)) {
      # i.e., character or logical
      x.data <- factor(x.data, exclude=exclude.factor)
   }

   if (u > max.factor.levels) {
      warning(sprintf('Factor variable %s has too many levels (%d), truncating to %d', x, u, max.factor.levels))
      tab <- table(x.data)
      data[[x]] <- x.data
      levels.ord <- order.factor.freq(data, x, w, if.ordered=TRUE, exclude.factor=exclude.factor)
      levels.reduced <- levels.ord[1:max.factor.levels]
      idx.trunc <- !(x.data %in% levels.reduced)
      levels(x.data) <- c(levels(x.data), '.other.')
      x.data[idx.trunc] <- '.other.'
      x.data <- factor(x.data, levels=c(levels.reduced, '.other.'), exclude=exclude.factor)
   } else if (u == 2) {
      # binary variables are trivially ordered
      x.data <- ordered(x.data)
   }

   # suppress unused levels
   # for NA treatment, see http://r.789695.n4.nabble.com/droplevels-inappropriate-change-td4723942.html
   # droplevels(x.data, exclude=exclude.factor) does not seem to be working correctly
   return(x.data[,drop=T])
}


# determine number of histogram bins
# note: this is a modification the nclass.FD function, which sometimes returns '1'
# for very uneven distributions.

nclass.FD.modified <- function(x, max.bins=200) {
   u <- unique(x)
   if (length(u) == 1) {
      # only one value
      return(1L)
   }
   # minimum reasonable bin size
   h <- 0.5 * min(abs(diff(u)), na.rm=TRUE)

   # for very few unique values, show all of them individually
   if (length(u) > 10) {
      h2 <- 2 * stats::IQR(x, na.rm=TRUE)
      if (h2 == 0)
         h2 <- 2 * stats::mad(x, constant=2, na.rm=TRUE)
      if (h2 > 0)
         h <- max(h, h2 * length(x)^(-1/3))
   }

   # capping the number of bins, since for some features the estimated bins are very high causing the hist() function to throw an error
   num.bins <- ceiling(diff(range(u, na.rm=TRUE))/h)
   num.bins <- min(max.bins,
                   max(min(5, length(u), max.bins, length(x)/log2(length(x))),
                       num.bins))
   num.bins
}

# break a variable into intervals
# - method=quantile - bins with equal number of samples
# - method=histogram - equidistant intervals
# - max.breaks - requested number of breaks
# - keep.breaks - if there are less than that many distinct values to
#     begin with, don't do additional binning
#
# for factor variables with more than keep.breaks many levels:
# - if y is not 'NULL', order levels by mean y, else by frequency
#   (y is not used for numeric variables)
# - bin the corresponding integer vector of level IDs
# - the transformed level names are lists of levels
quantize <- function(data, x,
                     y='NULL',
                     w='NULL',
                     max.breaks=5,
                     keep.breaks=10,
                     method=c('quantile', 'histogram'),
                     exclude.factor=NA) {
   x.data <- data[[x]]
   u <- unique(x.data)
   if (length(u) <= keep.breaks) {
      return(ordered(x.data))
   }

   method <- match.arg(method)

   max.breaks <- min(max.breaks, nrow(data))

   if (!is.numeric(data[[x]])) {
      # order levels
      data[[x]] <- as.factor(x.data)
      x.levels <- order.factor.value(data, x, y, w, decreasing=TRUE, exclude.factor=exclude.factor)
      x.data <- factor(x.data, levels=x.levels, exclude=exclude.factor)
      # treat as integers for binning
      x.data <- as.integer(x.data)
   }

   # compute number of breaks
   if (method == 'quantile') {
      probs <- seq(0, max.breaks) / max.breaks
      if (w=='NULL') {
         breaks <- quantile(x.data, probs=probs, na.rm=TRUE)
      } else {
         breaks <- Hmisc::wtd.quantile(x.data, weights=data[[w]],
                                probs=probs, na.rm=TRUE, normwt=TRUE)
      }
   } else {
      n.breaks <- nclass.FD.modified(x.data, max.breaks)
      breaks <- graphics::hist(x.data, breaks=n.breaks, plot=FALSE)$breaks
   }

   if (!is.numeric(data[[x]])) {
      breaks <- as.integer(round(breaks))
      # dig.lab: in cut(), avoids scientific formatting for integers
      dig.lab=50
   } else {
      dig.lab = 3L # same as usual default in cut()
   }

   # break into intervals
   x.data <- ordered(cut(x.data, breaks=unique(breaks), include.lowest=TRUE, dig.lab=dig.lab))

   # recreate level information for factors
   if (!is.numeric(data[[x]])) {
      # make lists of elements, e.g.:
      #  - original levels 'a', 'b', 'c', 'd'
      #  - binned level '(2,4]' ->  'c,d'
      x.data <- as.factor(x.data)
      reconstruct.levels <- function(interval, levels) {
         left.bracket <- substr(interval,1,1)
         left.idx  <- as.integer(gsub('[(\\[]([0-9]+),.*', '\\1', interval))
         right.idx <- as.integer(gsub('.*,([0-9]+)\\]', '\\1', interval))
         if (left.bracket == '(') {
            left.idx <- left.idx + 1
         }
         # catch if level doesn't fit interval pattern
         if (is.na(left.idx) || is.na(right.idx)) {
            return(interval);
         }
         return(paste(levels[left.idx:right.idx], sep='', collapse=','))
      }

      levels(x.data) <- sapply(levels(x.data), reconstruct.levels, levels=x.levels)
      x.data <- droplevels(x.data)
   }

   return(x.data)
}


# entropy: H(vars) = - \sum p(vars)log2(p(vars))
# assumption: vars are factors
entropy <- function(data, vars, w, exclude.factor=NA) {
   tab  <- marginalize(data, vars, w, exclude.factor=exclude.factor)
   log.tab <- log2(tab)
   log.tab[tab==0] <- 0 # guard for -Inf
   return(-sum(tab * log.tab))
}


# conditional entropy: H(y|x) = H(x,y) - H(x)
# assumption: x, y are factors
cond.entropy <- function(data, x, y, w, exclude.factor=NA) {
   if (x == y) {
      return(0);
   }
   hxy <- entropy(data, c(x, y), w, exclude.factor)
   hx <- entropy(data, x, w, exclude.factor)
   return(hxy - hx)
}


marginalize <- function(data, var.list, w='NULL', exclude.factor=NA) {
   rhs <- paste(var.list, sep='', collapse='+')
   prop.table(xtabs(as.formula(sprintf('%s ~ %s', ifelse(w=='NULL', c(''), w), rhs)), data, exclude=exclude.factor))
}


# heuristic to decide whether to plot a heat map:
# - axes x, y must be either numeric, or ordered factors
# - axes must have at least min.size distinct values
# - at least a fraction min.coverage of the grid points must have data
is.grid.like <- function(data, x, y, z, min.size=3, min.coverage=0.5) {
   if (!is.numeric(data[[z]]) ||
       (!(is.numeric(data[[x]]) || is.ordered(data[[x]]))) ||
       (!(is.numeric(data[[y]]) || is.ordered(data[[y]])))) {
      return(FALSE)
   }
   ux <- length(unique(data[[x]]))
   if (ux < min.size) {
      return(FALSE)
   }
   uy <- length(unique(data[[y]]))
   if (uy < min.size) {
      return(FALSE)
   }
   u <- nrow(unique(data[,c(x, y)]))
   return(u >= min.coverage * ux * uy)
}


# 2D scatter plot of numeric or factor variables
# - overlay smoothing line if both variables are numeric
# - convert.duplicates.to.weights=TRUE - unique repeated points, count frequency
# - convert.duplicates.to.weights=FALSE - plot each repeated point separately, add
#   jitter if there are more than min.points.jitter in one location
# - weights are drawn with a shaded circle with proportional area
# - w.jitter, h.jitter - the amount of jittering, relative to total display
# - factor.jitter - amount of jitterering for factors, absolute
gplt.scatter <- function(data, x, y='NULL', w='NULL',
                         convert.duplicates.to.weights=TRUE,
                         min.points.jitter=3,
                         w.jitter=0.01,
                         h.jitter=0.01,
                         factor.jitter=0.25, ...) {

   num.x <- is.numeric(data[[x]])
   flip <- FALSE

   if (y=='NULL') {
      # for line plots, make up a constant y coordinate
      data[['.y_const']] <- 1
      y <- '.y_const'
      # treat as factor
      num.y <- FALSE
      flip=TRUE
   } else {
      num.y <- is.numeric(data[[y]])
      if (num.x && (!num.y)) {
         # HACK: ggplot2 does not implement vertical dodging, therefore we
         # use coord_flip() as a workaround
         flip <- TRUE
         tmp      <- w.jitter
         w.jitter <- h.jitter
         h.jitter <- tmp
      }
   }

   if (convert.duplicates.to.weights) {
      # record frequency/total weight for each unique row of the data
      if (w == 'NULL') {
         data <- ddply(data, names(data), function(D) nrow(D))
         names(data)[length(data)] <- '.freq.'
         w <- '.freq.'
      } else {
         data <- ddply(data, setdiff(names(data), w), function(D) sum(D[[w]]))
         names(data)[length(data)] <- w
      }
   }

   use.weights <- (w != 'NULL' && length(unique(data[[w]])) > 1)

   if (!flip) {
      p <- ggplot(data, aes_string(x=x, y=y, weight=w,
                                   ymin=min(y, na.rm=TRUE),
                                   ymax=max(y, na.rm=TRUE),
                                   xmin=min(x, na.rm=TRUE),
                                   xmax=max(x, na.rm=TRUE)))
   } else {
      p <- ggplot(data, aes_string(x=y, y=x, weight=w,
                                   ymin=min(x, na.rm=TRUE),
                                   ymax=max(x, na.rm=TRUE),
                                   xmin=min(y, na.rm=TRUE),
                                   xmax=max(y, na.rm=TRUE)))
   }

   # dodging

   if (y == '.y_const') {
      pos <- 'identity'
   } else {
      pos <- position_dodge(width=0.9)
      # note: if aesthetics xmin/xmax left unspecified, dodge width affects
      # diagram area!
   }


   alpha <- 0.6
   if (nrow(data) > 1000) {
      # if there are a lot of points, better make them transparent
      alpha = 0.1
   }

   if (use.weights) {
      p <- p + geom_point(aes_string(size=w), alpha=alpha,
                          position=pos, ...) +
            scale_size_area(guide=FALSE)
   } else {

      if (max(table(data[, c(x, y)])) >= min.points.jitter) {
         # test for repeated points, optionally jitter
         # - if both x and y are numeric, jitter in both directions
         # - if one is a factor, only jitter in this direction
         #   (enough empty space between levels)
         if (num.x) {
            if (num.y) {
               w.jitter <- w.jitter * diff(range(data[[x]], na.rm=TRUE))
               h.jitter <- h.jitter * diff(range(data[[y]], na.rm=TRUE))
            } else {
               # num.x, !num.y
               w.jitter <- 0
               h.jitter <- factor.jitter
            }
         } else {
            # !num.x
            if (num.y) {
               w.jitter <- factor.jitter
               h.jitter <- 0
            } else {
               w.jitter <- factor.jitter
               h.jitter <- factor.jitter
            }
         }
      } else {
         w.jitter <- 0
         h.jitter <- 0
      }
      if (flip) {
         tmp      <- w.jitter
         w.jitter <- h.jitter
         h.jitter <- tmp
      }
      p <- p + geom_point(alpha=alpha, position=position_jitter(width=w.jitter, height=h.jitter), ...)
   }

   if (y == '.y_const') {
      # no labels and ticks necessary
      p <- p + ylab(NULL) +
         theme(axis.ticks.y=element_blank(),
               axis.text.y=element_blank())
   } else if (num.x && num.y) {
      if (w == 'NULL' || w == '.freq.') {
         # note: for more than 1000 points, ggplot2 uses mgcv package,
         # leads to weird dependency issue
            p <- p + geom_smooth(fill='cornflowerblue', color='cornflowerblue', alpha=0.2)
         } else {
            p <- p + geom_smooth(aes(weight=w), fill='cornflowerblue', color='cornflowerblue', alpha=0.2)
         }
   }

   if (flip) {
      p <- p + coord_flip()
   }

   return(p)
}

# simpler version of gplt.scatter without flipping
# gplt.scatter <- function(data, x, y='NULL', w='NULL',
#                          convert.duplicates.to.weights=TRUE,
#                          min.points.jitter=3,
#                          w.jitter=0.01,
#                          h.jitter=0.01,
#                          factor.jitter=0.25, ...) {
#
#    num.x <- is.numeric(data[[x]])
#    pos <- 'identity'
#
#    if (y=='NULL') {
#       # for line plots, make up a constant y coordinate
#       data[['.y_const']] <- 1
#       y <- '.y_const'
#       # treat as factor
#       num.y <- FALSE
#    } else {
#       num.y <- is.numeric(data[[y]])
#    }
#
#    if (convert.duplicates.to.weights) {
#       # record frequency/total weight for each unique row of the data
#       if (w == 'NULL') {
#          data <- ddply(data, names(data), function(D) nrow(D))
#          names(data)[length(data)] <- '.freq.'
#          w <- '.freq.'
#       } else {
#          data <- ddply(data, setdiff(names(data), w), function(D) sum(D[[w]]))
#          names(data)[length(data)] <- w
#          data[[w]] <- data[[w]] / min(data[[w]], na.rm=TRUE)
#       }
#
#    } else { # convert.duplicates.to.weights
#
#       if (max(table(data[, c(x, y)])) >= min.points.jitter) {
#          # test for repeated points, optionally jitter
#          # - if both x and y are numeric, jitter in both directions
#          # - if one is a factor, only jitter in this direction
#          #   (enough empty space between levels)
#          if (num.x) {
#             if (num.y) {
#                w.jitter <- w.jitter * diff(range(data[[x]], na.rm=TRUE))
#                h.jitter <- h.jitter * diff(range(data[[y]], na.rm=TRUE))
#             } else {
#                # num.x, !num.y
#                w.jitter <- 0
#                h.jitter <- factor.jitter
#             }
#          } else {
#             # !num.x
#             if (num.y) {
#                w.jitter <- factor.jitter
#                h.jitter <- 0
#             } else {
#                w.jitter <- factor.jitter
#                h.jitter <- factor.jitter
#             }
#          }
#          pos = position_jitter(width=w.jitter, height=h.jitter)
#       }
#    }
#
#    p <-  ggplot(data, aes_string(x=x, y=y, weight=w)) + geom_point(alpha=0.6, position=pos, ...)
#
#    if (w != 'NULL' && length(unique(data[[w]])) > 1) {
#       # add a shaded halo according to weights
#       p <- p + geom_point(aes_string(size=w), position=pos, alpha=0.1, ...) +
#          scale_size_area(max_size=10, guide=FALSE)
#    }
#
#    if (y == '.y_const') {
#       # no labels and ticks necessary
#       p <- p + ylab(NULL) +
#          theme(axis.ticks.y=element_blank(),
#                axis.text.y=element_blank())
#    } else if (num.x && num.y) {
#       p <- p + geom_smooth()
#    }
#
#    return(p)
# }


# heat map
gplt.heat <- function(data, x, y, z, w='NULL', ...) {

   # geom_tile can be directly applied without calling
   # ddply explicitly; however, slow and color scale
   # limits are based on original, not summarized data

   x.data <- data[[x]]
   y.data <- data[[y]]
   z.data <- data[[z]]

   # compute (weighted) means

   means <- group.center(data, z, c(x, y), w=w, method='mean')
   names(means) <- c(x, y, z)

   p <- ggplot(means, aes_string(x=x, y=y, fill=z)) +
      geom_tile(color=NA, ...) +
      scale_fill_gradientn(colours=rev(grDevices::rainbow(2)))
   return(p)
}


# spine plot of two variables
# note: We could use mosaicplot, but we want to consistentl stick with ggplot
#    another option would be the prodplots package, but currently seems
#    to be in development stage
#
# - plot.margin.x - horizontal gap between rectangles for x-values
# - no gaps between y-rectangles
# - exclude-factor - flag to apply in factor(...) calls, to preserve NA values
#   if required
gplt.spine <- function(data, x, y, w='NULL',
                       plot.margin.x=0.05,
                       exclude.factor=NA,
                       ...) {

   for (n in c(x, y)) {
      if (is.numeric(data[[n]])) {
         data[[n]] <- quantize(data, n, w=w)
      }
   }

   x.lev <- length(levels(data[[x]]))
   y.lev <- length(levels(data[[y]]))

   xy.tab  <- marginalize(data, c(x, y), w, exclude.factor=exclude.factor)
   x.tab   <- marginalize(data, x, w, exclude.factor=exclude.factor)
   y.cond  <- sweep(xy.tab, 1, x.tab, FUN='/')

   x.tab.df             <- as.data.frame(x.tab, responseName='x.mrg')
   x.tab.df$x.mrg.cum   <- c(0, cumsum(x.tab.df$x.mrg)[1:(x.lev-1)])
   x.tab.df$x.cnt       <- seq(0, x.lev-1)
   x.tab.df$x.center    <- with(x.tab.df, x.mrg.cum  + x.mrg/2 + x.cnt * plot.margin.x)

   y.cond.df            <- as.data.frame(aperm(y.cond), responseName='y.cond')
   y.cond.df$y.cond.cum <- as.vector(apply(y.cond, 1, function(x) {c(0, cumsum(x)[1:(y.lev-1)])}))

   plot.data <- y.cond.df
   plot.data <- merge(plot.data, x.tab.df)

   plot.data$left     <- plot.data$x.mrg.cum  + plot.data$x.cnt * plot.margin.x
   plot.data$right    <- plot.data$left       + plot.data$x.mrg
   plot.data$bottom   <- plot.data$y.cond.cum
   plot.data$top      <- plot.data$bottom     + plot.data$y.cond

   ggplot(plot.data, aes(xmin=left, xmax=right, ymin=bottom, ymax=top)) +
      # make outline color same as fill color - black outline will hide colors
      # for very narrow stripes
      geom_rect(aes_string(fill=y, col=y), ...) +
      scale_x_continuous(breaks=x.tab.df$x.center, labels=levels(data[[x]])) +
      theme(axis.ticks.x=element_blank())
}


# spine plot of three variables
# - plot.margin.x - horizontal gap between rectangles for x-values
# - plot.margin.y - vertical gap between rectangles for y-values
# - no gaps between z-rectangles
# - exclude-factor - flag to apply in factor(...) calls, to preserve NA values
#   if required
gplt.spine3 <- function(data, x, y, z, w='NULL',
                        plot.margin.x=0.05,
                        plot.margin.y=0.02,
                        exclude.factor=NA, ...) {

   for (n in c(x, y, z)) {
      if (is.numeric(data[[n]])) {
         data[[n]] <- quantize(data, n, w=w)
      }
   }

   x.lev <- length(levels(data[[x]]))
   y.lev <- length(levels(data[[y]]))
   z.lev <- length(levels(data[[z]]))

   xyz.tab <- marginalize(data, c(x, y, z), w, exclude.factor=exclude.factor)
   xy.tab  <- marginalize(data, c(x, y), w, exclude.factor=exclude.factor)
   x.tab   <- marginalize(data, x, w, exclude.factor=exclude.factor)
   y.tab   <- marginalize(data, y, w, exclude.factor=exclude.factor)
   y.cond  <- sweep(xy.tab, 1, x.tab, FUN='/')
   z.cond  <- sweep(xyz.tab, 1:2, xy.tab, FUN='/')

   x.tab.df             <- as.data.frame(x.tab, responseName='x.mrg')
   x.tab.df$x.mrg.cum   <- c(0, cumsum(x.tab.df$x.mrg)[1:(x.lev-1)])
   x.tab.df$x.cnt       <- seq(0, x.lev-1)
   x.tab.df$x.center    <- with(x.tab.df, x.mrg.cum  + x.mrg/2 + x.cnt * plot.margin.x)

   y.tab.df             <- as.data.frame(y.tab, responseName='y.mrg')
   y.tab.df$y.mrg.cum   <- c(0, cumsum(y.tab.df$y.mrg)[1:(y.lev -1)])
   y.tab.df$y.cnt       <- seq(0, y.lev-1)
   y.tab.df$y.center    <- with(y.tab.df, y.mrg.cum  + y.mrg/2 + y.cnt * plot.margin.y)

   y.cond.df            <- as.data.frame(aperm(y.cond), responseName='y.cond')
   y.cond.df$y.cond.cum <- as.vector(apply(y.cond, 1, function(x) {c(0, cumsum(x)[1:(y.lev-1)])}))

   z.cond.df            <- as.data.frame(aperm(z.cond), responseName='z.cond')
   z.cond.df$z.cond.cum <- as.vector(apply(z.cond, 2:1, function(x) {c(0, cumsum(x)[1:(z.lev-1)])}))

   plot.data <- z.cond.df
   plot.data <- merge(plot.data, y.cond.df)
   plot.data <- merge(plot.data, x.tab.df)
   plot.data <- merge(plot.data, y.tab.df)

   plot.data$left     <- plot.data$x.mrg.cum  + plot.data$z.cond.cum * plot.data$x.mrg + plot.data$x.cnt * plot.margin.x
   plot.data$right    <- plot.data$left       + plot.data$z.cond * plot.data$x.mrg
   plot.data$bottom   <- plot.data$y.cond.cum + plot.data$y.cnt * plot.margin.y
   plot.data$top      <- plot.data$bottom     + plot.data$y.cond

   ggplot(plot.data, aes(xmin=left, xmax=right, ymin=bottom, ymax=top)) +
      # make outline color same as fill color - black outline will hide colors
      # for very narrow stripes
      geom_rect(aes_string(fill=z, col=z), ...) +
      scale_y_continuous(breaks=y.tab.df$y.center, labels=levels(data[[y]])) +
      scale_x_continuous(breaks=x.tab.df$x.center, labels=levels(data[[x]])) +
      # write labels slanted to mitigate clutter
      theme(axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank())
}

# bar plot with stat_bin
gplt.bar <- function(data, x, w='NULL', ...) {
   p <- ggplot(data, aes_string(x=x, weight=w)) +
      geom_bar(...)
   return(p)
}

# bar plot with stat_identity
gplt.bar.identity <- function(data, x, y, z='NULL', w='NULL', ...) {
   p <- ggplot(data, aes_string(x=x, y=y, color=z, weight=w)) +
      geom_bar(stat='identity', position=position_dodge(), ...)
   return(p)
}

# violin (or box) plot with overlayed median points
gplt.box <- function(data, x, y, w='NULL', med='NULL', use.geom.violin=TRUE, ...) {

   # note: geom_violin throws an error if all values are equal
   # workaround: add very small jitter
   if (use.geom.violin) {
      data[[y]] <- data[[y]] +
         1E-10 * max(data[[y]], na.rm=TRUE) * (runif(nrow(data)) - 0.5)
   }
   p <- ggplot(data, aes_string(x=x, y=y, weight=w, ymax=max(y,na.rm=TRUE)))

   if (use.geom.violin) {
      p <- p + geom_violin(scale='width', ...)
      if (med != 'NULL') {
         # dodge does not work correctly when width is not specified
         # see https://github.com/hadley/ggplot2/issues/525
         p <- p + geom_point(mapping=aes_string(y=med),
                             position=position_dodge(width=0.9),
                             size=2, shape=1)
      }
   } else {
      p <- p + geom_boxplot(...)
   }

   return(p)
}


# density plot with overlayed vertical median lines
gplt.density <- function(data, x, w='NULL', med='NULL',
                         use.geom.density=TRUE,
                         n.breaks.histogram=NA,
                         ...) {
   p <- ggplot(data, aes_string(x=x, weight=w))
   if (use.geom.density) {
      p <- p + geom_density(aes(y=..scaled..), adjust=0.5, trim=TRUE, na.rm=TRUE, ...)

      # geom_vline cannot be used, see https://github.com/hadley/ggplot2/issues/426
      if (med != 'NULL') {
         p <- p + geom_vline(aes(xintercept=.center.), lwd=1)
      }

      # density numbers themselves not very meaningful
      p <- p + theme(axis.ticks.y=element_blank(),
                     axis.text.y=element_blank())

   } else {
      if (is.numeric(n.breaks.histogram)) {
         n.breaks <- n.breaks.histogram
      } else {
         n.breaks <- nclass.FD.modified(data[[x]])
      }
      bin.width <- diff(range(data[[x]], na.rm=TRUE))/n.breaks
      p <- p + geom_histogram(binwidth=bin.width, position='identity', ...)
   }

   return(p)
}


# hexbin plot with overlayed smoothing line
gplt.hex <- function(data, x, y, w='NULL', ...) {
   p <- ggplot(data, aes_string(x=x, y=y, weight=w)) +
      geom_hex(...) +
      geom_smooth() +
      # log scaling of color often reveals more details
      scale_fill_gradientn(colours=rev(grDevices::rainbow(2)), trans='log')
   return(p)
}


# blank plot with optional message text
gplt.blank <- function(text=NULL, ...) {
   p <- ggplot(data.frame(x=c(0,1), y=c(0,1)), aes(x=x, y=y)) +
      geom_blank(...) +
      labs(x=NULL, y=NULL) +
      xlim(0,1) +
      ylim(0,1) +
      theme(axis.ticks.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.y=element_blank(),
            panel.background=element_blank())

   if (!is.null(text)) {
      p <- p + annotate('text', x=0.5, y=0.5, label=text, hjust=0.5, vjust=0.5)
   }
   return(p)
}


#'Create options structure for \code{plotluck}
#'
#'@param ... parameters to override default settings
#'@return a named list of options, usable by function \code{plotluck}
#'
#'  \code{plotluck} accepts a list of options to modify its behavior. Calling
#'  \code{plotluck.options} without arguments produces a list with the default
#'  values. Specifying any number of attribute/value pairs overrides these
#'  values selectively.
#'
#'  Available options and their meaning is described inline in the documentation
#'  of \code{link{plotluck}}.
#'
#'@section Shortcut options:
#'\itemize{
#'\item By default, hex, density, box, and violin plots revert to
#'scatter plots for low number of points. \code{prefer.scatter=FALSE} (resp.
#'\code{TRUE}) sets \code{min.points.hex}, \code{min.points.density}, and
#'\code{min.points.box} to zero (resp. \code{Inf}).
#'\item If variable \code{z} is specified and the plot type, as well as its
#'number of levels, allows it, coloring is used to for representation in a single
#'plot; otherwise, \code{plotluck} reverts to faceting. Users can influence this
#'behavior by setting \code{prefer.color.for.z=TRUE} (resp. \code{FALSE}), which
#'causes \code{max.colors.scatter}, \code{max.colors.density},
#'\code{max.colors.box}, and \code{max.colors.bar} to be \code{Inf} (resp. 0).
#'}
#'
#'@note \code{plotluck}'s aim is to provide a function that is usable
#'  "out-of-the-box", with no or very little manual tweaking. If you find
#'  yourself needing to change option values repeatedly or find the presets to
#'  be suboptimal, please contact the author.
#'
#'@seealso \code{\link{plotluck}}
#'@export
#' @examples
#' # list all default options
#' plotluck.options()
#'
#' data(iris)
#' # default with violin plot
#' plotluck(iris, Species, Petal.Length)
#'
#' # use box-and-whiskers plot instead
#' plotluck(iris, Species, Petal.Length, opts=plotluck.options(use.geom.violin=FALSE))
#'
#'@export
plotluck.options <- function(...) {
   opts <- list(
      max.factor.levels=100,
      few.unique.as.factor=10,
      discretize.intervals.z=5,
      min.points.hex=5000,
      min.points.density=20,
      min.points.box=10,
      min.size.grid.heat=2,
      min.coverage.heat=0.5,
      convert.duplicates.to.weights=TRUE,
      min.points.jitter=3,
      width.jitter=0.01,
      height.jitter=0.01,
      factor.jitter=0.25,
      max.colors.scatter=3,
      max.colors.density=3,
      max.colors.box=3,
      max.colors.bar=3,
      trans.log.thresh=2,
      spine.plot.margin.x=0.05,
      spine.plot.margin.y=0.02,
      max.sample.rows=100000,
      use.geom.violin=TRUE,
      max.levels.violin=20,
      use.geom.density=TRUE,
      n.breaks.histogram=NA,
      max.facets.column=10,
      max.facets.row=10,
      exclude.factor=NULL,
      fill.default='deepskyblue',
      alpha.default=0.3,
      theme.axis.x.factor=theme(axis.text.x=element_text(angle=-45, hjust=0, vjust=1)))
   overrides <- list(...)

   # shortcut options
   if ('prefer.color.for.z' %in% names(overrides)) {
      use.color.for.z <- overrides[['prefer.color.for.z']]
      overrides['prefer.color.for.z'] <- NULL
      if (!is.na(use.color.for.z)) {
         if (use.color.for.z) {
            opts$max.colors.scatter <- Inf
            opts$max.colors.density <- Inf
            opts$max.colors.box     <- Inf
            opts$max.colors.bar     <- Inf
         } else {
            opts$max.colors.scatter <- 0
            opts$max.colors.density <- 0
            opts$max.colors.box     <- 0
            opts$max.colors.bar     <- 0
         }
      }
   }
   if ('prefer.scatter' %in% names(overrides)) {
      prefer.scatter <- overrides[['prefer.scatter']]
      overrides['prefer.scatter'] <- NULL
      if (!is.na(prefer.scatter)) {
         if (prefer.scatter) {
            opts$min.points.hex <- Inf
            opts$min.points.density <- Inf
            opts$min.points.box     <- Inf
         } else {
            opts$min.points.hex <- 0
            opts$min.points.density <- 0
            opts$min.points.box     <- 0
         }
      }
   }

   unknown <- setdiff(names(overrides), names(opts))
   if (length(unknown) > 0) {
      stop(sprintf('Unknown options: %s', paste(unknown, sep='',collapse =', ')))
   }
   opts[names(overrides)] <- overrides
   return(opts)
}


sample.data <- function(data, w='NULL', max.rows) {
   n.row <- nrow(data)
   if (n.row > max.rows) {
      warning(sprintf('Data set has %d rows, sampling down to %d rows', n.row, max.rows))
      if (w == 'NULL') {
         data <- data[sample(1:n.row, max.rows),, drop=FALSE]
      } else {
         # note: weighted sampling itself can be quite slow
         data <- data[sample(1:n.row, max.rows, prob=data[[w]]),, drop=FALSE]
      }
   }
   return(data);
}


#' "I'm feeling lucky" for ggplot
#'
#' The aim of \code{plotluck} is to let the user focus on \emph{what} to plot,
#' and automate the \emph{how}. It examines the data characteristics of one,
#' two, or three variables and accordingly creates a scatter, box, bar, density,
#' hex or spine plot, or a heat map. It also automates handling of observation
#' weights, log-scaling of axes, reordering of factor levels, and overlays of
#' smoothing curves and median lines.
#'
#' @param data a data frame
#' @param x,y,z column names
#' @param w weight column (optional)
#' @param opts a named list of options (optional)
#' @param ... additional parameters to be passed to the respective geom_* objects
#' @return a ggplot object
#' @export
#'@keywords hplot, aplot, dplot
#'@concept automation
#'@concept visualization
#'@concept plotting
#'@concept exploratory data analyis
#'@concept ggplot
#'@concept ggplot2
#'@concept heat map
#'@concept density plot
#'@concept violin plot
#'@concept hexbin
#'@concept histogram
#'@concept bar plot
#'@concept box plot
#'@concept spine plot
#'@concept scatter plot
#'@concept heat map
#'
#' @seealso \code{\link{ggplot}}, \code{\link{plotluck.options}}, \code{\link{plotluck.multi}}
#'
#' @section Determining the type of plot: One-dimensional plots (i.e., when
#'   \code{y} and \code{z} are \code{NULL}) are either density plots (for
#'   numeric variables) or bar plots (for factors). In the former case, if the
#'   number of data points is not sufficient (as determined by the option
#'   \code{min.points.density}) a strip chart (1D scatter plot) is drawn.
#'
#'   The following table summarizes the general heuristics for \emph{two} dimensions
#'   (\code{x, y} are not \code{NULL}):
#'
#' \tabular{lll}{
#' \strong{x}\tab \strong{y}\tab \strong{type}\cr
#' num\tab num\tab hexbin\cr
#' num\tab fact\tab density (histogram)\cr
#' fact\tab num\tab violin (box-and-whisker)\cr
#' fact\tab fact\tab spine\cr
#' }
#'
#' In the first three cases, if the number of data points is not sufficient (as
#' determined by the options \code{min.points.hex}, \code{min.points.density},
#' or \code{min.points.box}), the plot switches to a scatter plot.
#'
#' Density plots come with an overlaid vertical (weighted) median line; and
#' numeric scatter plots with a smoothing line including standard deviations.
#
#' For three dimensions, the algorithm first tries sue coloring to fit the
#' \code{z} dimension into a single graph. If this is not possible (e.g., in the
#' case of the hexbin plot, which encodes the number of points in a bin by the
#' color; or as determined by the options \code{max.colors.scatter},
#' \code{max.colors.density}, \code{max.colors.box}, or \code{max.colors.bar}),
#' facetting is used, in a direction (i.e., horizontally or vertically subject
#' to \code{max.facets.row} and \code{max.facets.column})  that facilitates
#' comparing the particular type of plot. To this end, a numerical \code{z}
#' variable is discretized into \code{discretize.intervals.z} many intervals of
#' equal count or weight.
#'
#'There are three notable exceptions from this rule:
#'\itemize{
#'\item If all three variables are categorical, a spine plot is drawn.
#'\item If the data is "grid-like", a heat map is drawn with the mean value of \code{z} providing the color.
#' The data is considered grid-like if
#' \itemize{
#' \item \code{x} and \code{y} are numeric or ordered factors, and \code{z} is numeric.
#' \code{x} and \code{y} have at least \code{min.size.grid.heat} distinct values each.
#' \item A sufficient fraction (as determined by option
#' \code{min.coverage.heat}) of the \code{(x,y)} grid points has at least one
#' data point.
#' }
#' \item If each factor combination occurs at most once in the data set, we resort to
#' bar plots.
#' }
#'
#'@section Reordering of factor levels: To better illustrate the relation
#'  between an independent factor variable and a dependent numerical variable
#'  (or an ordered factor), factor levels are reordered according to the
#'  (weighted) mean of the dependent variable. Typically, the dependent variable
#'  is \code{y}, but it is \code{z} for spine plots and heat maps. If there is
#'  no dependent numeric variable, ordering proceeds by frequency.
#'
#'  \code{plotluck} respects \emph{ordered} factors, so to prevent undesirable
#'  reordering, set the data type to \code{ordered}.
#'
#'@section Instance weights: Argument \code{w} allows to specify weights or
#'  frequency counts for each row of data. All plots and summary lines take this
#'  weight into account (well, unfortunately all except hexbin, as the
#'  underlying implementation does not support it. We could approximate it by
#'  resampling with replacement, but this seems a little crude).
#'
#'  In scatter plots, weights are indicated by a shaded disk whose area is
#'  proportional to the weight.
#'
#'  There is some redundancy with regard to \emph{repeated points}, since they
#'  are conceptually identical to unique points with an associated frequency
#'  weight. By default, the option \code{convert.duplicates.to.weights} is
#'  switched on, which does exactly this. If it is \code{FALSE}, horizontal and
#'  vertical jittering will be applied if the number of duplicated points
#'  exceeds \code{min.points.jitter}. The amount of jittering can be controlled
#'  with \code{width.jitter}, \code{height.jitter}, and \code{factor.jitter}.
#'
#'@section Axis scaling: \code{plotluck} supports logarithmic and log-modulus
#'  axis scaling. log-modulus is considered if values are both positive and
#'  negative; in this case, the transform function is \code{f(x) = sign(x) *
#'  log(1+abs(x))}.
#'
#'  The decision whether to apply scaling is based on the proportion of total
#'  display range that is occupied by the 'core' region of the distribution
#'  (between the lower and upper quartiles); namely, if the transform could
#'  enlarge this region by a factor of at least \code{trans.log.thresh}.
#'
#'@section Missing values: Typically, missing (\code{NA}) values in factors are
#'  ignored. To make these explicit, create a separate level using the
#'  \code{exclude} option in \code{\link{factor}}. In this case, \code{plotluck}
#'  will respect these levels. You can also explicitly set the exclude handling
#'  using the option \code{exclude.factor}, with the same semantics as in
#'  function \code{factor()}.
#'
#'  Missing values in numeric variables are ignored.
#'
#'@section Sampling: For very large data sets, plots can take a very long time
#'  (or even crash R). \code{plotluck} has a built-in stop-gap: If the data
#'  comprises more than \code{max.sample.rows}, it will be sampled down to that
#'  size (taking weights into account, if supplied).
#'
#'@section Factor preprocessing: Frequently, when numeric variables have very
#'  few values despite sufficient data size, it helps to treat these values as
#'  the levels of a factor; this is governed by option
#'  \code{few.unique.as.factor}.
#'
#'  If an unordered factor has too many levels, plots can get messy. In this
#'  case, only the \code{max.factor.levels} most frequent ones are retained,
#'  while the rest are merged into a default level '.other.'.
#'
#'@section Column name matching: Column names \code{x, y, z, w} are matched by
#'  unique prefix, and ignoring case.
#'
#'@section Remarks on the choice of plot types: By default, \code{plotluck} uses
#'  violin and density plots in place of the more traditional box-and-whisker
#'  plots and histograms; these modern graph types convey the shape of a
#'  distribution better. In the former case, summary statistics like mean and
#'  quantiles are less useful if the distribution is not unimodal; a wrong
#'  choice of the number of bins of a histogram can create misleading artifacts.
#'
#'  If the resulting graph would contain too many (more than
#'  \code{max.levels.violin}) violin plots in a row, the algorithm switches to
#'  box plots. The defaults can also be directly overriden by changing options
#'  \code{use.geom.violin} and \code{use.geom.density}. The number of bins of a
#'  histogram can be customized with \code{n.breaks.histogram}. The default
#'  setting, \code{NA}, applies a heuristic estimate.
#'
#'  Due to their well-documented problematic aspects, pie charts and stacked bar
#'  graphs are not supported.
#'
#'  Cleveland's Dot Plots (\code{\link{dotchart}}) can be produced as a special
#'  case of scatter plots.
#'
#'  With real-world data (as opposed to smooth mathematical functions),
#'  three-dimensional scatter, surface, or contour plots can often be hard to
#'  read if the shape of the distribution is not suitable, data coverage is
#'  uneven, or if the perspective is not carefully chosen depending on the data.
#'  Therefore, we have refrained from incorporating them.
#'
#'@section Remarks on the use of options: For completeness, we have included the
#'  description of option parameters in the current help page. However, the
#'  tenet of this function is to be usable "out-of-the-box", with no or very
#'  little manual tweaking required. If you find yourself needing to change
#'  option values repeatedly or find the presets to be suboptimal, please
#'  contact the author.
#
#'@section What \code{plotluck} does not: This function is designed for generic
#'  out-of-the-box plotting, and not suitable to produce more specialized types
#'  of plots that arise in specific application domains (e.g., association,
#'  stem-and-leaf, star plots, geo maps, etc). It is restricted to at most three
#'  variables. Parallel plots with variables on different scales (such as time
#'  series of multiple related signals) are not supported.
#'
#' @examples
#' #1D
#' #density plot
#' data(iris)
#' plotluck(iris, Petal.Length)
#'
#'  # bar plot
#' data(movies, package='ggplot2movies')
#' plotluck(movies, mpaa)
#'
#' # 2D scatter plot
#' plotluck(iris, Petal.Length, Petal.Width)
#' plotluck(movies,votes,rating)
#'
#' # density
#' plotluck(movies, rating, mpaa)         # using facets
#'
#' data(diamonds, package='ggplot2')
#' plotluck(diamonds, price, cut)
#'
#' plotluck(iris, Sepal.Length, Species)  # using colors
#'
#' # box/violin plot
#' plotluck(movies, mpaa, rating)
#' plotluck(movies, Documentary, budget) # implicit conversion of binary variables
#'
#' # spine plot
#' data(Titanic)
#' plotluck(as.data.frame(Titanic), Class, Survived, w=Freq)
#'
#' data(occupationalStatus)
#' df <- as.data.frame(occupationalStatus)
#' df$origin <- ordered(df$origin)
#' df$destination <- ordered(df$destination)
#' plotluck(df, origin, destination, w=Freq)
#'
#' # 3D
#' # heat map
#' plotluck(diamonds, cut, color, price)
#'
#' # scatter plat
#' plotluck(movies, length, rating, mpaa)             # using facets
#' plotluck(iris, Petal.Length, Petal.Width, Species) # using colors
#'
#' # density plot
#' plotluck(diamonds, price, cut, color)
#'
#' # bar plot
#' plotluck(movies, mpaa, rating, Action) # grouped
#'
#' # spine plot
#' plotluck(as.data.frame(Titanic), Class, Sex, Survived, w=Freq)

plotluck <- function(data, x, y=NULL, z=NULL, w=NULL,
                     opts=plotluck.options(),
                     ...) {

   x <- deparse(substitute(x))
   y <- deparse(substitute(y))
   z <- deparse(substitute(z))
   w <- deparse(substitute(w))

   # match partial column names, if unambiguous
   cols.lc <- tolower(names(data))

   for (n in c('x', 'y', 'z', 'w')) {
      n.val <- tolower(get(n))
      if (n.val != 'null') {
         idx <- pmatch(tolower(n.val), cols.lc)
         if (is.na(idx)) {
            # if several are matching, find all of them
            idx <- which(!is.na(sapply(cols.lc, function(tab, x) pmatch(x, tab), x=n.val)))

            if (length(idx) == 0) {
               stop(sprintf('No match for "%s" in data columns', n.val))
            } else {
               stop(sprintf('Ambiguous match for "%s": %s', n.val, paste(names(data)[idx], sep='', collapse=',')))
            }
         } else {
            assign(n, names(data)[idx])
            if (all(is.na(data[[idx]]))) {
               stop(sprintf('Variable %s is completely missing, giving up', names(data)[idx]))
            }
         }
      }
   }

   if (y == 'NULL' && z != 'NULL') {
      stop('z can only be specified together with y')
   }

   if (y == x) {
      y <- 'NULL'
   }

   if (z %in% c(x, y)) {
      z <- 'NULL'
   }

   vars.non.null   <- unique(c(x, y, z)[c(x, y, z) != 'NULL'])
   vars.w.non.null <- unique(c(x, y, z, w)[c(x, y, z, w) != 'NULL'])

   # reminder: all data manipulations have to be done before ggplot is called

   data <- data[,vars.w.non.null, drop=FALSE]

   # process weights
   if (w != 'NULL') {
      if  (!is.numeric(data[[w]])) {
         stop('Weight must be numeric')
      }
      if  (any(data[[w]]<0)) {
         stop('Weight must be non-negative')
      }
      weight.na <- is.na(data[[w]])
      if  (any(weight.na)) {
         warning('Weight is NA for %d instances, deleting', length(which(weight.na)))
         data <- data[!weight.na,]
      }
      # if weights are integer, Hmisc::wtd.quantile() can lead to NA due to overflow
      data[[w]] <- as.double(data[[w]])
   }

   # if data size too large, apply sampling
   data <- sample.data(data, w, opts$max.sample.rows)

   # technicality: if the user has created a factor with 'exclude=FALSE', i.e.
   # is including NA as a separate level, we have to preserve it in all subsequent
   # factor(..) transformations.
   # For simplicity, we assume this is the same for all factors in the data set.
   # Unfortunately, we have to pass this option down to all functions that internally
   # manipulate factors.

   exclude.factor <- opts$exclude.factor
   if (is.null(exclude.factor)) {
      exclude.factor <- !any(sapply(1:length(data),
                                    function(x) {is.factor(data[[x]]) && any(is.na(levels(data[[x]])))}))
      if (exclude.factor) {
         exclude.factor <- NA
      }
   }

   trans       <- list('NULL'=scales::identity_trans())
   is.num      <- list('NULL'=FALSE)
   is.ord      <- list('NULL'=FALSE)
   uniq        <- list('NULL'=0)
   non.na      <- list('NULL'=FALSE)
   title.label <- list('NULL'=NULL)

   for (n in vars.non.null) {
      data[[n]] <- preprocess.factors(data, n, w,
                                     max.factor.levels=opts$max.factor.levels,
                                     few.unique.as.factor=opts$few.unique.as.factor,
                                     exclude.factor=exclude.factor)
      trans[[n]]       <- get.trans.fun(data[[n]], opts$trans.log.thresh)
      title.label[[n]] <- format.label.trans(n, trans[[n]])
      is.num[[n]]      <- is.numeric(data[[n]])
      is.ord[[n]]      <- is.numeric(data[[n]]) || is.ordered(data[[n]])
      uniq[[n]]        <- length(unique(data[[n]]))
   }

   title.label[[z]] <- NULL
   vars.numeric     <- vars.non.null[unlist(is.num[vars.non.null])]
   vars.non.numeric <- setdiff(vars.non.null, vars.numeric)

   # special case heat map
   grid.xy <-is.grid.like(data, x, y, z,
                          min.size=opts$min.size.grid.heat,
                          min.coverage=opts$min.coverage.heat)

   # order factor levels to best reflect dependent variable

   order.by <- structure(list(x, y, z), names=list(x, y, z))

   if (z != 'NULL' &&
          (grid.xy ||
              ((!is.num[[x]]) && (!is.num[[y]])))) {
      # heat map or spine plot: dependent variable is z
      order.by[[x]] <- z
      order.by[[y]] <- z
   } else {
      if (is.ord[[y]]) {
         order.by[[x]] <- y
      }
      if (is.ord[[x]]) {
         order.by[[y]] <- x
      }

      # if possible, order z by numeric variable, preferably y
      if (is.ord[[y]]) {
         order.by[[z]] <- y
      } else if (is.ord[[x]]) {
         order.by[[z]] <- x
      }
   }

   for (n in vars.non.numeric) {
      data[[n]] <- factor(data[[n]],
                          levels=order.factor.value(data, n, order.by[[n]], w, exclude.factor=exclude.factor),
                          exclude=exclude.factor)
   }

   vars.covered <- intersect(c(x, y), vars.non.null) # the variables reflected in the plot
   color.usable <- TRUE  # can we color the plot to reflect additional variables?
   flip         <- FALSE # for num/factor scatter plots, we need to apply coord_flip(), see below

   if (grid.xy) {
      p <- gplt.heat(data, x, y, z, w,
                     ...)
      type.plot <- 'heat'
      color.usable <- FALSE
      vars.covered <- c(x, y, z)
   } else {

      # for 3D, treat z as factor except for heat map
      if (is.num[[z]]) {
         data[[z]] <- quantize(data, z, w=w, max.breaks=opts$discretize.intervals.z)
         data[[z]] <- factor(data[[z]],
                             levels=order.factor.value(data, z, order.by[[z]], w,
                                                      exclude.factor=exclude.factor),
                             exclude=exclude.factor)
         trans[[z]]       <- get.trans.fun(data[[z]], opts$trans.log.thresh)
         is.num[[z]]      <- is.numeric(data[[z]])
         is.ord[[z]]      <- is.numeric(data[[z]]) || is.ordered(data[[z]])
         uniq[[z]]        <- length(unique(data[[z]]))
         vars.non.numeric <- unique(c(vars.non.numeric, z))
         vars.numeric     <- setdiff(vars.non.null, vars.non.numeric)
      }

      facet.label <- list('NULL'='NULL')
      for (n in vars.non.null) {
         data[[sprintf('facet.label.%s', n)]] <- format.facets(data, n)
      }

      # precompute medians
      # note: stat_summary(xintercept=..., geom='vline) does not work inside facets; see
      # https://groups.google.com/forum/#!topic/ggplot2/z_wHIzwSSzM

      if (length(vars.numeric) == 1) {
         grp.med <- group.center(data, vars.numeric[1], vars.non.numeric, w=w)
         data <- merge(data, grp.med)
      }

      # compute joint number of factor combinations, and
      # and the maximum/median number of instances in them
      if (length(vars.non.numeric) > 0) {
         t <- table(data[,vars.non.numeric], exclude=exclude.factor)
         num.level   <- length(t)
         n.med.level <- median(t)
         n.max.level <- max(t)
      } else {
         num.level <- 1
         n.med.level <- nrow(data)
         n.max.level <- nrow(data)
      }

      # determine type of plot

      if (is.num[[x]] && is.num[[y]]) {
         if (nrow(data) >= opts$min.points.hex) {
            type.plot <- 'hex'
            color.usable <- FALSE
            p <- gplt.hex(data, x, y, w, alpha=opts$alpha.default, ...)
         } else {
            type.plot <- 'scatter.num.num'
            p <- gplt.scatter(data, x, y, w,
                              convert.duplicates.to.weights=opts$convert.duplicates.to.weights,
                              min.points.jitter=opts$min.points.jitter,
                              w.jitter=opts$width.jitter,
                              h.jitter=opts$height.jitter,
                              ...)
         }
      } else if (is.num[[x]] && (!is.num[[y]])) {

         if (n.med.level > opts$min.points.density) {
            type.plot <- 'density'
            vars.covered <- x
            if (opts$use.geom.density == TRUE) {
               if (w == 'NULL') {
                  title.label[[y]] <- 'density'
               } else {
                  title.label[[y]] <- sprintf('%s density', w)
               }
            } else {
               # histogram
               if (w == 'NULL') {
                  title.label[[y]] <- 'count'
               } else {
                  title.label[[y]] <- w
               }
            }

            p <- gplt.density(data, x, w=w, med='.center.',
                              use.geom.density=opts$use.geom.density,
                              n.breaks.histogram=opts$n.breaks.histogram,
                              alpha=opts$alpha.default, ...)
         } else {
            # not enough points for density plot
            type.plot <- 'scatter.num.fact'
            p <- gplt.scatter(data, x, y, w,
                              convert.duplicates.to.weights=opts$convert.duplicates.to.weights,
                              min.points.jitter=opts$min.points.jitter, ...)
            # HACK: ggplot2 does not implement vertical dodging, therefore we
            # use coord_flip() as a workaround; see code for gplt.scatter
            flip <- TRUE

            tmp              <- title.label[[x]]
            title.label[[x]] <- title.label[[y]]
            title.label[[y]] <- tmp
         }
      } else if ((!is.num[[x]]) && y == 'NULL') {
         # no log transforms for bars; see http://www.perceptualedge.com/articles/b-eye/dot_plots.pdf
         type.plot <- 'bar'
         color.usable <- FALSE
         vars.covered <- x
         is.num[[y]] <- TRUE

         if (w == 'NULL') {
            title.label[[y]] <- 'count'
         } else {
            title.label[[y]] <- w
         }

         p <- gplt.bar(data, x, w, alpha=opts$alpha.default, ...)
      } else if ((!is.num[[x]]) && is.num[[y]]) {
         if (n.med.level > opts$min.points.box) {
            type.plot <- 'box'

            # if there are too many violin plots in a horizontal row, they
            # just look like black lines
            use.geom.violin <- opts$use.geom.violin
            if (use.geom.violin && num.level > opts$max.levels.violin) {
               use.geom.violin <- FALSE
            }
            p <- gplt.box(data, x, y, w=w, med='.center.',
                          use.geom.violin=use.geom.violin,
                          alpha=opts$alpha.default, ...)
         } else if (n.max.level > 1) {
            type.plot <- 'scatter.fact.num'
            p <- gplt.scatter(data, x, y, w=w,
                              convert.duplicates.to.weights=opts$convert.duplicates.to.weights,
                              min.points.jitter=opts$min.points.jitter, ...)
         } else {
            # special case: single point each - use bar plot
            # no log transform
            type.plot <- 'bar'
            color.usable <- TRUE
            vars.covered <- c(x, y)
            p <- gplt.bar.identity(data, x=x, y=y, w=w, alpha=opts$alpha.default, ...)
         }
      } else if ((!is.num[[x]]) & (!is.num[[y]])) {
         type.plot <- 'spine'
         color.usable <- FALSE
         if (z == 'NULL') {
            is.num[[y]] <- TRUE      # display y grid lines
            title.label[[y]] <- NULL # y reflected in color legend
            p <- gplt.spine(data, x, y, w,
                            plot.margin.x=opts$spine.plot.margin.x,
                            exclude.factor=exclude.factor,
                            alpha=opts$alpha.default, ...)
            vars.covered <- c(x, y)
         } else {
            p <- gplt.spine3(data, x, y, z, w,
                             plot.margin.x=opts$spine.plot.margin.x,
                             plot.margin.y=opts$spine.plot.margin.y,
                             exclude.factor=exclude.factor,
                             alpha=opts$alpha.default, ...)
            vars.covered <- c(x, y, z)
         }
      }
   }

   # y, z-coordinate: represent as color or wrap?
   vars.remaining <- setdiff(vars.non.null, vars.covered)

   color.guides <- TRUE # set to FALSE if guide was redundant

   if (type.plot == 'density') {
      if (length(vars.remaining) == 0) {
         # color 1D density plots with default color
         p <- p + aes(fill='something') +
            aes(color='something') +
            scale_fill_manual(values=opts$fill.default) +
            scale_color_manual(values=opts$fill.default)
         color.guides <- FALSE
      } else {
         # might be overwritten in 1734, but no harm in coloring in either case,
         # even for wrapping
         p <- p + aes_string(fill=y) +
            aes_string(color=y)
            color.guides <- FALSE
      }
   }

   if (type.plot %in% c('box', 'bar', 'scatter.fact.num')) {
      p <- p + aes_string(fill=x) +
         aes_string(color=x)
      color.guides <- FALSE
   }

   if (type.plot == 'scatter.num.fact') {
      p <- p + aes_string(fill=y) +
         aes_string(color=y)
      color.guides <- FALSE
   }

   color.used <- FALSE

   if (length(vars.remaining) > 0) {

      var.next <- vars.remaining[1]

      # try to use color to express the remaining variables

      if (color.usable) {

         u <- uniq[[var.next]]

         if (type.plot == 'scatter.num.num' && u <= opts$max.colors.scatter) {
            p <- p + aes_string(color=var.next)
            color.used <- TRUE
         } else if (( type.plot == 'density'          && u <= opts$max.colors.density) ||
                     (type.plot == 'scatter.num.fact' && u <= opts$max.colors.scatter) ||
                     (type.plot == 'scatter.fact.num' && u <= opts$max.colors.scatter) ||
                     (type.plot == 'box'              && u <= opts$max.colors.box) ||
                     (type.plot == 'bar'              && u <= opts$max.colors.bar)) {
            p <- p + aes_string(fill=var.next) + aes_string(color=var.next)
            color.used <- TRUE
         }
      }
   }

   if (color.used == TRUE) {
      vars.remaining <- vars.remaining[-1]
   } else if (!color.guides) {
      p <- p + guides(fill=FALSE, color=FALSE)
   }

   # the remaining variables will be faceted
   if (length(vars.remaining) == 2) {
      p <- p + opts$theme.axis.x.factor # axis text can easily overlap in small diagrams
      p <- p + facet_grid(
         as.formula(sprintf('facet.label.%s~facet.label.%s', vars.remaining[1], vars.remaining[2])))
   } else if (length(vars.remaining) == 1) {

      preferred.order <- 'none'
      if (type.plot == 'density' ||
          type.plot == 'scatter.num.fact') {
         preferred.order <- 'col'
      }
      if (type.plot == 'box' ||
          type.plot == 'bar' ||
          type.plot == 'scatter.fact.num') {
         preferred.order <- 'row'
      }
      ncol <- NULL
      nrow <- NULL
      switch <- NULL
      if (preferred.order == 'col') {
         nrow <- opts$max.facets.row
         switch <- 'y'
      } else  if (preferred.order == 'row') {
         p <- p + opts$theme.axis.x.factor  # axis text can easily overlap in small diagrams
         ncol <- opts$max.facets.column
         switch <- 'x'
      }
      p <- p + facet_wrap(as.formula(sprintf('~facet.label.%s', vars.remaining[1])), nrow=nrow, ncol=ncol, switch = switch)
   }

   # labels
   p <- p + xlab(title.label[[x]]) +
      ylab(title.label[[y]])
   if (z != 'NULL') {
      labs(title=title.label[[z]])
   }

   # for factors, write horizontal axis labels at an angle to avoid clutter
   if (type.plot %in% c('box', 'bar', 'scatter.fact.num', 'spine')) {
      p <- p + opts$theme.axis.x.factor
   }

   # note: coord_flip() changes the scales and labels along with the coordinates,
   # but not the grid orientation and axis text properties

   # axis scaling
   if (!flip) {
      if (is.num[[x]]) {
         p <- p + scale_x_continuous(trans=trans[[x]])
      }
      if (is.num[[y]]) {
         p <- p + scale_y_continuous(trans=trans[[y]])
      }
   } else {
      if (is.num[[x]]) {
         p <- p + scale_y_continuous(trans=trans[[x]])
      }
      if (is.num[[y]]) {
         p <- p + scale_x_continuous(trans=trans[[y]])
      }
   }

   # points on the grid line are hard to see; remove grid line for factors
   if (is.num[[x]]) {
      p <- p + theme(panel.grid.major.x=element_line(color='black', linetype='dotted'),
               panel.grid.minor.x=element_line(color='black',linetype='dotted'))
   } else {
      p <- p + theme(panel.grid.major.x=element_blank(),
                     panel.grid.minor.x=element_blank())
   }

   if (is.num[[y]]) {
      p <- p + theme(panel.grid.major.y=element_line(color='black', linetype='dotted'),
               panel.grid.minor.y=element_line(color='black',linetype='dotted'))
   } else {
      p <- p + theme(panel.grid.major.y=element_blank(),
                     panel.grid.minor.y=element_blank())
   }

   p <- p + theme(panel.background=element_blank(), # get rid of gray background
                  strip.background = element_blank(), # no background for facet labels
                  axis.line.x=element_line(),
                  axis.line.y=element_line(),
                  panel.grid=element_line(color='black'),
                  legend.position='bottom', # use maximum area for actual plot
                  legend.text=element_text(angle=-45), # avoid overlapping tick text
                  legend.text.align=0) # use instead of hjust/vjust (these are not doing anything above))

   return(p)
}

# return a list of conditional entropies H(target|x), for each column in the data frame
# for comparability, all variables are quantized into the same number of bins
cond.entropy.data <- function(data, target, w='NULL', exclude.factor=NA) {
   n.row <- nrow(data)
   n.levels <- floor(max(2, min(log2(n.row), n.row/10))) # heuristic
   data[[target]] <- quantize(data, target, y='NULL', w=w,
                              max.breaks=n.levels, keep.breaks=n.levels, exclude.factor=exclude.factor);

   cond.entropy.quantized <- function(data, x, y, w, n.levels, exclude.factor) {
      if (x == y) {
         return(0);
      }
      x.data <- quantize(data, x, y, max.breaks=n.levels, keep.breaks=n.levels, exclude.factor=exclude.factor);
      # note: we need to make a copy for the special case of x=w
      return(cond.entropy(cbind(data, .x=x.data), '.x', y, w, exclude.factor=exclude.factor))
   }

   sapply(names(data), function(x) cond.entropy.quantized(data=data,
                                          x=x, y=target, w=w, n.levels=n.levels, exclude.factor=exclude.factor))
}


#'Create several plots at once with \code{plotluck}
#'
#'Create a number of 1D or 2D ggplots for the columns a data frame by calling
#'\code{plotluck} repeatedly.
#'@param data a data frame
#'@param x,y column names, or \code{all}
#'@param w weight column (optional)
#'@param in.grid flag whether a grid of plots should be produced
#'@param entropy.order order dependency plots by conditional entropy
#'@param max.cols,max.rows maximum number of plots to put on one page. If
#'       necessary, multiple pages are generated.
#'@param opts a named list of \code{plotluck} options (optional)
#'@param ... additional parameters to be passed to \code{plotluck}
#'@return an object of class plotluck_multi.
#'
#'  With \code{in.grid=TRUE}, produces a grid of plots, rendered as minimal
#'  "spark lines" without annotations; otherwise, opens full plots in a separate
#'  window each.
#'
#'  With \code{x=all} and \code{y=NULL}, a 1D plot for the distribution of each
#'  variable is drawn. When \code{x} (resp. \code{y}) is set to a column name
#'  and \code{y=all} (resp. \code{x=all}), a 2D plot is created with each
#'  variable in \code{data} on the \code{y}-axis (resp. \code{x}-axis), in turn.
#'  Finally, for \code{x=y=all}, a chart for each pair of variables is drawn,
#'  with a diagonal of 1D distribution plots; this is analogous to the behavior
#'  of the default plot method for data frames, see
#'  \code{\link{plot.data.frame}}.
#'
#'  With \code{entropy.order=TRUE}, when one target variable \code{t} is drawn
#'  against all other columns \code{z}, the plots are sorted by an estimate of
#'  empirical conditional entropy \code{H(t|z)}; the goal is to prioritize the
#'  more predictive independent variables. For large data sets the calculation
#'  can sometimes be time consuming; it can be suppressed by setting
#'  \code{entropy.order=FALSE}. Entropy ordering is never applied for 1-D
#'  distributions or a complete matrix of all variable pairs.
#'
#'@note The class \code{plotluck_multi} does not have any functionality; its
#'  sole purpose is to make this function work in the same way as \code{ggplot}
#'  and \code{plotluck}, namely, do the actual drawing if and only if the return
#'  value is not assigned.
#'
#'@seealso \code{\link{plotluck}}
#'@export
#' @examples
#' data(iris)
#' # All 1D distributions
#' plotluck.multi(iris)
#'
#' # 2D dependencies with one fixed variable on horizontal axis
#' plotluck.multi(iris, Species)
#'
#' # 2D dependencies with one fixed variable on vertical axis
#' plotluck.multi(iris, all, Species)
#'
#' # All pairs of variables
#' plotluck.multi(iris, all, all)

plotluck.multi <- function(data, x=NULL, y=NULL, w=NULL,
                           in.grid=TRUE, entropy.order=TRUE,
                           max.rows=10, max.cols=10,
                           opts=plotluck.options(),
                           ...) {

   x <- deparse(substitute(x))
   y <- deparse(substitute(y))
   w <- deparse(substitute(w))
   main <- deparse(substitute(data))

   if ((x == 'all' || x == 'NULL') && y == 'NULL') {
      # 1D
      x <- 'all'
   } else if (x == 'NULL' && y=='all') {
      # switch
      x <- 'all'
      y <- 'NULL'
   } else if (x != 'NULL' && x != 'all' && y != 'NULL' && y != 'all') {
      stop('At least one of x or y must be "all"')
   } else {
      if (x == 'NULL') {
         x <- 'all'
      }
      if (y == 'NULL') {
         y <- 'all'
      }
   }

   # if data size too large, apply sampling here; expensive to repeat
   data <- sample.data(data, w, opts$max.sample.rows)

   # if exactly one variable is specified, order plots by conditional entropy
   if (entropy.order &&
          length(which(c(x, y) == 'NULL')) == 0 &&
          length(which(c(x, y) == 'all')) == 1) {
      if (y == 'all') {
         target <- x
      } else {
         target <- y
      }
      cond.ent <- cond.entropy.data(data, target, w=w)
      data <- data[,order(cond.ent)]
   }

   vars <- list()

   for (n in c('x','y')) {
      n.val <- get(n)
      if (n.val == 'all') {
         vars[[n]] <- names(data)
      } else {
         vars[[n]] <- n.val
      }
   }

   combinations <- expand.grid(vars[['x']], vars[['y']], stringsAsFactors=FALSE)
   names(combinations) <- c('x','y')
   if (x != 'NULL') {
      combinations$xlab <- sprintf(' + xlab("%s")', combinations$x)
   } else {
      combinations$xlab <- ''
   }
   if (y != 'NULL') {
      combinations$ylab <- sprintf(' + ylab("%s")', combinations$y)
   } else {
      combinations$ylab <- ''
   }

   # try to make a square layout
   cols <- ceiling(sqrt(nrow(combinations)))
   rows <- ceiling(nrow(combinations)/cols)

   suppress.xlab <- FALSE
   suppress.ylab <- FALSE

   if (!in.grid) {
      theme.multi <- ''
   } else {
      # grid of small 'sparklines':
      # remove all annotation text
      theme.multi <- " + theme_minimal() +
         theme(legend.position='none',
         axis.ticks.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.y=element_blank(),
         axis.text.y=element_blank(),
         strip.background=element_blank(),
         strip.text=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         axis.title.x=element_text(size=rel(0.7)),
         axis.title.y=element_text(size=rel(0.7)))"

      theme.multi <- gsub('\\s','', theme.multi)

      # slow and hardly visible
      opts$convert.duplicates.to.weights <- FALSE

      # use the limited space to make subgraphs
      # relatively rectangular
      opts$max.facets.row    <- NULL
      opts$max.facets.column <- NULL

      # does everything fit on one page?
      is.square <- TRUE
      if (cols > max.cols || cols > max.rows) {
         is.square <- FALSE
         cols <- max.cols
         rows <- min(max.rows, ceiling(nrow(combinations) / cols))
      }
      if (x == 'all' && y == 'all' && is.square) {
         # if the full cross product fits on one page,
         # write the axis labels only on the margins
         suppress.xlab <- 'margin'
         suppress.ylab <- 'margin'
      }
      if (x != 'all') {
         # do not repeat the axis label for the constant dimension
         if (x != 'NULL') {
            main <- x
            suppress.xlab <- 'margin'
         }
         else {
            # distribution/density plot
            suppress.xlab <- 'all'
         }
      }
      if (y != 'all') {
         if (y != 'NULL') {
            main <- y
            suppress.ylab <- 'margin'
         } else {
            suppress.ylab <- 'all'
         }
      }
   }

   call.strs <- sprintf('plotluck(data,x=%s, y=%s, w=%s, opts=opts, ...)%s%s%s',
      combinations$x, combinations$y, w, combinations$xlab, combinations$ylab, theme.multi)

   call.str <- sprintf('list(%s)', paste(call.strs, collapse=','))
   structure(list(plots=eval(parse(text=call.str)),
                  in.grid=in.grid,
                  cols=cols,
                  rows=rows,
                  suppress.xlab=suppress.xlab,
                  suppress.ylab=suppress.ylab,
                  main=main),
             class='plotluck_multi')
}

# plot multiple graphs in a grid layout, possibly over multiple pages
mplot <- function(plots, rows=ceiling(sqrt(length(plots))),
                  cols=ceiling(sqrt(length(plots))),
                  suppress.xlab=FALSE, suppress.ylab=FALSE) {

   num.plots <- length(plots)
   if (cols > num.plots) {
      cols <- num.plots
      rows <- 1
   }

   size.page <- rows * cols

   layout <- matrix(seq(1, size.page),
               ncol=cols, nrow=rows, byrow=TRUE)

   for (p in 1:ceiling(num.plots/size.page)) {

      grid::grid.newpage()
      grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

      # Make each plot, in the correct location
      for (i in 1:min(size.page, num.plots-(p-1)*size.page)) {
         # Get the i,j matrix positions of the regions that contain this subplot
         matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
         plot.idx <- (p-1)*size.page+i
         plot.current <- plots[[plot.idx]]

         # suppress.*lab == 'margin' means:
         # no x-axis labels except at the bottom;
         # no y-axis labels except in the leftmost plots
         if (suppress.xlab == 'all' ||
             (suppress.xlab == 'margin' && matchidx$row != rows
              && plot.idx + cols <= num.plots)) { # last complete row
            plot.current <- plot.current + xlab(NULL)
         }
         if (suppress.ylab == 'all' ||
             (suppress.ylab == 'margin' && matchidx$col != 1)) {
            plot.current <- plot.current + ylab(NULL)
         }

         # do not fail if a single subplot isn't well-defined
         tryCatch(
            print(plot.current, vp = grid::viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col)),
            error = function(e) print(gplt.blank(e), vp = grid::viewport(layout.pos.row = matchidx$row,
                                                                layout.pos.col = matchidx$col)))
      }
   }
}

# auxiliary class to achieve consistent behavior of plotluck.multi with ggplot
# and plotluck: draw the plot if an only if the return value is not assigned.
#'@export
print.plotluck_multi <- function(x, ...) {
   if (x$in.grid) {
       mplot(x$plots, rows=x$rows, cols=x$cols,
             suppress.xlab=x$suppress.xlab, suppress.ylab=x$suppress.ylab)
   } else {
      lapply(x$plots, print)
   }
}

