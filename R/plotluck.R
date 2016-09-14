#'@import ggplot2
#'@importFrom plyr ddply summarise .
#'@importFrom stats as.formula median quantile runif xtabs


safe_log <- function(x,...) {
   ifelse(x<=0,0,log(x,...))
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
   y
}


decode.int.power <- function(x, base=10, offset=0) {
   if (offset == 0) {
      # regular log transform
      base^x
   } else {
      # log-modulus transform
      ifelse(x <= 0, - base^abs(x),
             ifelse(x == 1, 0,
                    base^abs(x-2)))
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

      sapply(seq(min.pos, max.pos, by),
             function(x) decode.int.power(x,base=base, offset=offset))
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


# calculate axis scale transform based on value distribution
# either:
# - identity
# - log-10
# - log-10-modulus: sign(x) * log(1+|x|)
# John & Draper 1980; see http://blogs.sas.com/content/iml/2014/07/14/log-transformation-of-pos-neg/

get.trans.fun <- function(data, x, expansion.thresh=2, label=x) {

   if (!is.numeric(data[[x]]) ||
       length(unique(data[[x]])) <= 2) {
      return(list(scales::identity_trans(), label))
   }

   q <- quantile(data[[x]], c(0, 0.25, 0.5, 0.75, 1), na.rm=TRUE)

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
      list(fun.trans,  sprintf('%s (log scale)', label))
   } else {
      list(scales::identity_trans(), label)
   }
}


Mode <- function(x) {
   ux <- unique(x)
   # if not returning a data frame, ddply will create an integer result column
   data.frame(ux[which.max(tabulate(match(x, ux)))])
}

weighted.mode <- function(data, x, w) {
   wtd <- plyr::ddply(data, x, function(d) sum(d[[w]], na.rm=TRUE))
   idx <- which.max(wtd[[2]])
   # need to rename to make later merge possible
   names(wtd)[1] <- '.wtd.mode.'
   wtd[idx, 1, drop=FALSE]
}

# calculate measure of central tendency per group.
# for numeric data[[x]], apply (weighted or unweighted) method as specified.
# for ordered factors, compute the (weighted or unweighted) mean or median over
#    the sorted integer levels, then round at the end.
# for unordered factors, always compute the mode (most frequent value)

group.central <- function(data, x, group.vars, w='NULL', method=c('median', 'mean'), col.name='.center.') {

   method <- match.arg(method)

   is.num <- is.numeric(data[[x]])
   is.ord <- (!is.num) && is.ordered(data[[x]])

   lev <- NULL
   grp.center <- NULL

   if (is.ord) {
      # treat ordinals as integer
      lev <- levels(data[[x]])
      data[[x]] <- as.integer(data[[x]])
   }

   if (w == 'NULL') {
      if (is.num || is.ord) {
         fun <- ifelse(method=='median', median, mean)
         grp.center <- plyr::ddply(data, group.vars, function(d) fun(d[[x]], na.rm=TRUE))
      } else {
         grp.center <- plyr::ddply(data, group.vars, function(d) Mode(d[[x]]))
      }
   } else {
      if (is.num || is.ord) {
         if (method == 'median') {
            # workaround for problems in wtd.quantile:
            # - if all values are NA, wtd.mean throws error: 'zero non-NA points'
            # - if weights are integer, sum can lead to NA due to overflow
            f <- function(d) if (all(is.na(d[[x]]))) { NA } else { Hmisc::wtd.quantile(d[[x]],
                                                                                       weights=as.numeric(d[[w]]), probs=0.5,
                                                                                       type='i/(n+1)', na.rm=TRUE, normwt=TRUE) }
            grp.center <- plyr::ddply(data, group.vars, f)
         } else {
            grp.center <- plyr::ddply(data, group.vars, function(d) Hmisc::wtd.mean(d[[x]], d[[w]], na.rm=TRUE))
         }
      } else { # !(is.num || is.ord)
         grp.center <- plyr::ddply(data, group.vars, function(d) weighted.mode(d, x, w))
      }
   }

   # for ordinals, reattach levels
   if (is.ord) {
      grp.center[[length(grp.center)]] <- ordered(lev[round(grp.center[[length(grp.center)]])], exclude=NULL)
   }

   names(grp.center)[length(grp.center)] <- col.name
   grp.center
}

# make all values of data[[x]] identically equal to the group mean
replace.by.central <- function(data, x, g, w) {
   x.avg <- group.central(data, x, g, w, 'mean', '.replace.center.')
   data <- merge(data, x.avg)
   data[[x]] <- data[[length(data)]]
   data[[length(data)]] <- NULL
   data
}

# return data frame with the levels of unordered factor x resorted by (weighted) frequency
order.factor.by.freq <- function(data, x, w='NULL', decreasing=FALSE) {
   x.data <- data[[x]]
   if (is.numeric(x.data) || is.ordered(x.data)) {
       # do not change ordered factors
      return(data)
   }

   if (w=='NULL') {
      tab <- table(x.data, useNA='ifany')
   } else {
      tab <- xtabs(as.formula(sprintf('%s ~ %s', w, x)), data, exclude=NULL, na.action=na.pass)
   }
   ord <- order(tab, decreasing=decreasing)
   data[[x]] <- factor(x.data, levels=levels(x.data)[ord],
                       exclude=NULL)
   return(data)
}


# return data frame with the levels of unordered factor x resorted by (weighted) central tendency of dependent variable y
order.factor.by.value <- function(data, x, y, w='NULL', decreasing=FALSE) {

   x.data <- data[[x]]
   if (is.numeric(x.data) || is.ordered(x.data)) {
      # do not change ordered factors
      return(data)
   }

   if (is.null(y) || y == 'NULL' || y == x) {
      return(order.factor.by.freq(data, x, w, decreasing))
   }

   y.data <- data[[y]]

   # we are displaying median lines for distributions, so the sorting should agree with this
   centers <- group.central(data, y, x, w, method='median')

   # if medians are equal for two or more classes, resolve ties by using the mean
   if (length(unique(centers$.center.)) != nrow(centers)) {
      centers <- group.central(data, y, x, w, method='mean')
   }
   if (!is.numeric(centers$.center.)) {
      centers$.center. <- as.integer(centers$.center.)
   }
   ord <- order(centers$.center., decreasing=decreasing)

   data[[x]] <- factor(x.data,
                       levels= centers[ord, x],
                       exclude=NULL)
   data
}

# reorder unordered factor levels of x,y,z to highlight dependencies

order.factors <- function(data, x, y, z='NULL', w='NULL') {

   # if at least one factor is ordered, sort the other factors accordingly
   # if none is ordered, sort according to z (if it exists), otherwise y

   if (z == 'NULL') {
      v.all <- c(y,x)
   } else {
      v.all <- c(z,y,x)
   }
   v.sort <- NULL
   for (v in v.all) {
      if (is.numeric(data[[v]]) || is.ordered(data[[v]])) {
         v.sort <- v
         break
      }
   }
   if (is.null(v.sort)) {
      v.sort <- v.all[1]
      data <- order.factor.by.freq(data, v.sort, w)
   }

   for (v in setdiff(v.all, v.sort)) {
      if ((!is.numeric(data[[v]])) && (!is.ordered(data[[v]]))) {
         # v is sortable
         data <- order.factor.by.value(data, v, v.sort, w)
      }
   }
   data
}

# limit the number of levels of a factor so as to not overcrowd the display
# - for unordered factors, keep the highest ones, group lower levels into '.other.'
# - for ordered factors, form about equidistant subgroups, in order
# - If there a slightly more levels than desired (according to tol), still don't
#   keep it, as quantizing makes the factor less intelligible
# - by.frequency keep the most frequent levels (according to count or weight)
#
limit.factor.levels <- function(data, x, w='NULL',
                                 max.levels=30,
                                 tol=max.levels*1.5,
                                 by.frequency=TRUE,
                                 reverse=FALSE) {
   x.data <- data[[x]]

   if (is.numeric(x.data)) {
      return(data)
   }

   u <- length(levels(x.data))

   if (u <= tol) {
      return(data)
   }
   warning(sprintf('Factor variable %s has too many levels (%d), truncating to %d', x, u, max.levels))

   if (!is.ordered(x.data)) {
      ord <- seq(length(levels(x.data)))
      if(by.frequency) {
         # keep the most frequent levels in the same order, combine the rest into '.other.',
         data.tmp <- data
         data.tmp[[x]] <- as.integer(data.tmp[[x]])
         tab <- NULL
         if (w=='NULL') {
            tab <- table(data.tmp[[x]], useNA='ifany')
         } else {
            tab <- xtabs(as.formula(sprintf('%s ~ %s', w, x)), data.tmp,  exclude=NULL, na.action=na.pass)
         }
         # make sure that ties are broken according to previous order!
         ord <- order(tab, as.numeric(names(tab)), decreasing=FALSE)
      }
      levels.ord <- levels(x.data)[ord]
      if (!reverse) {
         levels.reduced <- levels.ord[seq(u-max.levels+1,u)]
      } else {
         levels.reduced <- levels.ord[seq(max.levels)]
      }
      # note: we want to maintain the original order of levels!
      idx.trunc <- !(x.data %in% levels.reduced)
      x.data <- factor(x.data, levels=c('.other.', levels(x.data)),exclude=NULL)
      x.data[idx.trunc] <- '.other.'
      data[[x]] <- factor(x.data, exclude=NULL)

   } else {  # is.ordered(data[[x]])

     #  Note: The weighted case is very hard to get right for very unequal weights.
     #  Ignoring for now.

      breaks <- 1 + (seq(0,max.levels)/max.levels) * (length(levels(x.data))-1)

      # dig.lab: in cut(), avoids scientific formatting for integers
      dig.lab=50
      x.data <- cut(as.integer(x.data), breaks=breaks,
                    right=FALSE, include.lowest=TRUE, ordered.result=TRUE, dig.lab=dig.lab)

      # recreate level information by parsing those generated by 'cut'
      # make lists of elements, e.g.:
      #  - original levels 'a', 'b', 'c', 'd'
      #  - binned level '[2,4)' ->  'c,d'
      x.data <- as.factor(x.data)

      reconstruct.levels <- function(interval, levels) {
         left.bracket <- substr(interval,1,1)
         right.bracket <- substr(interval,nchar(interval),nchar(interval))
         eps <- 0
         if (left.bracket == '(') {
            eps <- 1e-10
         }
         if (right.bracket == ')') {
            eps <- -1e-10
         }
         left.idx  <- ceiling(eps + as.numeric(gsub('[(\\[]([0-9.]+),.*', '\\1', interval)))
         right.idx <- floor(eps + as.numeric(gsub('.*,([0-9.]+).', '\\1', interval)))

         # catch if level doesn't fit interval pattern - this shouldn't happen!
         if (is.na(left.idx) || is.na(right.idx)) {
            return(interval);
         }
         return(paste(levels[left.idx:right.idx], sep='', collapse=','))
      }

      levels(x.data) <- sapply(levels(x.data), reconstruct.levels, levels=levels(data[[x]]))
      x.data <- x.data[,drop=T]
      data[[x]] <- x.data
   }
   data
}

# if one numerical axis has very few values, it sometimes makes sense to treat them as a factor
discretize.few.unique <- function(data, x, few.unique.as.factor=5) {
   if (!is.numeric(data[[x]])) {
      return(data)
   }
   non.na <- !is.na(data[[x]])
   u <- length(unique(data[non.na,x]))

   if (u <= few.unique.as.factor &&
       nrow(data) > u * u) { # heuristic similar to histogram binning
      data[[x]] <- ordered(data[[x]], exclude=NULL)
   }
   data
}


# determine number of histogram bins
# note: this is a modification the 'nclass.FD' function, which sometimes returns '1'
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


# break a numeric variable into intervals
# - method=quantile - bins with equal number of samples
# - method=histogram - equidistant intervals
# - max.breaks - requested number of breaks
# - tol - if there are less than that many distinct values to
#     begin with, don't do additional binning
# - estimate.breaks=FALSE - use exactly max.breaks many breaks

discretize <- function(data, x,
                       w='NULL',
                       max.breaks,
                       tol=max.breaks*1.5,
                       estimate.breaks=TRUE,
                       method=c('quantile', 'histogram')) {
   x.data <- data[[x]]

   if (!is.numeric(x.data)) {
      return(data)
   }

   u <- length(unique(x.data))

   if (missing(max.breaks)) {
      max.breaks = ceiling(min(u-1, sqrt(nrow(data))))
   }
   max.breaks = floor(min(max.breaks, u-1, sqrt(nrow(data))))
   max.breaks <- min(max.breaks, nrow(data))


   if (u - 1 <= max(tol, max.breaks)) {
      return(data)
   }

   method <- match.arg(method)

   breaks <- NULL
   # compute number of breaks
   if (method == 'quantile') {
      probs <- seq(0, max.breaks) / max.breaks
      if (w=='NULL') {
         breaks <- quantile(x.data, probs=probs, na.rm=TRUE)
      } else {
         breaks <- Hmisc::wtd.quantile(x.data, weights=data[[w]],
                                       probs=probs, na.rm=TRUE, normwt=TRUE,
                                       type='i/(n+1)')
      }
      data[[x]] <- ordered(cut(x.data, breaks=unique(breaks), include.lowest=TRUE),
                           exclude=NULL)

   } else { # method == 'histogram'
      n.breaks <- max.breaks
      if (estimate.breaks) {
         n.breaks <- nclass.FD.modified(x.data, max.breaks)
      }
      breaks <- graphics::hist(x.data, breaks=n.breaks, plot=FALSE)$breaks
      # leave type numeric
      step<-breaks[2] - breaks[1]
      data[[x]] <- round((x.data-breaks[1])/step)*step+breaks[1]
   }

   data
}

# - convert character vectors into factors
# - with na.rm=TRUE, remove rows containing NA
# - from here on, don't exclude any levels from factors any more
# - check for all-NA and constant values
# - convert character vectors into factors
# - convert logicals into ordered factors

preprocess.factors <- function(data, na.rm) {

   for (x in names(data)) {

      x.data <- data[[x]]
      non.na <- !is.na(x.data)
      if (length(which(non.na)) == 0) {
         stop(sprintf('Variable %s has only missing values, giving up', x))
      }

      u <- length(unique(x.data[non.na]))
      if (u == 1) {
         stop(sprintf('Variable %s has only single level ("%s"), giving up', x, unique(x.data[non.na])))
      }

      if (u == 2) {
         # binary variables are trivially ordered
         data[[x]] <- ordered(x.data)
      } else if (!is.numeric(x.data)) {

         if (!is.factor(x.data)) {
            # i.e., character or logical
            x.data <- factor(x.data, exclude=NULL)
         }

         # suppress unused levels
         # for NA treatment, see http://r.789695.n4.nabble.com/droplevels-inappropriate-change-td4723942.html
         # droplevels(x.data, exclude=exclude.factor) does not seem to be working correctly
         data[[x]] <- x.data[,drop=TRUE]
      }
   }

  if (na.rm) {
    factor.with.na.level <- any(sapply(seq(length(data)),
                                       function(x) {is.factor(data[[x]]) && any(is.na(levels(data[[x]])))}))
    if (factor.with.na.level) {
      warning('Option "na.rm=TRUE", but data contains factor(s) with an NA level. Keeping these ...')
    }
    data <- na.omit(data)

    if (nrow(data) == 0) {
      stop(sprintf('Variable(s) have only missing values, giving up'))
    }

  }

  # if na.rm==TRUE, all NA's have been removed. In either case, from now on all levels should be
  # included

  for (i in seq(length(data))) {
    if (is.factor(data[[i]])) {
      data[[i]] <- factor(data[[i]], exclude=NULL, ordered=is.ordered(data[[i]]))
      # note: rename <NA> - else ggplot doesn't show a color legend
      idx<-which(is.na(levels(data[[i]])))
      if (length(idx) > 0) {
        levels(data[[i]])[idx] <- '?'
      }
    }
  }

   data
}


marginalize <- function(data, var.list, w='NULL') {
   rhs <- paste(var.list, sep='', collapse='+')
   prop.table(xtabs(as.formula(sprintf('%s ~ %s', ifelse(w=='NULL', c(''), w), rhs)), data, exclude=NULL, na.action=na.pass))
}

# entropy: H(vars) = - \sum p(vars)log2(p(vars))
# assumption: vars are factors
entropy <- function(data, vars, w) {
   tab  <- marginalize(data, vars, w)
   log.tab <- log2(tab)
   log.tab[tab==0] <- 0 # guard for -Inf
   -sum(tab * log.tab)
}


# conditional entropy: H(y|x) = H(x,y) - H(x)
# assumption: x, y are factors
cond.entropy <- function(data, x, y, w) {
   if (x == y) {
      0
   } else {
      hxy <- entropy(data, c(x, y), w)
      hx <- entropy(data, x, w)
      hxy - hx
   }
}

quantize <- function(data, v, target='NULL', w, n.levels, estimate.breaks=FALSE, method='histogram') {
  if (is.numeric(data[[v]])) {
    data <- discretize(data, v, w=w,
                       max.breaks=n.levels, tol=n.levels, estimate.breaks, method=method)
  } else {
    if (target != 'NULL') {
      data <- order.factor.by.value(data, v, target, w=w)
    }
    data <- limit.factor.levels(data, v, w,
                                max.levels=n.levels,
                                tol=n.levels)
  }
  data
}


# return a list of conditional entropies H(target|x), for each column in the data frame.
# for comparability, all variables are quantized into the same number of bins
cond.entropy.data <- function(data, target, given, w='NULL') {
   n.row <- nrow(data)
   n.levels <- floor(max(2, min(log2(n.row), n.row/10))) # heuristic

   cond.entropy.quantized <- function(data, x, y, w, n.levels) {
      if (x == y) {
         0
      } else {
         x.data <- quantize(data, x, y, w, n.levels)
         # note: we need to make a copy for the special case of x=w
         names(x.data)[names(x.data)==x]<-'.x'
         if (w != 'NULL') {
            x.data[[w]] <- data[[w]]
         }
         cond.entropy(x.data, '.x', y, w)
      }
   }

   entropies <- matrix(nrow=length(data), ncol=length(data))
   colnames(entropies) <- names(data)
   rownames(entropies) <-names(data)

   if (!missing(target)) {
      # single row/column
      data <- quantize(data, target, w=w, n.levels=n.levels)
      entropies[target,] <- sapply(names(data), function(x)
         cond.entropy.quantized(data, x, target, w, n.levels))
   } else {
      for (target in names(data)) {
         data <- quantize(data, target, w=w, n.levels=n.levels)
         if (!missing(given)) {
            entropies[target, given] <- cond.entropy.quantized(data, x=given, y=target, w, n.levels)
         } else {
            # all vs all
            entropies[target,] <- sapply(names(data), function(x)
               cond.entropy.quantized(data, x, target, w, n.levels))
         }

      }
   }
   entropies
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
   u >= min.coverage * ux * uy
}

# some common themes
theme_panel_num_x <-
   theme(panel.grid.major.x=element_line(color='black', linetype='dotted'),
         panel.grid.minor.x=element_line(color='black',linetype='dotted'))
theme_panel_num_y <-
   theme(panel.grid.major.y=element_line(color='black', linetype='dotted'),
         panel.grid.minor.y=element_line(color='black',linetype='dotted'))
# points on the grid line are hard to see; remove grid line for factors
theme_panel_fac_x <- theme(panel.grid.major.x=element_blank(),
                           panel.grid.minor.x=element_blank())
theme_panel_fac_y <- theme(panel.grid.major.y=element_blank(),
                           panel.grid.minor.y=element_blank())
# for factors, write horizontal axis labels at an angle to avoid overlap
theme_slanted_text_x <- theme(axis.text.x=element_text(angle=-45, hjust=0, vjust=1))

# detect if x axis text might overlap
estimate.label.overlap <- function(labels, breaks=seq(length(labels))) {
   if (length(breaks) < 2) {
      return(FALSE)
   }
   rg <- breaks[length(breaks)] - breaks[1]
   for (i in seq(length(breaks)-1)) {
      space <- (breaks[i+1] - breaks[i]) / rg
      label.width <- as.numeric(grid::convertX(unit(1, 'strwidth', labels[i]),'npc'))
      if (label.width > space) {
         return(TRUE)
      }
   }
   return(FALSE)
}

# add color and/or fill scale layers to plot
# - use qualitative (sequential) brewer scale for (un)ordered factors, if
#   palette size supports it.
#   Otherwise, create a colorRampPalette for ordered factors, or use the default if unordered.
# - use a gradient fill for numeric vectors
# - if too few colors, use default (first two colors of brewer palette might not look good by themselves!)

add.color.fill <- function(p, data, x, aesth=c('color', 'fill'),
                           colors.gradient=c(hcl(66,60,85), hcl(128,100,35)),
                           palette.brewer.seq='YlGn',
                           palette.brewer.qual='Set1') {

   na.value <- 'darkgray'
   if ('color' %in% aesth) {
     p <- p + aes_string(color=x)
   }
   if ('fill' %in% aesth) {
     p <- p + aes_string(fill=x)
   }
   if (length(unique(data[[x]])) <= 2) {
      # just use default
   } else if (is.numeric(data[[x]])) {
      if ('color' %in% aesth) {
         p <- p + scale_color_gradientn(colors=colors.gradient)
      }
      if ('fill' %in% aesth) {
         p <- p + scale_fill_gradientn(colors=colors.gradient)
      }
   } else if (is.ordered(data[[x]])) {
      if (length(levels(data[[x]])) <= 8) {
         if ('color' %in% aesth) {
            p <- p + scale_color_brewer(palette=palette.brewer.seq, na.value=na.value)
         }
         if ('fill' %in% aesth) {
            p <- p + scale_fill_brewer(palette=palette.brewer.seq, na.value=na.value)
         }
      } else { # length(levels(data[[x]])) > 8
         pal <- colorRampPalette(colors.gradient)(length(levels(data[[x]])))
         names(pal) <- levels(data[[x]])
         if ('color' %in% aesth) {
            p <- p + scale_color_manual(values=pal, na.value=na.value)
         }
         if ('fill' %in% aesth) {
            p <- p + scale_fill_manual(values=pal, na.value=na.value)
         }
      }
   } else { # unordered factor
      if (length(levels(data[[x]])) <= 9) {
         if ('color' %in% aesth) {
            p <- p + scale_color_brewer(palette=palette.brewer.qual, na.value=na.value)
         }
         if ('fill' %in% aesth) {
            p <- p + scale_fill_brewer(palette=palette.brewer.qual, na.value=na.value)
         }
      } else { # length(levels(data[[x]])) > 9
         # use default
         if ('color' %in% aesth) {
            p <- p + scale_color_hue()
         }
         if ('fill' %in% aesth) {
            p <- p + scale_fill_hue()
         }
      }
   }
   p
}

# place legend either at the bottom or on the right, wherever it occupies less space
# pack as tightly as possible
add.color.legend <- function(p, data, x, aesth=c('color','fill')) {

   aesth=match.arg(aesth)
   p <- p + theme(legend.key.width=unit(0.5,'lines'))
   if (is.numeric(data[[x]])) {
      # color bars always on the right margin, to avoid number label overwriting
      if (aesth == 'color') {
         p <- p + guides(color=guide_colorbar(direction='vertical'), fill=FALSE)
      } else {
         p <- p + guides(fill=guide_colorbar(direction='vertical'), color=FALSE)
      }
      p <-  p + theme(legend.position='right', legend.background=element_blank())
      return(p)
   }

   # rough estimate of key symbol plus spacing; in the default theme.gray(), net key size is 1.2 lines
   size.key <- as.numeric(grid::convertX(unit(1.6, 'lines'), 'npc'))
   title.length <- as.numeric(grid::convertX(unit(1, 'strwidth', x), 'npc'))
   legend.margin <- as.numeric(unit(0.2, 'cm'))

   # vertical placement
   ll <- length(levels(data[[x]]))
   npc.length <- sapply(levels(data[[x]]), function(x) grid::convertX(unit(1, 'strwidth', x), 'npc'))
   npc.length <- sort(npc.length, decreasing=TRUE)
   ncol.vert <- ceiling(ll*size.key)
   nrow.vert <- ceiling(ll/ncol.vert)
   npc.width.vert <- legend.margin + max(title.length, ncol.vert * (size.key + npc.length[1]))

   # horizontal placement
   # try out how many legend columns we can fit underneath the graph
   ncol.horz <- 0
   npc.horz <- title.length
   while (ncol.horz < ll && (npc.horz + npc.length[ncol.horz+1] < 1)) {
      ncol.horz <- ncol.horz + 1
      npc.horz <- npc.horz + npc.length[ncol.horz] + size.key
   }
   nrow.horz <- ceiling(ll/ncol.horz)
   npc.height.horz <- legend.margin + nrow.horz * size.key
   aspect.ratio <- as.numeric(grid::convertX(unit(1, 'npc'),'cm')) /
      as.numeric(grid::convertY(unit(1, 'npc'),'cm'))
   if (npc.width.vert < aspect.ratio * npc.height.horz) {
      # vertical - by column, starting with first level at the bottom
      if (aesth == 'color') {
         p <- p + guides(color=guide_legend(direction='vertical', byrow=FALSE,
                                            reverse=TRUE, ncol=ncol.vert), fill=FALSE)
      } else {
         p <- p + guides(fill=guide_legend(direction='vertical', byrow=FALSE,
                                           reverse=TRUE, ncol=ncol.vert), color=FALSE)
      }
      p + theme(legend.position='right', legend.background=element_blank())
   } else {
      # horizontal - by row
      if (aesth == 'color') {
         p <- p + guides(color=guide_legend(direction='horizontal', byrow=TRUE,
                                       reverse=FALSE, nrow=nrow.horz), fill=FALSE)
      } else {
         p <- p + guides(fill=guide_legend(direction='horizontal', byrow=TRUE,
                                      reverse=FALSE, nrow=nrow.horz), color=FALSE)
      }
      p + theme(legend.position='bottom', legend.background=element_blank())
   }
}

add.smooth.line <- function(p, w='NULL', fill.smooth='NULL') {

   if (w == 'NULL') {
      # note: for more than 1000 points, ggplot2 uses mgcv package,
      # leads to weird dependency issue
      if (fill.smooth == 'NULL') {
         p + geom_smooth(alpha=0.2, na.rm=TRUE)
      } else {
         p + geom_smooth(alpha=0.2, na.rm=TRUE, color=fill.smooth, fill=fill.smooth)
      }
   } else {
      if (fill.smooth == 'NULL') {
         p + geom_smooth(aes_string(weight=w), alpha=0.2, na.rm=TRUE)
      } else {
         p + geom_smooth(aes_string(weight=w), alpha=0.2, na.rm=TRUE, color=fill.smooth, fill=fill.smooth)
      }
   }
}

add.axis.transform <- function(p, data, x, ax=c('x','y'), trans.log.thresh=2) {
   if (!is.numeric(data[[x]])) {
      p
   } else {
      ax <- match.arg(ax)
      trans <- get.trans.fun(data, x, trans.log.thresh)
      if (ax == 'x') {
         # note: "scale_x_continuous(trans=trans.x[[1]], name=trans.x[[2]])" doesn't allow
         # to overwrite the label later
         p + scale_x_continuous(trans=trans[[1]]) + xlab(trans[[2]])
      } else {
         p + scale_y_continuous(trans=trans[[1]]) + ylab(trans[[2]])
      }
   }
}

# 2D scatter plot of numeric or factor variables
# - overlay smoothing line if both variables are numeric
# - scatter.dedupe='area' - unique repeated points, count frequency
# - scatter.dedupe='jitter' - plot each repeated point separately, add
#   jitter if there are more than min.points.jitter with identical coordinates
# - counts/weights are drawn with a shaded circle of proportional area
# - jitter.x, jitter.y - the amount of jittering as a multiple of resolution
# - trans.log.thresh - threshold to decide on log-transform
# - fill.smooth - fill color for smoothing line
gplt.scatter <- function(data, x, y, w='NULL',
                         scatter.dedupe=c('area','jitter'),
                         min.points.jitter=3,
                         jitter.x=0.4,
                         jitter.y=0.4,
                         trans.log.thresh=2,
                         max.factor.levels=30,
                         fill.smooth='NULL',
                         ...) {

   scatter.dedupe <- match.arg(scatter.dedupe)
   flip <- FALSE
   if (is.numeric(data[[x]]) && (!is.numeric(data[[y]]))) {
      # HACK: ggplot2 does not implement vertical dodging, therefore we
      # use coord_flip() as a workaround
      flip <- TRUE
      tmp <- x
      x <- y
      y <- tmp
   }

   num.x <- is.numeric(data[[x]])
   num.y <- is.numeric(data[[y]])
   ord.x <- is.ordered(data[[x]])
   ord.y <- is.ordered(data[[y]])

   data <- order.factors(data, x, y, w=w)
   for (v in c(x,y)) {
      data <- limit.factor.levels(data, v, w=w,
                                  max.levels=max.factor.levels)
   }

   if (scatter.dedupe == 'area') {
      # record frequency/total weight for each unique row of the data
      if (w == 'NULL') {
         # equivalent to 'uniq -c'
         data <- plyr::ddply(data, names(data), function(D) nrow(D))
         names(data)[length(data)] <- '.freq.'
         w <- '.freq.'
      } else {
         data <- plyr::ddply(data, setdiff(names(data), w), function(D) sum(D[[w]], na.rm=TRUE))
         names(data)[length(data)] <- w
      }
   }

   # dodging

   if (!num.x) {
      # example involving dodging:
      # i2<-iris
      # i2$c<-factor(cut(i2$Petal.Width,3))
      # names(i2)[6]<-'pw.quant'
      # plotluck(Petal.width~Species|pw.quant,i2, opts=plotluck.options(geom='scatter'))
      pos <- position_dodge(0.5)
   } else {
      pos <- position_identity()
   }

   # if there are a lot of points, better make them transparent
   alpha <- max(0.3, 0.8 - 0.8/2000 * nrow(data))

   p <- ggplot(data, aes_string(x=x, y=y, weight=w))

   if (w != 'NULL' && length(unique(data[[w]])) > 1) {
      # use point size to represent count/weight
      p <- p + geom_point(aes_string(size=w), alpha=alpha,
                          position=pos, na.rm=TRUE, ...) +
         scale_size(guide=FALSE)
   } else {
      # use jittering
      if (max(table(data[, c(x, y)], useNA='ifany')) >= min.points.jitter) {
         # test for repeated points, optionally jitter
         # - if both x and y are numeric, jitter in both directions
         # - if one is a factor, only jitter in this direction
         #   (enough empty space between levels)
        if (num.x && !num.y) {
          jitter.x <- 0
          # resolution(data[[y]], zero=FALSE) == 1
        } else if (!num.x && num.y) {
          # resolution(data[[x]], zero=FALSE) == 1
          jitter.y <- 0
        } else {
          jitter.x <- jitter.x * resolution(data[[x]], zero=FALSE)
          jitter.y <- jitter.y * resolution(data[[y]], zero=FALSE)
        }
      } else {
         jitter.x <- 0
         jitter.y <- 0
      }

      p <- p + geom_point(alpha=alpha, position=position_jitter(width=jitter.x, height=jitter.y),
                          na.rm=TRUE, ...)
   }

   # draw smoothing line
   if (num.x && num.y) {
      p <- add.smooth.line(p, w, fill.smooth)
   }

   # axis transformation
   p <- add.axis.transform(p, data, x, 'x', trans.log.thresh)
   p <- add.axis.transform(p, data, y, 'y', trans.log.thresh)

   if (!flip) {
      if (num.x) {
         p <- p + theme_panel_num_x + theme_panel_num_y
      } else {
         p <- p + theme_panel_num_y
         if (length(levels(data[[x]])) < 10) {
            # grid lines for factors only necessary when very many
            p <- p + theme_panel_fac_x
         } else {
            p <- p + theme(panel.grid.major.x=element_line(color='black', linetype='dotted'),
                           panel.grid.minor.x=element_blank())
         }
         if (estimate.label.overlap(levels(data[[x]]))) {
            p <- p + theme_slanted_text_x
         }
      }
   } else {
      p <- p + coord_flip()
      if (num.x) {
         p <- p + theme_panel_num_y + theme_panel_num_x
      } else {
         p <- p + theme_panel_num_x
         if (length(levels(data[[x]])) < 10) {
            p <- p + theme_panel_fac_y
         } else {
            p <- p + theme(panel.grid.major.y=element_line(color='black', linetype='dotted'),
                           panel.grid.minor.y=element_blank())
         }
      }
   }
   p
}


# hexbin plot with overlayed smoothing line
gplt.hex <- function(data, x, y, w='NULL', trans.log.thresh=2,
                     fill.smooth='NULL', ...) {

   p <- ggplot(data, aes_string(x=x, y=y)) +
      geom_hex(na.rm=TRUE, ...)

   p <- add.smooth.line(p, w, fill.smooth)

   # axis transformation
   p <- add.axis.transform(p, data, x, 'x', trans.log.thresh)
   p <- add.axis.transform(p, data, y, 'y', trans.log.thresh)

   p +
      # log scaling of color often reveals more details
      scale_fill_gradientn(colors=c(hcl(66,60,95), hcl(128,100,35)), trans='log', guide=FALSE) +
      theme(legend.position='right') +
      theme_panel_num_x + theme_panel_num_y
}


# raster plot
# point grid of discretized values, with optional discretized color (z)
# the color is a representation of an appropriate central tendency measure
# for the point bin (mode for factors, median for ordinal and numeric vectors)
gplt.raster <- function(data, x, y, z='NULL', w='NULL', trans.log.thresh=2,
                        max.factor.levels=30,
                        raster.resolution=30,
                        colors.gradient=c(hcl(66,60,85), hcl(128,100,35)),
                        palette.brewer.seq='YlGn',
                        palette.brewer.qual='Set1',
                        ...) {

   num.x <- is.numeric(data[[x]])
   ord.x <- is.ordered(data[[x]])
   num.y <- is.numeric(data[[y]])
   ord.y <- is.ordered(data[[y]])

   num.z <- FALSE
   ord.z <- FALSE
   ex.z <- (!is.null(z)) && (z != 'NULL')
   if (ex.z) {
      num.z <- is.numeric(data[[z]])
      ord.z <- is.ordered(data[[z]])
      if (!num.z && !ord.z) {
         data <- order.factor.by.freq(data, z, w)
      }
   }

   data <- order.factors(data, x, y, z, w)

   # quantize variables
   for (v in c(x,y)) {
      if (is.numeric(data[[v]])) {
         data <- discretize(data, v, w=w, max.breaks=raster.resolution, estimate.breaks=FALSE,
                            method='histogram')
      } else {
         data <- limit.factor.levels(data, v, w=w,
                                     max.levels=raster.resolution, tol=max.factor.levels)
      }
   }

   g.vars <- names(data)
   if (ex.z) {
      g.vars <- setdiff(g.vars, z)
   }

   if (w != 'NULL') {
      g.vars <- setdiff(g.vars, w)
   }

   # HACK: remove precomputed columns
   g.vars <- setdiff(g.vars, '.center.')

   # replace z values with group center
   if (ex.z) {
      data <- replace.by.central(data, z, g.vars, w)

      if (num.z) {
         # if z distribution is very skewed, better to quantize than to histogram
         trans.z <- get.trans.fun(data, z, trans.log.thresh)
         method <- 'histogram'
         if (trans.z[[1]]$name != 'identity') {
            method <= 'quantile'
         }

         # limited by color palette
         data <- discretize(data, z, w=w, max.breaks=9, tol=9,
                            estimate.breaks=FALSE, method=method)
      } else {
         data <- limit.factor.levels(data, z, w,
                                      max.levels=9,
                                      tol=9)
      }
      g.vars <- c(g.vars, z)
   }

   # sum up counts or weights
   if (w == 'NULL') {
      data <- plyr::ddply(data, g.vars, function(D) nrow(D))
      names(data)[length(data)] <- '.freq.'
      w <- '.freq.'
   } else {
      data <- plyr::ddply(data, g.vars, function(D) sum(D[[w]], na.rm=TRUE))
      names(data)[length(data)] <- w
   }

   p <- ggplot(data, aes_string(x=x, y=y, size=w, color=z, fill=z, weight=w)) +
      geom_point(na.rm=TRUE, ...) +
      scale_size(guide=FALSE)

   # axis transformation
   p <- add.axis.transform(p, data, x, 'x', trans.log.thresh)
   p <- add.axis.transform(p, data, y, 'y', trans.log.thresh)

   if (ex.z) {
      p <- add.color.fill(p, data, z,
                          colors.gradient=colors.gradient,
                          palette.brewer.seq=palette.brewer.seq,
                          palette.brewer.qual=palette.brewer.qual)
      p <- add.color.legend(p, data, z, 'color')
   }

   # show gridlines
   p <- p + theme_panel_num_x + theme_panel_num_y
   if (!num.x && estimate.label.overlap(levels(data[[x]]))) {
      p <- p + theme_slanted_text_x
   }

   p
}

# cleveland dot plot or bar plot with stat_bin
# no log transforms for bars; see http://www.perceptualedge.com/articles/b-eye/dot_plots.pdf
gplt.dot <- function(data, x, w='NULL', vertical=TRUE,
                     max.factor.levels=30, geom=c('dot', 'bar'), ...) {

   geom <- match.arg(geom)

   if (!is.ordered(data[[x]])) {
      data <- order.factor.by.freq(data, x, w)
   }

   data <- limit.factor.levels(data, x, w=w,
                               max.levels=max.factor.levels)

   ylab <- 'count'
   if (w != 'NULL') {
      ylab <- w
   }

   p <- ggplot(data, aes_string(x=x, weight=w))
   if (geom == 'dot') {
      # if there are only few points, plot direction is not obvious to the viewer
      if (!vertical) {
         shp <- '^'
      } else {
         shp <- '>'
      }
      p <- p + geom_point(stat='count', shape=shp, size=6, na.rm=TRUE, ...)
   } else {
      p <- p + geom_bar(position=position_dodge(), alpha=0.8, width=0.2, na.rm=TRUE, ...)
   }

   if (!vertical) {
      p <- p + theme_panel_num_y + ylab(ylab)

      if (length(levels(data[[x]])) < 10) {
         p <- p + theme_panel_fac_x
      } else {
         p <- p + theme(panel.grid.major.x=element_line(color='black', linetype='dotted'),
                        panel.grid.minor.x=element_blank())
      }

      if (estimate.label.overlap(levels(data[[x]]))) {
         p <- p + theme_slanted_text_x
      }

   } else { # horizontal
      p <- p + coord_flip() + theme_panel_num_x + ylab(ylab)
      if (length(levels(data[[x]])) < 10) {
         p <- p + theme_panel_fac_y
      } else {
         p <- p + theme(panel.grid.major.y=element_line(color='black', linetype='dotted'),
                        panel.grid.minor.y=element_blank())
      }
   }
   p
}


# violin plot with overlayed median point
# assumption: x is factor, y is numeric
gplt.violin <- function(data, x, y, w='NULL',
                        trans.log.thresh=2,
                        max.factor.levels=30,
                        ...) {

   flip <- FALSE
   if (is.numeric(data[[x]])) {
      if (is.numeric(data[[y]])) {
         stop('violin plot requires one factor variabe')
      }
      flip <- TRUE
      tmp <- x
      x <- y
      y <- tmp
   }
   # note: geom_violin throws an error if all values are equal
   # workaround: add very small jitter
   data[[y]] <- data[[y]] +
      1E-8 * min(data[[y]], na.rm=TRUE) * (runif(nrow(data)) - 0.5)

   if (!is.ordered(data[[x]])) {
      data <- order.factor.by.value(data, x, y, w)
   }

   data <- limit.factor.levels(data, x, w=w,
                               max.levels=max.factor.levels)

   p <- ggplot(data, aes_string(x=x, y=y, weight=w, ymax=max(y,na.rm=TRUE))) +
      geom_violin(scale='width', alpha=0.8, na.rm=TRUE, ...)

   if ('.center.' %in% names(data)) {
      # dodge does not work correctly when width is not specified
      # see https://github.com/hadley/ggplot2/issues/525
      p <- p + geom_point(mapping=aes(y=.center.),
                          position=position_dodge(width=0.9),
                          size=2, shape=1, na.rm=TRUE)
   }

   # axis transformation
   p <- add.axis.transform(p, data, y, 'y', trans.log.thresh)

   if (!flip) {
      p <- p + theme_panel_fac_x + theme_panel_num_y

      if (estimate.label.overlap(levels(data[[x]]))) {
         p <- p + theme_slanted_text_x
      }

   } else {
      p <- p + coord_flip() + theme_panel_fac_y + theme_panel_num_x
   }
   p
}


# box plot
# assumption: x is factor, y is numeric
gplt.box <- function(data, x, y, w='NULL',
                     trans.log.thresh=2,
                     max.factor.levels=30,
                     ...) {

   flip <- FALSE
   if (is.numeric(data[[x]])) {
      if (is.numeric(data[[y]])) {
         stop('box plot requires one factor variabe')
      }
      flip <- TRUE
      tmp <- x
      x <- y
      y <- tmp
   }

   if (!is.ordered(data[[x]])) {
      data <- order.factor.by.value(data, x, y, w)
   }

   data <- limit.factor.levels(data, x, w=w,
                               max.levels=max.factor.levels)

   p <- ggplot(data, aes_string(x=x, y=y, weight=w, ymax=max(y,na.rm=TRUE))) +
      geom_boxplot(alpha=0.8, na.rm=TRUE, ...)

   # axis transformation
   p <- add.axis.transform(p, data, y, 'y', trans.log.thresh)

   if (!flip) {
      p <- p + theme_panel_fac_x + theme_panel_num_y

      if (estimate.label.overlap(levels(data[[x]]))) {
         p <- p + theme_slanted_text_x
      }
   } else {
      p <- p + coord_flip() + theme_panel_fac_y + theme_panel_num_x
   }
   p
}


# density plot with overlayed vertical median lines
gplt.density <- function(data, x, w='NULL',
                         trans.log.thresh=2,
                         ...) {

   p <- ggplot(data, aes_string(x=x, weight=w)) +
      geom_density(aes(y=..scaled..), alpha=0.6, adjust=0.5, trim=TRUE, na.rm=TRUE, ...) +
      geom_rug(na.rm=TRUE)

   # geom_vline cannot be used, see https://github.com/hadley/ggplot2/issues/426
   # .center. explicitly precomputed before
   # known issue: color cannot be set here, dependent on conditional layer
   # (does not automatically inherit)
   if ('.center.' %in% names(data)) {
      p <- p + geom_vline(aes(xintercept=.center.), linetype='dashed', size=0.7, na.rm=TRUE)
   }

   # axis transformation
   p <- add.axis.transform(p, data, x, 'x', trans.log.thresh)

   ylab <- 'density'
   if (w != 'NULL') {
      ylab <- sprintf('%s density', w)
   }

   p + ylab(ylab) +
      # density numbers themselves not very meaningful
      theme(axis.ticks.y=element_blank(),
            axis.text.y=element_blank()) +
      theme_panel_num_x
}

# histogram plot with overlayed vertical median lines
gplt.histogram <- function(data, x, w='NULL',
                           n.breaks=NA,
                           trans.log.thresh=2,
                           ...) {

   if (is.numeric(n.breaks)) {
      n.breaks <- n.breaks
   } else {
      n.breaks <- nclass.FD.modified(data[[x]])
   }
   bin.width <- diff(range(data[[x]], na.rm=TRUE))/n.breaks
   p <- ggplot(data, aes_string(x=x, weight=w)) +
      geom_histogram(binwidth=bin.width, position='identity', alpha=0.6, na.rm=TRUE, ...) +
      geom_rug(na.rm=TRUE)

   if ('.center.' %in% names(data)) {
      p <- p + geom_vline(aes(xintercept=.center.), linetype='dashed', size=0.7, na.rm=TRUE)
   }

   # axis transformation
   p <- add.axis.transform(p, data, x, 'x', trans.log.thresh)

   ylab <- 'count'
   if (w != 'NULL') {
      ylab <- w
   }

   p + ylab(ylab) +
      theme_panel_num_x + theme_panel_num_y
}


# spine plot of two variables
# note: We could use mosaicplot, but we want to consistently stick with ggplot
#    another option would be the prodplots package, but currently seems
#    to be in development stage
#
# - plot.margin.x - horizontal gap between rectangles for x-values
# - no gaps between y-rectangles
#
# note: additional columns other than x,y,w are dropped; this precludes later
# faceting on them. They could be preserved, however this function is already
# complicated as it is.

gplt.spine <- function(data, x, y, w='NULL',
                       plot.margin.x=0.05,
                       max.factor.levels=10,
                       colors.gradient=c(hcl(66,60,85), hcl(128,100,35)),
                       palette.brewer.seq='YlGn',
                       palette.brewer.qual='Set1',
                       ...) {

   # ordering and factor truncation
   data <- order.factors(data, x, y, w=w)
   for (v in c(x, y)) {
      data <- quantize(data, v, w=w, n.levels=max.factor.levels)
   }

   x.lev <- length(levels(data[[x]]))
   y.lev <- length(levels(data[[y]]))

   xy.tab  <- marginalize(data, c(x, y), w)
   x.tab   <- marginalize(data, x, w)
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

   # remove empty rects
   idx <- plot.data$bottom < plot.data$top & plot.data$left < plot.data$right
   plot.data <- plot.data[idx,]

   p <- ggplot(plot.data, aes(xmin=left, xmax=right, ymin=bottom, ymax=top)) +
      # make outline color same as fill color - black outline will hide colors
      # for very narrow stripes
      geom_rect(aes_string(color=y, fill=y), alpha=0.6, na.rm=TRUE, ...) +
      scale_x_continuous(breaks=x.tab.df$x.center, labels=levels(data[[x]])) +
      theme(axis.ticks.x=element_blank()) +
      # y reflected in color legend; however, we do want the label when part of a multiplot!
      xlab(x)  + ylab(y)
      theme_panel_num_y

   p <- add.color.fill(p, data, y,
                       colors.gradient=colors.gradient,
                       palette.brewer.seq=palette.brewer.seq,
                       palette.brewer.qual=palette.brewer.qual)
   p <- add.color.legend(p, data, y, 'fill')

   if (estimate.label.overlap(labels=levels(data[[x]]), breaks=x.tab.df$x.center)) {
      p <- p + theme_slanted_text_x
   }
   p
}


# spine plot of three variables
# - plot.margin.x - horizontal gap between rectangles for x-values
# - plot.margin.y - vertical gap between rectangles for y-values
# - no gaps between z-rectangles
gplt.spine3 <- function(data, x, y, z, w='NULL',
                        plot.margin.x=0.05,
                        plot.margin.y=0.02,
                        max.factor.levels=10,
                        colors.gradient=c(hcl(66,60,85), hcl(128,100,35)),
                        palette.brewer.seq='YlGn',
                        palette.brewer.qual='Set1',
                        ...) {

   # ordering and factor truncation
   data <- order.factors(data, x, y, z, w=w)
   for (v in c(x, y, z)) {
      data <- quantize(data, v, w=w, n.levels=max.factor.levels)
   }

   x.lev <- length(levels(data[[x]]))
   y.lev <- length(levels(data[[y]]))
   z.lev <- length(levels(data[[z]]))

   xyz.tab <- marginalize(data, c(x, y, z), w)
   xy.tab  <- marginalize(data, c(x, y), w)
   x.tab   <- marginalize(data, x, w)
   y.tab   <- marginalize(data, y, w)
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

   # remove empty rects
   idx<-plot.data$bottom < plot.data$top & plot.data$left < plot.data$right
   plot.data <- plot.data[idx,]

   p <- ggplot(plot.data, aes(xmin=left, xmax=right, ymin=bottom, ymax=top)) +
      # make outline color same as fill color - black outline will hide colors
      # for very narrow stripes
      geom_rect(aes_string(fill=z, col=z), alpha=0.6, na.rm=TRUE, ...) +
      scale_y_continuous(breaks=y.tab.df$y.center, labels=levels(data[[y]])) +
      scale_x_continuous(breaks=x.tab.df$x.center, labels=levels(data[[x]]))

   p <- add.color.fill(p, data, z,
                       colors.gradient=colors.gradient,
                       palette.brewer.seq=palette.brewer.seq,
                       palette.brewer.qual=palette.brewer.qual)
   p <- add.color.legend(p, data, z, 'fill')

   p <- p + xlab(x) + ylab(y) +
      theme_panel_num_y +
      theme(axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank())

   # write labels slanted to mitigate overlap
   if (estimate.label.overlap(labels=levels(data[[x]]), breaks=x.tab.df$x.center)) {
      p <- p + theme_slanted_text_x
   }
   p
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
   p
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
plotluck.options <- function(opts,...) {
   if (missing(opts)) {
   opts <- list(
      na.rm=FALSE,
      geom='auto',
      sample.max.rows=100000,
      trans.log.thresh=2,
      n.breaks.histogram=NA,
      min.points.hex=5000,
      min.points.density=20,
      min.points.box=20,
      raster.resolution=30,
      scatter.dedupe='area',
      min.points.jitter=3,
      jitter.x=0.4,
      jitter.y=0.4,
      max.factor.levels=30,
      few.unique.as.factor=5,
      max.color.factors=3,
      max.factor.levels.violin=20,
      max.factor.levels.spine.x=20,
      max.factor.levels.spine.y=10,
      max.factor.levels.spine.z=5,
      spine.plot.margin.x=0.05,
      spine.plot.margin.y=0.02,
      facet.max.cols=10,
      facet.max.rows=10,
      facet.num.wrap=6,
      facet.num.grid=3,
      prefer.factors.vert=TRUE,
      fill.default='deepskyblue',
      colors.gradient=c(hcl(66,60,85), hcl(128,100,35)),
      palette.brewer.seq='YlGn',
      palette.brewer.qual='Set1',
      multi.entropy.order=TRUE,
      multi.max.rows=6,
      multi.max.cols=6,
      multi.in.grid=TRUE)
   }
   overrides <- list(...)
   unknown <- setdiff(names(overrides), names(opts))
   if (length(unknown) > 0) {
      stop(sprintf('Unknown options: %s', paste(unknown, sep='',collapse =', ')))
   }
   opts[names(overrides)] <- overrides
   opts
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
   data
}

# expected input formula: 'x', 'x*y', 'x+y', 'x+1'
# output: list of occurring variables, ignoring constants
parse.formula.term <- function(form) {
   if (class(form) == 'numeric') {
      return(character())
   }
   if (class(form) == 'name') {
      return(as.character(form))
   }
   if (as.character(form[[1]]) %in% c('+', '*')) {
      res <- character()
      if (class(form[[2]]) == 'name') {
         res <- as.character(form[[2]])
      } else if (class(form[[2]]) != 'numeric') {
         stop('invalid formula: at most two dependent or conditional variables allowed')
      }
      if (class(form[[3]]) == 'name') {
         res <- c(res, as.character(form[[3]]))
      } else if (class(form[[3]]) != 'numeric') {
         stop('invalid formula: at most two dependent or conditional variables allowed')
      }
      return(res)
   }
   stop('invalid formula')
}

# expected input: conditional formula with up to three variables
# output: list of lists (response, independent variable, conditionals)
parse.formula <- function(form) {
   if (as.character(form[[1]]) != '~') {
      stop('invalid formula')
   }
   node <- form[[2]]
   if (class(node) != 'name') {
      stop('invalid formula: only one dependent variable allowed')
   }
   resp <- as.character(node)
   indep <- character()
   cond <- character()

   node <- form[[3]]
   if (class(node) != 'name' && as.character(node[[1]])  == '|') {
      indep <- parse.formula.term(node[[2]])
      cond <- parse.formula.term(node[[3]])
   } else {
      indep <- parse.formula.term(node)
   }

   if (resp != '.' && (resp %in% indep || resp %in% cond)) {
     stop('invalid formula: variable cannot be both dependent and independent')
   }
   if (length(intersect(indep,cond))>0) {
     stop('variables can only be used once')
   }
   if (length(indep)+length(cond) > 2) {
     stop('invalid formula: at most 3 variables allowed')
   }
   return(list(resp, indep, cond))
}

# try to match user input to column names, if not exact
correct.varnames <- function(x, data) {
   if (length(x) == 0 || x == '.') {
      return(x)
   }
   for (i in seq(length(x))) {
      idx <- pmatch(tolower(x[i]), tolower(names(data)))
      if (!is.na(idx)) {
         x[i] <- names(data)[idx]
      } else {
         stop(sprintf('no such variable: %s', x[i]))
      }
   }
   x
}

# process weights
process.weights <- function(data, weights) {
   if (weights == 'NULL') {
      return(data)
   }

   if  (!is.numeric(data[[weights]])) {
      stop('Weight must be numeric')
   }
   if  (any(data[[weights]]<0)) {
      stop('Weight must be non-negative')
   }
   weight.na <- is.na(data[[weights]])
   if  (any(weight.na)) {
      warning('Weight is NA for %d instances, deleting', length(which(weight.na)))
      data <- data[!weight.na,]
   }
   # if weights are integer type, Hmisc::wtd.quantile() can lead to NA due to overflow
   data[[weights]] <- as.double(data[[weights]])
   data
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
#' to \code{facet.max.rows} and \code{facet.max.cols})  that facilitates
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
#'  weight. By default, the option \code{scatter.dedupe} is
#'  switched on, which does exactly this. If it is \code{FALSE}, horizontal and
#'  vertical jittering will be applied if the number of duplicated points
#'  exceeds \code{min.points.jitter}. The amount of jittering can be controlled
#'  with \code{jitter.x} and \code{jitter.y}.
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
#'  comprises more than \code{sample.max.rows}, it will be sampled down to that
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
#'  \code{max.factor.levels.violin}) violin plots in a row, the algorithm switches to
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

plotluck <- function(formula, data, weights,
                          opts=plotluck.options(),
                          ...) {
   parsed <- parse.formula(formula)
   response <- correct.varnames(parsed[[1]], data)
   indep <- correct.varnames(parsed[[2]], data)
   cond <- correct.varnames(parsed[[3]], data)

   vars <- c(response, indep, cond)
   if (!missing(weights)) {
      weights <- correct.varnames(deparse(substitute(weights)), data)
      vars <- c(vars, weights)
   } else {
      # note: we set missing aesthetics to 'NULL' instead of NULL, because
      # then ggplot interprets it correctly as a constant [e.g., aes(weights='NULL')]
      weights <- 'NULL'
   }
   vars <- unique(vars)

   if ('.' %in% vars) {
      # trellis of plots
      return(plotluck.multi(response, indep, data, w=weights,
                            in.grid=opts$multi.in.grid,
                            max.rows=opts$multi.max.rows,
                            max.cols=opts$multi.max.cols,
                            opts=opts,
                            ...))
   }

   # validation and preprocessing
   data <- data[,vars, drop=FALSE]
   data <- process.weights(data, weights)
   data <- preprocess.factors(data, opts$na.rm)
   data <- sample.data(data, weights, opts$sample.max.rows)

   for (v in c(response, indep)) {
      if (is.numeric(data[[v]])) {
         data <- discretize.few.unique(data, v, opts$few.unique.as.factor)
      }
   }

   # conditionals: discretize if numeric, and order
   # note: to have a single data frame containing all plot and facet variables,
   # we have to preprocess the conditionals beforehand.

   if (length(cond) > 0) {

      # if we have to discretize, we want that many facets
      intervals <- opts$facet.num.wrap
      if (length(cond) == 2) {
         intervals <- opts$facet.num.grid
      }

      order.by <- NULL
      if (is.numeric(data[[response]]) || is.ordered(data[[response]])) {
         order.by <- response
      } else if (length(indep) == 1 &&
                 (is.numeric(data[[indep]]) || is.ordered(data[[indep]]))) {
         order.by <- indep
      }

      for (i in seq(length(cond))) {
         if (is.numeric(data[[cond[i]]])) {
            data <- discretize(data, cond[i], w=weights, max.breaks=intervals,
                                          estimate.breaks=FALSE, method='histogram')
         } else {
            data <- limit.factor.levels(data, cond[i], w=weights,
                                         max.levels=opts$max.factor.levels)
            if (!is.ordered(data[[cond[i]]])) {
               data <- order.factor.by.value(data, cond[i], order.by, weights)
            }
         }
      }
   }

   vars.non.numeric <- cond
   vars.numeric <- NULL
   for (v in c(response, indep)) {
      if (!is.numeric(data[[v]])) {
         vars.non.numeric <- c(vars.non.numeric, v)
      } else {
         vars.numeric <- c(vars.numeric, v)
      }
   }

   # compute joint number of factor combinations, and
   # the maximum/median number of instances in them
   if (length(vars.non.numeric) > 0) {
      t <- table(data[,vars.non.numeric], useNA='ifany')
      num.groups   <- length(t)
      med.points.per.group <- median(t)
      n.max.level <- max(t)
   } else {
      num.groups <- 1
      med.points.per.group <- nrow(data)
      n.max.level <- nrow(data)
   }

   # precompute medians
   if (length(vars.numeric) == 1) {
      grp.med <- group.central(data, vars.numeric, vars.non.numeric, w=weights, method='median')
      data <- merge(data, grp.med)
   }

   # note: factor levels cannot be limited before ordering is known!

   # determine type of plot
   if (length(indep) == 0) {
      # distribution plot
      res <- determine.plot.type.0(data, response, weights,
                                   med.points.per.group, num.groups, opts, ...)
   } else if (length(indep) == 1) {
      res <- determine.plot.type.1(data, response, indep, cond, weights,
                                   med.points.per.group, opts, ...)
   } else { # length(indep) == 2
      res <- determine.plot.type.2(data, response, indep[[1]], indep[[2]],
                                   weights, opts, ...)
   }

   p <- res[[1]]
   type.plot <- res[[2]]

   if (length(indep) < 2) {
      p <- add.conditional.layer(p, data, response, indep, cond, type.plot, opts)
   }
   p <- p + theme(panel.background=element_blank(), # get rid of gray background
               strip.background=element_blank(), # no background for facet labels
               axis.line.x=element_line(),
               axis.line.y=element_line(),
               panel.grid=element_line(color='black'))
   p
}

determine.plot.type.0 <- function(data, response, w, med.points.per.group,
                                  num.groups, opts,
                                  ...) {
   p <- NULL
   type.plot <- NULL
   if (is.numeric(data[[response]])) {
      opts.geom <- match.arg(opts$geom, c('density', 'histogram', 'scatter', 'auto'))
      if (opts.geom == 'auto') {
         # enough points for density plot?
         opts.geom <- ifelse(med.points.per.group > opts$min.points.density, 'density', 'scatter')
      }

      if (opts.geom == 'density') {
         type.plot <- 'density'
         p <- gplt.density(data, response, w=w,
                           trans.log.thresh=opts$trans.log.thresh,
                           ...)
      } else if (opts.geom == 'histogram') {
         type.plot <- 'histogram'
         p <- gplt.histogram(data, response, w=w,
                             n.breaks=opts$n.breaks.histogram,
                             trans.log.thresh=opts$trans.log.thresh,
                             ...)
      } else { # opts.geom == 'scatter'

         # to plot a distribtion as a line plots,
         # make up a constant y coordinate
         indep <- '.y_const'
         data[[indep]] <- 1
         data[[indep]] <- as.factor(data[[indep]])
         type.plot <- 'scatter.num.1'
         if (!opts$prefer.factors.vert) {
            p <- gplt.scatter(data, indep, response, w=w,
                              scatter.dedupe=opts$scatter.dedupe,
                              min.points.jitter=opts$min.points.jitter,
                              trans.log.thresh=opts$trans.log.thresh,
                              max.factor.levels=opts$max.factor.levels,
                              ...) +
               theme(axis.ticks.x=element_blank(),
                     axis.text.x=element_blank()) +
               xlab('')
            # known issue: the conditioning layer can add theme_slanted_text_x,
            # in which case spurious '1's appear on the x-axis.
         } else {
            p <- gplt.scatter(data, response, indep, w=w,
                              scatter.dedupe=opts$scatter.dedupe,
                              min.points.jitter=opts$min.points.jitter,
                              trans.log.thresh=opts$trans.log.thresh,
                              max.factor.levels=opts$max.factor.levels,
                              ...) +
               theme(axis.ticks.y=element_blank(),
                     axis.text.y=element_blank()) +
               xlab('')
         }
      }
   } else { # !is.numeric(data[[response]])
      opts.geom <- match.arg(opts$geom, c('dot', 'bar', 'auto'))
      if (opts.geom == 'auto') {
         opts.geom <- ifelse(num.groups < 6, 'bar', 'dot')
      }
      type.plot <- opts.geom
      p <- gplt.dot(data, response, w=w, vertical=opts$prefer.factors.vert,
                    max.factor.levels=opts$max.factor.levels, geom=opts.geom, ...)
   }
   list(p, type.plot)
}

determine.plot.type.1 <- function(data, response, indep, cond,
                                  w, med.points.per.group, opts,
                                  ...) {
   p <- NULL
   type.plot <- NULL
   if (is.numeric(data[[indep]]) && is.numeric(data[[response]])) {
      opts.geom <- match.arg(opts$geom, c('hex', 'scatter', 'auto'))
      if (opts.geom == 'auto') {
         # if a lot of points, choose hex
         opts.geom <- ifelse(nrow(data) >= opts$min.points.hex, 'hex', 'scatter')
      }
      # HACK: the smoothing line can only be colored with a default color if it
      # wouldn't be used to distinguish multiple groups.
      fill.smooth <- ifelse(length(cond) == 0, opts$fill.default, 'NULL')
      if (opts.geom == 'hex') {
         type.plot <- 'hex'
         p <- gplt.hex(data, indep, response, w=w, trans.log.thresh=opts$trans.log.thresh,
                       # hex plot already has color, make the smoothing line neutral
                       fill.smooth='grey', ...)
      } else {
         type.plot <- 'scatter.num'
         p <- gplt.scatter(data, indep, response, w,
                           scatter.dedupe=opts$scatter.dedupe,
                           min.points.jitter=opts$min.points.jitter,
                           jitter.x=opts$jitter.x,
                           jitter.y=opts$jitter.y,
                           trans.log.thresh=opts$trans.log.thresh,
                           max.factor.levels=opts$max.factor.levels,
                           fill.smooth=fill.smooth,
                           ...)
      }
   } else if (!is.numeric(data[[indep]]) && !is.numeric(data[[response]])) {
      opts.geom <- match.arg(opts$geom, c('spine', 'raster', 'auto'))

      if (opts.geom == 'auto') {
         opts.geom <- ifelse(length(cond) == 0 && # HACK: spine plot cannot be faceted, due to difficulty of implementation
                                length(levels(data[[response]])) <= opts$max.factor.levels.spine.y &&
                                length(levels(data[[indep]])) <= opts$max.factor.levels.spine.x, 'spine', 'raster')
      }

      if (opts.geom == 'spine') {
         type.plot <- 'spine'
         p <- gplt.spine(data, indep, response, w=w,
                         plot.margin.x=opts$spine.plot.margin.x,
                         max.factor.levels=opts$max.factor.levels.spine.x,
                         colors.gradient=opts$colors.gradient,
                         palette.brewer.seq=opts$palette.brewer.seq,
                         palette.brewer.qual=opts$palette.brewer.qual,
                         ...)
      } else {
         type.plot <- 'raster'
         p <- gplt.raster(data, indep, response, z='NULL', w=w,
                          trans.log.thresh=opts$trans.log.thresh,
                          max.factor.levels=opts$max.factor.levels,
                          raster.resolution=opts$raster.resolution,
                          colors.gradient=opts$color.gradient,
                          palette.brewer.seq=opts$palette.brewer.seq,
                          palette.brewer.qual=opts$palette.brewer.qual,
                          ...)
      }
   } else {
      # mixed types: one factor, one numeric
      # assignment of variables to axes depends only on opts$prefer.factors.vert

      opts.geom <- match.arg(opts$geom, c('violin', 'box', 'scatter', 'auto'))

      var.fac <- response
      var.num <- indep
      switch <- is.numeric(data[[response]]) != opts$prefer.factors.vert
      if (switch) {
         var.num <- response
         var.fac <- indep
      }
      if (opts.geom == 'auto') {
         # if there are too many violin plots in a horizontal row, they
         # just look like black lines
         opts.geom <- ifelse(med.points.per.group > opts$min.points.box,
                             ifelse(length(levels(var.fac)) <= opts$max.factor.levels.violin, 'violin', 'box'),
                             'scatter')
      }

      if (opts.geom == 'violin') {
         type.plot <- 'violin'
         p <- gplt.violin(data, var.fac, var.num, w=w,
                          trans.log.thresh=opts$trans.log.thresh,
                          max.factor.levels=opts$max.factor.levels,
                          ...)
      } else if (opts.geom == 'box') {
         type.plot <- 'box'
         p <- gplt.box(data, var.fac, var.num, w=w,
                       trans.log.thresh=opts$trans.log.thresh,
                       max.factor.levels=opts$max.factor.levels,
                       ...)
      } else { # opts.geom == 'scatter'
         # not enough points for violin or box plot
         type.plot <- 'scatter.mixed'
         p <- gplt.scatter(data, var.fac, var.num, w=w,
                           scatter.dedupe=opts$scatter.dedupe,
                           min.points.jitter=opts$min.points.jitter,
                           trans.log.thresh=opts$trans.log.thresh,
                           max.factor.levels=opts$max.factor.levels,
                           ...)
      }
   }
   list(p, type.plot)
}

determine.plot.type.2 <- function(data, response, indep1, indep2, w, opts,
                                  ...) {
   p <- NULL
   type.plot <- NULL
   opts.geom <- match.arg(opts$geom, c('spine', 'raster', 'auto'))
   if (opts.geom == 'auto') {

      opts.geom <- ifelse((!is.numeric(data[[response]])) && (!is.numeric(data[[indep1]])) &&
                             (!is.numeric(data[[indep2]])) &&
                             length(levels(data[[response]])) <= opts$max.factor.levels.spine.z &&
                             length(levels(data[[indep2]])) <= opts$max.factor.levels.spine.y &&
                             length(levels(data[[indep1]])) <= opts$max.factor.levels.spine.x, 'spine', 'raster')
   }

   if (opts.geom == 'spine') {
      type.plot <- 'spine3'
      p <- gplt.spine3(data, indep1, indep2, response, w=w,
                       plot.margin.x=opts$spine.plot.margin.x,
                       plot.margin.y=opts$spine.plot.margin.y,
                       max.factor.levels=opts$max.factor.levels.spine.x,
                       colors.gradient=opts$color.gradient,
                       palette.brewer.seq=opts$palette.brewer.seq,
                       palette.brewer.qual=opts$palette.brewer.qual,
                       ...)
   } else {
      type.plot <- 'raster3'
      p <- gplt.raster(data, indep1, indep2, response, w=w,
                       trans.log.thresh=opts$trans.log.thresh,
                       max.factor.levels=opts$max.factor.levels,
                       raster.resolution=opts$raster.resolution,
                       colors.gradient=opts$colors.gradient,
                       palette.brewer.seq=opts$palette.brewer.seq,
                       palette.brewer.qual=opts$palette.brewer.qual,
                       ...)
   }
   list(p, type.plot)
}


# format factor levels with the first or last value replaced by <var>=<val>
# this helps save space on the grid display
format.facets <- function(data, x, show.var='first') {
   lev <- NULL
   if (is.numeric(data[[x]])) {
      lev <- as.character(sort(unique(data[[x]])))
   } else {
      lev <- levels(data[[x]])
   }
   names(lev) <- lev
   idx <- ifelse(show.var == 'first', 1, length(lev))
   lev[idx] <- sprintf('%s = %s', x, lev[idx])
   lev
}

add.facet.wrap <- function(p, data, cond, preferred.order, opts) {
   nrow <- NULL
   ncol <- NULL
   switch <- NULL
   show.var <- 'first'
   if (preferred.order == 'row') {
      ncol <- opts$facet.max.cols
      switch <- 'x'
   } else if (preferred.order == 'col') {
      nrow <- opts$facet.max.rows
      switch <- 'y'
      show.var <- 'last'
   }
   facet.labels <- format.facets(data, cond, show.var)

   p <- p + theme_slanted_text_x +  # axis text might overlap in small diagrams
      facet_wrap(cond,
                 labeller=as_labeller(facet.labels),
                 nrow=nrow, ncol=ncol, switch=switch) +
      theme(panel.border=element_rect(fill=NA))
   # known issue 1: always slanting text is not ideal, but no quick fix
   # known issue 2: for it would be good to order vertical facets like axis
   # (i.e., the lowest value at the bottom), but doesn't seem to be easy to configure.
}


add.facet.grid <- function(p, data, cond) {
   facet.labels.1 <- format.facets(data, cond[1], 'first')
   facet.labels.2 <- format.facets(data, cond[2], 'first')

   p + theme_slanted_text_x + # axis text can easily overlap in small diagrams
      facet_grid(sprintf('%s~%s', cond[1], cond[2]),
                 labeller=labeller(.rows=as_labeller(facet.labels.1),
                                   .cols=as_labeller(facet.labels.2))) +
      theme(panel.border=element_rect(fill=NA))
}


# when we get here, color will not be used for further distinction;
# no harm in using previous factor for (redundant) coloring
redundant.factor.color <- function(p, data, response, indep, type.plot, opts) {
   if (type.plot %in% c('spine', 'spine3', 'hex')) {
      # these types are already colored
      return(p)
   }
   fac <- NULL
   if (!is.numeric(data[[response]])) {
      fac <- response
   } else if  (length(indep) > 0 && !is.numeric(data[[indep]])) {
      fac <- indep
   }
   if (!is.null(fac)) {
      p <- add.color.fill(p, data, fac,
                          colors.gradient=opts$colors.gradient,
                          palette.brewer.seq=opts$palette.brewer.seq,
                          palette.brewer.qual=opts$palette.brewer.qual)  +
         guides(fill=FALSE, color=FALSE)

      return (p)
   } else if (type.plot %in% c('density', 'histogram')) {
      # color 1D density plots with default color
      p + aes(fill='something') +
         aes(color='something') +
         scale_fill_manual(values=opts$fill.default) +
         scale_color_manual(values=opts$fill.default) +
         guides(fill=FALSE, color=FALSE)
   } else {
      p
   }
}

# represent conditional variables either by color or facets
add.conditional.layer <- function(p, data, response, indep, cond, type.plot, opts) {

   if (length(cond) == 0) {
      return(redundant.factor.color(p, data, response, indep, type.plot, opts))
   }

   # heuristic for facet layout
   preferred.order <- 'none'
   if (type.plot %in% c('density', 'histogram')) {
      preferred.order <- 'col'
   } else if (!(type.plot %in% c('hex', 'scatter.num', 'raster'))) { # pure numeric -> no preference
      # mixed num/factor -> opposite to graph layout
      preferred.order <- ifelse(opts$prefer.factors.vert, 'row', 'col')
   }

   if (length(cond) == 1) {
      u <- length(unique(data[[cond]]))

      if (u <= opts$max.color.factors && !(type.plot %in% c('spine', 'raster', 'hex'))) {
         p <- add.color.fill(p, data, cond,
                             colors.gradient=opts$colors.gradient,
                             palette.brewer.seq=opts$palette.brewer.seq,
                             palette.brewer.qual=opts$palette.brewer.qual)

         aesth.legend <- ifelse(type.plot %in% c('density', 'histogram', 'bar'), 'fill', 'color')
         p <- add.color.legend(p, data, cond, aesth.legend)

      } else {
         p <- redundant.factor.color(p, data, response, indep, type.plot, opts)
         add.facet.wrap(p, data, cond, preferred.order, opts)
      }
   } else { #length(cond) == 2
      #u1 <- length(unique(data[[cond[1]]]))
      #u2 <- length(unique(data[[cond[2]]]))
      # if (min(u1,u2) <= opts$max.color.factors) {
      #    if (u1 < u2) {
      #       p + aes_string(fill=cond[1]) + aes_string(color=cond[1]) +
      #          add.facet.wrap(p, data, cond[2], preferred.order, opts)
      #    } else {
      #       p + aes_string(fill=cond[2]) + aes_string(color=cond[2]) +
      #          add.facet.wrap(p, data, cond[1], preferred.order, opts)
      #    }
      # } else {
      # example: plotluck(price~1|cut+clarity, diamonds)
      p <- redundant.factor.color(p, data, response, indep, type.plot, opts)
      add.facet.grid(p, data, cond)
      #}
   }
}


# plot multiple graphs in a grid layout, possibly over multiple pages
# suppress.*lab == 'all': no axis labels for this dimension
# suppress.*lab == 'margin' means:
# no x-axis labels except at the bottom;
# no y-axis labels except in the leftmost plots
mplot <- function(plots, rows=ceiling(sqrt(length(plots))),
                  cols=ceiling(sqrt(length(plots))),
                  suppress.xlab=FALSE, suppress.ylab=FALSE) {

   num.plots <- length(plots)
   if (cols > num.plots) {
      cols <- num.plots
      rows <- 1
   }

   size.page <- rows * cols

   layout <- matrix(seq( size.page),
                    ncol=cols, nrow=rows, byrow=TRUE)

   for (p in 1:ceiling(num.plots/size.page)) {

      grid::grid.newpage()
      grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

      # make each plot, in the correct location
      for (i in 1:min(size.page, num.plots-(p-1)*size.page)) {
         # get the i,j matrix positions of the regions that contain this subplot
         matchidx <- as.data.frame(which(layout==i, arr.ind=TRUE))
         plot.idx <- (p-1)*size.page+i
         plot.current <- plots[[plot.idx]]

         if (suppress.xlab == 'all' ||
             (suppress.xlab == 'margin' && matchidx$row != rows
              && plot.idx + cols <= num.plots)) { # last complete row

            # note: the following doesn't work if a coord_flip() was applied:
            # plot.current <- plot.current + xlab('')
            # axis labels move with the coordinates, themes don't!
            plot.current <- plot.current + theme(axis.title.x=element_blank())
         }
         if (suppress.ylab == 'all' ||
             (suppress.ylab == 'margin' && matchidx$col != 1)) {
            plot.current <- plot.current + theme(axis.title.y=element_blank())
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


plotluck.multi <- function(response, indep, data, w='NULL',
                           in.grid=TRUE,
                           max.rows=10, max.cols=10,
                           opts=plotluck.options(),
                           ...) {

   main <- deparse(substitute(data))

   # if data size too large, apply sampling here; expensive to repeat
   data <- sample.data(data, w, opts$sample.max.rows)

   vars.x <- NULL
   vrs.y <- NULL
   if (response == '.') {
      vars.y <- names(data)
      if (length(indep) == 0) {
         vars.x <- '1'
      } else if (indep == '.') {
         vars.x <- names(data)
      }  else {
         vars.x <- indep
      }
   } else { # response != '.'
      vars.y <- response
      if(indep != '.') {
         stop('error: invalid formula')
      }
      vars.x <- names(data)
   }

   # order variables by conditional entropy
   cond.ent <- NULL
   if (opts$multi.entropy.order) {
      if (length(vars.y) > 1 && length(vars.x) == 1 && vars.x != '1') {
         cond.ent <- cond.entropy.data(data, given=vars.x, w=w)
         vars.y <- vars.y[order(cond.ent[,vars.x])]
      } else if (length(vars.y) == 1 && length(vars.x) > 1) {
         cond.ent <- cond.entropy.data(data, target=vars.y, w=w)
         vars.x <- vars.x[order(cond.ent[vars.y,])]
      } else if (length(vars.y) > 1 && length(vars.x) > 1) {
         cond.ent <- cond.entropy.data(data, w=w)
         var.ent <- (apply(cond.ent, 1, mean) + apply(cond.ent, 1, mean))/2
         vars.x <- vars.x[order(var.ent)]
         vars.y <- vars.x
      }
   }

   combi <- expand.grid(seq(length(vars.x)),
                               seq(length(vars.y)), stringsAsFactors=FALSE)
   names(combi) = c('grid.x', 'grid.y')
   combi$x <- vars.x[combi$grid.x]
   combi$y <- vars.y[combi$grid.y]
   combi$fac.x <- sapply(combi$x, function(v) !is.numeric(data[[v]]))
   combi$fac.y <- sapply(combi$y, function(v) !is.numeric(data[[v]]))

   # For mixed factor/num combinations, plotluck() assigns the axes only
   # depending on opts$prefer.factors.vert. In the output lattice for
   # plotluck.multi, we need to control the direction (and show opposites
   # at mirror position)
   combi$opts <- 'prefer.factors.vert=TRUE'
   combi$labs <- ''
   if (length(vars.x) > 1 && length(vars.y) == 1) {
      if (is.numeric(data[[vars.y]])) {
         combi$opts <-
            'prefer.factors.vert=FALSE'
      }
   } else  if (length(vars.x) == 1 && length(vars.y) > 1) {
      if (!is.numeric(data[[vars.x]])) {
         combi$opts <-
            'prefer.factors.vert=FALSE'
      }
   } else {  # length(vars.x) > 1 && length(vars.y) > 1
      below.diag <- combi$grid.x > combi$grid.y
      combi$opts[below.diag & combi$fac.x] <-
         'prefer.factors.vert=FALSE'
      combi$opts[(!below.diag) & (!combi$fac.y)] <-
         'prefer.factors.vert=FALSE'

      # label the corner elements with the variable name, instead
      # of 'density' or 'count' (diagonal contains distribution plots)
      idx.id <- combi$x == combi$y
      combi$labs[idx.id] <-
         sprintf('+xlab("%s")+ylab("%s")', combi$x[idx.id],
                combi$y[idx.id])
   }

   # try to make a square layout
   cols <- ceiling(sqrt(nrow(combi)))
   rows <- ceiling(nrow(combi)/cols)

   # figure out which plots are not at the margins, and
   # don't show axis labels for them if redundant
   suppress.xlab <- FALSE
   suppress.ylab <- FALSE

   theme.multi <- ''

   is.square <- TRUE

   if (in.grid) {

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

      # area scale slow and hardly visible
      opts$scatter.dedupe <- 'jitter'

      # use the limited space to make subgraphs
      # relatively rectangular
      opts$facet.max.rows <- NULL
      opts$facet.max.cols <- NULL

      # do all plots fit on a single page?
      if (cols > max.cols || cols > max.rows) {
         is.square <- FALSE
         cols <- max.cols
         rows <- min(max.rows, ceiling(nrow(combi)/cols))
      }

      if (length(vars.x) > 1 && length(vars.y) > 1 && is.square) {
          # if the full cross product fits on one page,
          # write the axis labels only on the margins
          suppress.xlab <- 'margin'
          suppress.ylab <- 'margin'
      }
      if (length(vars.x) == 1) {
         # do not repeat the axis label for the constant dimension
         if (vars.x != '1') {
            main <- vars.x
            suppress.xlab <- 'margin'
         }
         else {
            # distribution/density plot
            suppress.ylab <- 'all'
         }
      }
      if (length(vars.y) == 1) {
         if (vars.y != '1') {
            main <- y
            suppress.ylab <- 'margin'
         } else {
            suppress.ylab <- 'all'
         }
      }
   }

   if (length(vars.x) > 1 && length(vars.y) > 1 &&
       ((!in.grid) || (!is.square))) {
     # plots with the axis transposed are redundant;
     # if not part of grid, no need to show them
     combi <- combi[combi$grid.x >= combi$grid.y,]
   }

   # plot distributions for x~x
   combi$x[combi$x == combi$y] <- '1'

   if (w != 'NULL')
   {
      call.strs <- sprintf('plotluck(%s~%s, data, w=%s, opts=plotluck.options(opts,%s), ...)%s%s',
                           combi$y, combi$x, w, combi$opts, combi$labs, theme.multi)
   } else {
      call.strs <- sprintf('plotluck(%s~%s, data, opts=plotluck.options(opts,%s), ...)%s%s',
                           combi$y, combi$x, combi$opts, combi$labs, theme.multi)
   }

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


#'Plot a dataset with \code{plotluck} for a randomly generated formula
#'
#'@param data a data frame
#'@param ... additional parameters to be passed to \code{plotluck}, such as
#'       \code{weights} and \code{opts}.
#'@return a ggplot2 object.
#'
#'  \code{sample.plotluck} samples the formula by (1) uniformly drawing the
#'  number of variables (1-3) and then, for each variable, (2) uniformly
#'  choosing one of the existing types (numeric, ordered or unordered factor),
#'  and (3) uniformly selecting one of the columns of that type.
#'
#'@seealso \code{\link{plotluck}}
#'@export
#'@examples
#'set.seed(42)
#' data(iris)
#' sample.plotluck(iris)
sample.plotluck <- function(data, ...) {
   idx.ord <- which(sapply(data, is.ordered))
   idx.fac <- which(sapply(data, function(x) {is.factor(x) && (!is.ordered(x))}))
   idx.num <- setdiff(setdiff(seq(length(data)),idx.fac), idx.ord)
   w <- match.call()$weights
   if (!is.null(w)) {
      idx.w <- which(names(data) == w)
      idx.num <- setdiff(idx.num, idx.w)
   }

   types <- list()
   if (length(idx.fac) > 0) {
      types[[length(types)+1]] <- idx.fac
   }
   if (length(idx.ord) > 0) {
      types[[length(types)+1]] <- idx.ord
   }
   if (length(idx.num) > 0) {
      types[[length(types)+1]] <- idx.num
   }
   form.length <- sample(c(1,2,3),1)
   cond.pos <- sample(seq(form.length), 1) + 0.5
   form <- NULL
   op <- NULL
   for (pos in seq(form.length)) {
      form <- c(form, op)
      if (pos > cond.pos) {
         cond.pos <- 1000
         op <- '|'
         if (pos == 1) {
            form <- c(form, '~', '1')
         }
      } else if (pos == 1) {
         op <- '~'
      } else {
         op <- '+'
      }
      while (TRUE) {
         var.type <- sample(seq(length(types)), 1)
         var <- sample(seq(length(types[[var.type]])),1)
         var <- names(data)[(types[[var.type]][var])]
         if (!(var %in% form)) {
            # avoid duplication of variables
            break
         }
      }
      form <- c(form, var)
   }
   if (length(form) == 1) {
      form <- c(form, '~', '1')
   }
   form <-do.call(paste, as.list(form))

   # collect extra arguments
   arg <- sapply(match.call()[seq(2,length(match.call()))], deparse)
   s <- list()
   for (i in seq(length(arg))) {
      s[[length(s)+1]] <- sprintf('%s = %s', names(arg)[i], arg[i])
   }
   s <- do.call(paste,c(s,sep=', '))
   s <- sprintf('plotluck(%s, %s)', form, s)
   print(s)
   plotluck(as.formula(form), data, ...) + labs(title=s)
}
