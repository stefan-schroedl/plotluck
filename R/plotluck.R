#' @import ggplot2
#' @importFrom plyr ddply summarise .
#' @importFrom stats as.formula median quantile runif xtabs
#' @importFrom grDevices colorRampPalette hcl
#' @importFrom stats na.omit na.pass
#'
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

get.trans.fun <- function(data, x, expansion.thresh=2, label=x, verbose=FALSE) {

  if (!is.numeric(data[[x]]) ||
      length(unique(data[[x]])) <= 3) {
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

  if (ratio.before == 0) {
     ratio.after <- 1
  } else {
     range.core.after     <- fun.trans$trans(range.core)
     range.complete.after <- fun.trans$trans(range.complete)
     ratio.after          <- diff(range.core.after) / diff(range.complete.after)
  }
  if (verbose) {
    if (ratio.after <= expansion.thresh * ratio.before) {
      cat(sprintf('Not applying logarithmic axis scaling for %s; expansion ratio is %f, trans.log.thresh = %f\n',
                  x, ratio.after/ratio.before, expansion.thresh))
    } else {
      cat(sprintf('Applying logarithmic axis scaling for %s; expansion ratio is %f, trans.log.thresh = %f\n',
                  x, ratio.after/ratio.before, expansion.thresh))
    }
  }
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
# ordered.as.num - treat ordered factors as integers, and return a numeric average
# (in the case that the centers are used for sorting only, this prevents the rounding loss)
group.central <- function(data, x, group.vars, w='NULL', method=c('median', 'mean'),
                          col.name='.center.', ordered.as.num=FALSE) {

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

  if (w == 'NULL' || length(unique(data[[w]])) == 1) {
    if (is.num || is.ord) {
      fun <- ifelse(method=='median', median, mean)
      grp.center <- plyr::ddply(data, group.vars, function(d) fun(d[[x]], na.rm=TRUE))
    } else {
      # known issue: for ordering, it would be better to take also the frequency of the mode
      # into account
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

  # for ordinals, reattach levels, but keep order
  if (is.ord && (!ordered.as.num)) {
    grp.center[[length(grp.center)]] <- ordered(lev[round(grp.center[[length(grp.center)]])], levels=lev, exclude=NULL)
  }

  names(grp.center)[length(grp.center)] <- col.name
  grp.center
}

# make all values of data[[x]] identically equal to the group mean
replace.by.central <- function(data, x, g, w) {
  x.avg <- group.central(data, x, g, w, method='mean', col.name='.replace.center.')
  data <- merge(data, x.avg)
  data[[x]] <- data[[length(data)]]
  data[[length(data)]] <- NULL
  data
}

# return data frame with the levels of unordered factor x resorted by (weighted) frequency
order.factor.by.freq <- function(data, x, w='NULL', decreasing=FALSE, verbose=FALSE) {
  x.data <- data[[x]]
  if (is.numeric(x.data) || is.ordered(x.data)) {
    # do not change ordered factors
    return(data)
  }

  if (w=='NULL') {
    if (verbose) {
      cat(sprintf('Ordering %s levels by number of lines\n', x))
    }
    tab <- table(x.data, useNA='ifany')
  } else {
    if (verbose) {
      cat(sprintf('Ordering %s levels by sum of %s\n', x, w))
    }
    tab <- xtabs(as.formula(sprintf('%s ~ %s', w, x)), data, exclude=NULL, na.action=na.pass)
  }
  ord <- order(tab, decreasing=decreasing)
  data[[x]] <- factor(x.data, levels=levels(x.data)[ord],
                      exclude=NULL)
  return(data)
}


# return data frame with the levels of unordered factor x resorted by (weighted) central tendency of dependent variable y
order.factor.by.value <- function(data, x, y, w='NULL', decreasing=FALSE, verbose=FALSE) {

  x.data <- data[[x]]
  if (is.numeric(x.data) || is.ordered(x.data) || length(unique(data[[x]])) == 1) {
    # do not change ordered factors
    return(data)
  }

  if (is.null(y) || y == 'NULL' || y == x) {
    return(order.factor.by.freq(data, x, w, decreasing))
  }

  y.data <- data[[y]]

  if (verbose) {
    cat(sprintf('Ordering %s levels by %s\n', x, y))
  }
  # we are displaying median lines for distributions, so the sorting should agree with this
  centers <- group.central(data, y, x, w, method='median', ordered.as.num=TRUE)

  # if medians are equal for two or more classes, resolve ties by using the mean
  if (length(unique(centers$.center.)) != nrow(centers)) {
    centers <- group.central(data, y, x, w, method='mean', ordered.as.num=TRUE)
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

order.factors <- function(data, x, y, z='NULL', w='NULL', verbose=FALSE) {

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
      data <- order.factor.by.value(data, v, v.sort, w, verbose=verbose)
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

  # note: there can be unused levels!
  u <- length(levels(x.data))

  if (u <= tol) {
    return(data)
  }
  message(sprintf('Factor variable %s has too many levels (%d), truncating to %d', x, u, max.levels))

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
      tab <- tab[order(tab, as.numeric(names(tab)), decreasing=FALSE)]
    }
    levels.ord <- levels(x.data)[as.integer(names(tab))]
    ll <- length(levels.ord)
    if (!reverse) {
      levels.reduced <- levels.ord[seq(max(1,ll-max.levels+1),ll)]
    } else {
      levels.reduced <- levels.ord[seq(min(max.levels,ll))]
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
    x.data <- ordered(x.data)

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
discretize.few.unique <- function(data, x, few.unique.as.factor=5, verbose=FALSE) {
  if (!is.numeric(data[[x]])) {
    return(data)
  }
  non.na <- !is.na(data[[x]])
  # works for both data frames and tibbles
  u <- nrow(unique(data[non.na,x,drop=FALSE]))

  if (u <= few.unique.as.factor &&
      nrow(data) > u * u) { # heuristic similar to histogram binning
    data[[x]] <- ordered(data[[x]], exclude=NULL)
    info.threshold(verbose, '%s has only %d levels, treating as factor', few.unique.as.factor, x, u)
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

preprocess.factors <- function(data, vars, na.rm=FALSE) {

  for (x in vars) {

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
      message('Option "na.rm=TRUE", but data contains factor(s) with an NA level. Keeping these ...')
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
   lb <- length(breaks)

   # account for some empty space at the boundaries
   mrg.left <- (breaks[2]-breaks[1])/2
   mrg.right <- (breaks[lb]-breaks[lb-1])/2
   rg <- breaks[lb] - breaks[1] + mrg.left + mrg.right

   for (i in seq(lb)) {
      space.left  <- ifelse(i==1,  2*mrg.left,  (breaks[i]-breaks[i-1]))/rg
      space.right <- ifelse(i==lb, 2*mrg.right, (breaks[i+1]-breaks[i]))/rg
      label.width <- as.numeric(grid::convertX(unit(1, 'strwidth', labels[i]),'npc'))
      if (label.width > min(space.left, space.right)) {
         return(TRUE)
      }
   }
   return(FALSE)
}


# add color and/or fill scale layers to plot
# - use qualitative (sequential) brewer scale for (un)ordered factors, if
#   palette size supports it.
#   Otherwise, extend it using colorRampPalette.
# - if too few colors, use qualitative (first few colors of brewer palette might not look good by themselves!)

get.palette <- function(x, palette.brewer.seq='YlGn', palette.brewer.qual='Set1',
                        adjust.size=FALSE) {

   is.num <- is.numeric(x)
   is.ord <- is.ordered(x)
   if (is.num) {
      u <- length(unique(x))
   } else {
      # note: there can be unused levels
      u <- length(levels(x))
   }
   if (u <= 2) {
      # treat as qualitative
      is.ord <- FALSE
      is.num <- FALSE
   }
   name <- ifelse(is.num || is.ord, palette.brewer.seq, palette.brewer.qual)

   sz <- RColorBrewer::brewer.pal.info[name, 'maxcolors']
   if (is.na(sz)) {
      stop('invalid palette name: ', name)
   }

   pal <- RColorBrewer::brewer.pal(sz, name)

   if (adjust.size) {
      if ((is.num || is.ord) && u < sz) {
         # note: the first few colors of the brewer palettes are very light, sometimes
         # hard to see. If less values are needed than the palette size, we want to
         # rather drop the lightest ones instead of the darkest ones.
         pal <- pal[seq(sz-u+1, sz)]
      }
      if ((!is.num) && u > sz) {
         # not enough colors, interpolate
         # note: maybe do something different for qualitative scale?
         pal <- colorRampPalette(pal)(u)
      }
      names(pal) <- levels(x)
   }
   pal
}

add.color.fill <- function(p, data, x, aesth=c('color', 'fill'),
                           palette.brewer.seq='YlGn',
                           palette.brewer.qual='Set1') {

   na.value <- 'darkgray'
   if ('color' %in% aesth) {
      p <- p + aes_string(color=x)
   }
   if ('fill' %in% aesth) {
      p <- p + aes_string(fill=x)
   }

   pal <- get.palette(data[[x]], palette.brewer.seq, palette.brewer.qual, TRUE)
   lp <- length(pal)

   if (is.numeric(data[[x]])) {
      if ('color' %in% aesth) {
         p <- p + scale_color_gradientn(colors=pal[c(1,lp)])
      }
      if ('fill' %in% aesth) {
         p <- p + scale_fill_gradientn(colors=pal[c(1,lp)])
      }
   } else { # !is.num

      if ('color' %in% aesth) {
         p <- p + scale_color_manual(values=pal, na.value=na.value)
      }
      if ('fill' %in% aesth) {
         p <- p + scale_fill_manual(values=pal, na.value=na.value)
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

add.axis.transform <- function(p, data, x, ax=c('x','y'), trans.log.thresh=2, verbose=FALSE) {
  if (!is.numeric(data[[x]])) {
    p
  } else {
    ax <- match.arg(ax)
    trans <- get.trans.fun(data, x, trans.log.thresh, verbose=verbose)
    if (ax == 'x') {
      # note: "scale_x_continuous(trans=trans.x[[1]], name=trans.x[[2]])" doesn't allow
      # to overwrite the label later
      p + scale_x_continuous(trans=trans[[1]]) + xlab(trans[[2]])
    } else {
      p + scale_y_continuous(trans=trans[[1]]) + ylab(trans[[2]])
    }
  }
}

# note: these layers work, but in case a log transform is applied,
# the transformed data is received here - we need the original one.
StatCenterX <- ggproto("StatCenterX", Stat,
                       required_aes = c("x"),
                       default_aes = aes(xintercept = ..x..),

                       setup_params = function(data, params) {
                          if (is.null(params$weight)) {
                             params$weight <- 1
                          }
                          params
                       },

                       compute_group = function(data, scales, weight) {
                          grp <- setdiff(names(data), c('x','weight'))
                          data <- group.central(data, 'x', grp, w='weight', method='median', col.name='x')
                          data
                       }
)

StatCenterY <- ggproto("StatCenterY", Stat,
                      required_aes = c("y", "group"),

                      setup_params = function(data, params) {
                         if (is.null(params$weight)) {
                            params$weight <- 1
                         }
                         params
                      },

                      compute_group = function(data, scales, weight) {
                         grp <- setdiff(names(data), c('y','weight'))
                         data <- group.central(data, 'y', grp, w='weight', method='median', col.name='y')
                         data$x <- data$group
                         data
                      }
)

geom_point_center <- function(mapping = NULL, data = NULL,
                              position = position_dodge(width=0.9),
                              shape = 1, color = "black", size = 2, fill = NA,  alpha = NA,
                              show.legend = NA, inherit.aes = TRUE, na.rm = TRUE, ...)
{
   layer(data = data, mapping = mapping, stat = StatCenterY, geom = GeomPoint,
         position = position, show.legend = show.legend, inherit.aes = inherit.aes,
         params = list(shape=shape, color=color, size=size, fill=fill, alpha=alpha,
                       na.rm = na.rm, ...))
}

# different from base geom_vline(), adapt color according to group, because of inherit.aes=TRUE
geom_vline_center <- function(mapping = NULL, data = NULL,
                              position = 'identity',
                              color = NULL, size = 0.7, alpha = NA,
                              linetype='dashed',
                              show.legend = NA, inherit.aes = TRUE, na.rm = TRUE, ...)
{
   layer(data = data, mapping = mapping, stat = StatCenterX, geom = GeomVline,
         position = position, show.legend = show.legend, inherit.aes = inherit.aes,
         params = list(size=size, alpha=alpha, linetype=linetype,
                       na.rm = na.rm, ...))
}

geom_vline_inheritable <- function (mapping = NULL, data = NULL, ..., xintercept, na.rm = FALSE,
          show.legend = NA)
{
   if (!missing(xintercept)) {
      data <- data.frame(xintercept = xintercept)
      mapping <- aes(xintercept = xintercept)
      show.legend <- FALSE
   }
   layer(data = data, mapping = mapping, stat = StatIdentity,
         geom = GeomVline, position = PositionIdentity, show.legend = show.legend,
         inherit.aes = TRUE, params = list(na.rm = na.rm, ...))
}


# 2D scatter plot of numeric or factor variables
# - overlay smoothing line if both variables are numeric
# - dedupe.scatter='area' - unique repeated points, count frequency
# - dedupe.scatter='jitter' - plot each repeated point separately, add
#   jitter if there are more than min.points.jitter with identical coordinates
# - counts/weights are drawn with a shaded circle of proportional area
# - jitter.x, jitter.y - the amount of jittering as a multiple of resolution
# - trans.log.thresh - threshold to decide on log-transform
# - fill.smooth - fill color for smoothing line
gplt.scatter <- function(data, x, y, w='NULL',
                         dedupe.scatter=c('area','jitter'),
                         min.points.jitter=3,
                         jitter.x=0.4,
                         jitter.y=0.4,
                         trans.log.thresh=2,
                         max.factor.levels=30,
                         fill.smooth='NULL',
                         verbose=FALSE,
                         ...) {

  dedupe.scatter <- match.arg(dedupe.scatter)
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

  data <- order.factors(data, x, y, w=w, verbose=verbose)
  for (v in c(x,y)) {
    data <- limit.factor.levels(data, v, w=w,
                                max.levels=max.factor.levels)
  }

  if (dedupe.scatter == 'area') {
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
  p <- add.axis.transform(p, data, x, 'x', trans.log.thresh, verbose=verbose)
  p <- add.axis.transform(p, data, y, 'y', trans.log.thresh, verbose=verbose)

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
  list(plot=p, data=data)
}


# hexbin plot with overlayed smoothing line
gplt.hex <- function(data, x, y, w='NULL', trans.log.thresh=2,
                     fill.smooth='NULL', verbose=FALSE, ...) {

  p <- ggplot(data, aes_string(x=x, y=y)) +
    geom_hex(na.rm=TRUE, ...)

  p <- add.smooth.line(p, w, fill.smooth)

  # axis transformation
  p <- add.axis.transform(p, data, x, 'x', trans.log.thresh, verbose=verbose)
  p <- add.axis.transform(p, data, y, 'y', trans.log.thresh, verbose=verbose)

  p <- p +
    # log scaling of color often reveals more details
    scale_fill_gradientn(colors=c(hcl(66,60,95), hcl(128,100,45)), trans='log', guide=FALSE) +
    theme(legend.position='right') +
    theme_panel_num_x + theme_panel_num_y

  list(plot=p, data=data)
}


# heat map
# point grid of discretized values, with optional discretized color (z)
# the color is a representation of an appropriate central tendency measure
# for the point bin (mode for factors, median for ordinal and numeric vectors)
gplt.heat <- function(data, x, y, z='NULL', w='NULL', trans.log.thresh=2,
                        max.factor.levels=30,
                        resolution.heat=30,
                        palette.brewer.seq='YlGn',
                        palette.brewer.qual='Set1',
                        verbose=FALSE,
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

  data <- order.factors(data, x, y, z, w, verbose=verbose)

  # quantize variables
  for (v in c(x,y)) {
    if (is.numeric(data[[v]])) {
      data <- discretize(data, v, w=w, max.breaks=resolution.heat, estimate.breaks=FALSE,
                         method='histogram')
    } else {
      data <- limit.factor.levels(data, v, w=w,
                                  max.levels=min(resolution.heat,max.factor.levels),
                                  tol=max.factor.levels)
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

    breaks <- length(get.palette(data[[z]], palette.brewer.seq,
                               palette.brewer.qual)) - 1

    if (num.z) {
      # if z distribution is very skewed, better to quantize than to histogram
      trans.z <- get.trans.fun(data, z, trans.log.thresh)
      method <- 'histogram'
      if (trans.z[[1]]$name != 'identity') {
        method <= 'quantile'
      }

      # limited by color palette
      data <- discretize(data, z, w=w, max.breaks=breaks, tol=breaks,
                         estimate.breaks=FALSE, method=method)
    } else {
      data <- limit.factor.levels(data, z, w,
                                  max.levels=breaks,
                                  tol=breaks)
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

  # scale the area of the rectangle such that:
  # - the largest one is equal to the grid length
  # - the smallest one is still visible
  # - the area is proportional to the weight

  res.x <- resolution(as.numeric(data[[x]]))
  res.y <- resolution(as.numeric(data[[y]]))

  w.max <- sqrt(max(data[[w]]))
  w.min <- sqrt(min(data[[w]]))
  w.rg <- w.max - w.min
  if (w.rg == 0) {
     sz <- rep(1, nrow(data))
  } else {
     sz.max <- 1
     sz.min <- max(0.1, w.min/w.max)
     sz <- sz.min + (sqrt(data[[w]]) - w.min) * (sz.max - sz.min)  / w.rg
  }

  data$width <- res.x * sz
  data$height <- res.y * sz

  p <- ggplot(data, aes_string(x=x, y=y, height='height', width='width', color=z, fill=z, weight=w)) +
    geom_tile(na.rm=TRUE, ...)

  # axis transformation
  # known issue: logarithmic scaling can make the rectangles overlap!
  # plotluck(name ~ sleep_total + brainwt, data = msleep)
  # disabling for now ...
  #p <- add.axis.transform(p, data, x, 'x', trans.log.thresh, verbose=verbose)
  #p <- add.axis.transform(p, data, y, 'y', trans.log.thresh, verbose=verbose)

  if (ex.z) {
    p <- add.color.fill(p, data, z,
                        palette.brewer.seq=palette.brewer.seq,
                        palette.brewer.qual=palette.brewer.qual)
    p <- add.color.legend(p, data, z, 'fill')
  }

  # show gridlines
  p <- p + theme_panel_num_x + theme_panel_num_y
  if (!num.x && estimate.label.overlap(levels(data[[x]]))) {
    p <- p + theme_slanted_text_x
  }

  list(plot=p, data=data)
}

# cleveland dot plot or bar plot with stat_bin
# no log transforms for bars; see http://www.perceptualedge.com/articles/b-eye/dot_plots.pdf
gplt.dot <- function(data, x, w='NULL', vertical=TRUE,
                     max.factor.levels=30, geom=c('dot', 'bar'), verbose=FALSE, ...) {

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
  list(plot=p, data=data)
}


# violin plot with overlayed median point
# assumption: x is factor, y is numeric
gplt.violin <- function(data, x, y, w='NULL',
                        trans.log.thresh=2,
                        max.factor.levels=30,
                        verbose=FALSE,
                        ...) {

   flip <- FALSE
   if (is.numeric(data[[x]])) {
      if (is.numeric(data[[y]])) {
         stop('Violin plot requires one factor variable')
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
      data <- order.factor.by.value(data, x, y, w=w, verbose=verbose)
   }

   data <- limit.factor.levels(data, x, w=w,
                               max.levels=max.factor.levels)

   p <- ggplot(data, aes_string(x=x, y=y, weight=w, ymax=max(y,na.rm=TRUE))) +
      geom_violin(scale='width', alpha=0.8, na.rm=TRUE, ...)

   if ('.center.' %in% names(data)) {
      # dodge does not work correctly when width is not specified
      # see https://github.com/hadley/ggplot2/issues/525
      p <- p + geom_point(mapping=aes_(y=~.center.),
                          position=position_dodge(width=0.9),
                          size=2, shape=1, na.rm=TRUE)
   }

   # '+ geom_point_center()' is unfortunately not working correctly for the case that a
   # group (y-value) AND a color is specified, e.g.,
   # plotluck(Petal.Length~Species|Petal.Width,iris,opts=plotluck.options(geom='violin',max.factor.levels.color=10))
   # In this case, compute_group() in StatCenterX is only called once per y-value.

   # axis transformation
   p <- add.axis.transform(p, data, y, 'y', trans.log.thresh, verbose=verbose)

   if (!flip) {
      p <- p + theme_panel_fac_x + theme_panel_num_y

      if (estimate.label.overlap(levels(data[[x]]))) {
         p <- p + theme_slanted_text_x
      }

   } else {
      p <- p + coord_flip() + theme_panel_fac_y + theme_panel_num_x
   }
   list(plot=p, data=data)
}


# box plot
# assumption: x is factor, y is numeric
gplt.box <- function(data, x, y, w='NULL',
                     trans.log.thresh=2,
                     max.factor.levels=30,
                     verbose=FALSE,
                     ...) {

  flip <- FALSE
  if (is.numeric(data[[x]])) {
    if (is.numeric(data[[y]])) {
      stop('Box plot requires one factor variable')
    }
    flip <- TRUE
    tmp <- x
    x <- y
    y <- tmp
  }

  if (!is.ordered(data[[x]])) {
    data <- order.factor.by.value(data, x, y, w, verbose=verbose)
  }

  data <- limit.factor.levels(data, x, w=w,
                              max.levels=max.factor.levels)

  p <- ggplot(data, aes_string(x=x, y=y, weight=w, ymax=max(y,na.rm=TRUE))) +
    geom_boxplot(alpha=0.8, na.rm=TRUE, ...)

  # axis transformation
  p <- add.axis.transform(p, data, y, 'y', trans.log.thresh, verbose=verbose)

  if (!flip) {
    p <- p + theme_panel_fac_x + theme_panel_num_y

    if (estimate.label.overlap(levels(data[[x]]))) {
      p <- p + theme_slanted_text_x
    }
  } else {
    p <- p + coord_flip() + theme_panel_fac_y + theme_panel_num_x
  }
  list(plot=p, data=data)
}


# density plot with overlayed vertical median lines
gplt.density <- function(data, x, w='NULL',
                         trans.log.thresh=2,
                         verbose=FALSE,
                         ...) {

   ylab <- 'density'
   if (w != 'NULL') {
      ylab <- sprintf('%s density', w)
   }

   p <- ggplot(data, aes_string(x=x, weight=w)) +
      geom_density(aes_(y=~..scaled..), alpha=0.6, adjust=0.5, trim=TRUE, na.rm=TRUE, ...) +
      geom_rug(na.rm=TRUE)
      # unfortunately the following does not work for log transformation -
      # the layer received the transformed data + geom_vline_center()

   if ('.center.' %in% names(data)) {
      p <- p + geom_vline_inheritable(aes_(xintercept=~.center.),
                                      linetype='dashed', size=0.7, na.rm=TRUE)
   }

   # axis transformation
   p <- add.axis.transform(p, data, x, 'x', trans.log.thresh, verbose=verbose) +
      ylab(ylab) +
      # density numbers themselves not very meaningful
      theme(axis.ticks.y=element_blank(),
            axis.text.y=element_blank()) +
      theme_panel_num_x
   list(plot=p, data=data)
}

# histogram plot with overlayed vertical median lines
gplt.histogram <- function(data, x, w='NULL',
                           n.breaks=NA,
                           trans.log.thresh=2,
                           verbose=FALSE,
                           ...) {

  if (is.numeric(n.breaks)) {
    n.breaks <- n.breaks
  } else {
    n.breaks <- nclass.FD.modified(data[[x]])
  }
  bin.width <- diff(range(data[[x]], na.rm=TRUE))/n.breaks
  p <- ggplot(data, aes_string(x=x, weight=w)) +
    geom_histogram(binwidth=bin.width, position='identity', alpha=0.6, na.rm=TRUE, ...) +
    geom_rug(na.rm=TRUE) + geom_vline_center()

  # axis transformation
  p <- add.axis.transform(p, data, x, 'x', trans.log.thresh, verbose=verbose)

  ylab <- 'count'
  if (w != 'NULL') {
    ylab <- w
  }

  p <- p + ylab(ylab) +
    theme_panel_num_x + theme_panel_num_y
  list(plot=p, data=data)
}


# helper function
norm.sum <- function(x) {
   x[is.na(x)] <- 1e-20
   x <- x / sum(x)
   x
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
                       palette.brewer.seq='YlGn',
                       palette.brewer.qual='Set1',
                       verbose=FALSE,
                       ...) {

  # ordering and factor truncation
  data <- order.factors(data, x, y, w=w, verbose=verbose)
  for (v in c(x, y)) {
    data <- quantize(data, v, w=w, n.levels=max.factor.levels)
    if (is.numeric(data[[v]])) {
      data[[v]] <- ordered(data[[v]], exclude=NULL)
    }
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
  #y.cond.df$y.cond.cum <- as.vector(apply(y.cond, 1, function(x) {c(0, cumsum(x)[1:(y.lev-1)])}))

  # note: if there is no data for an x/y combination, make y coordinates non-NA for displaying
  # 'zero-area' rects
  y.cond.df <- plyr::ddply(y.cond.df, x, plyr::here(transform), y.cond=norm.sum(y.cond))
  y.cond.df <- plyr::ddply(y.cond.df, x, transform, y.cond.cum=cumsum(y.cond))
  y.cond.df$y.cond.cum <- y.cond.df$y.cond.cum - y.cond.df$y.cond

  plot.data <- y.cond.df
  plot.data <- merge(plot.data, x.tab.df)

  plot.data$left     <- plot.data$x.mrg.cum  + plot.data$x.cnt * plot.margin.x
  plot.data$right    <- plot.data$left       + plot.data$x.mrg
  plot.data$bottom   <- plot.data$y.cond.cum
  plot.data$top      <- plot.data$bottom     + plot.data$y.cond

  # remove empty rects
  #idx <- plot.data$bottom < plot.data$top & plot.data$left < plot.data$right
  #plot.data <- plot.data[idx,]

  p <- ggplot(plot.data, aes_(xmin=~left, xmax=~right, ymin=~bottom, ymax=~top)) +
    # make outline color same as fill color - black outline will hide colors
    # for very narrow stripes
    geom_rect(aes_string(color=y, fill=y), size=1, alpha=0.6, na.rm=TRUE, ...) +
    scale_x_continuous(breaks=x.tab.df$x.center, labels=levels(data[[x]])) +
    theme(axis.ticks.x=element_blank()) +
    # y reflected in color legend; however, we do want the label when part of a multiplot!
    xlab(x)  + ylab(y) +
    theme_panel_fac_y

  p <- add.color.fill(p, data, y,
                      palette.brewer.seq=palette.brewer.seq,
                      palette.brewer.qual=palette.brewer.qual)
  p <- add.color.legend(p, data, y, 'fill')

  if (estimate.label.overlap(labels=levels(data[[x]]), breaks=x.tab.df$x.center)) {
    p <- p + theme_slanted_text_x
  }
  list(plot=p, data=data)
}

# spine plot of three variables
# - plot.margin.x - horizontal gap between rectangles for x-values
# - plot.margin.y - vertical gap between rectangles for y-values
# - no gaps between z-rectangles
gplt.spine3 <- function(data, x, y, z, w='NULL',
                        plot.margin.x=0.05,
                        plot.margin.y=0.02,
                        max.factor.levels=10,
                        palette.brewer.seq='YlGn',
                        palette.brewer.qual='Set1',
                        verbose=FALSE,
                        ...) {

  # ordering and factor truncation
  data <- order.factors(data, x, y, z, w=w, verbose=verbose)
  for (v in c(x, y, z)) {
    data <- quantize(data, v, w=w, n.levels=max.factor.levels)
    if (is.numeric(data[[v]])) {
      data[[v]] <- ordered(data[[v]], exclude=NULL)
    }
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
  # note: if there is no data for an x/y combination, make z coordinates
  # non.NA and evenly spaced for displaying 'zero-area' rects
  z.cond.df <- plyr::ddply(z.cond.df, c(y,x), plyr::here(transform), z.cond=norm.sum(z.cond))
  z.cond.df <- plyr::ddply(z.cond.df, c(y,x), transform, z.cond.cum=cumsum(z.cond))
  z.cond.df$z.cond.cum <- z.cond.df$z.cond.cum - z.cond.df$z.cond

  #z.cond.df$z.cond.cum <- as.vector(apply(z.cond, 2:1,
  #                                         function(x) {c(0, cumsum(x)[1:(z.lev-1)])}))

  plot.data <- z.cond.df
  plot.data <- merge(plot.data, y.cond.df)
  plot.data <- merge(plot.data, x.tab.df)
  plot.data <- merge(plot.data, y.tab.df)

  plot.data$left     <- plot.data$x.mrg.cum  + plot.data$z.cond.cum * plot.data$x.mrg + plot.data$x.cnt * plot.margin.x
  plot.data$right    <- plot.data$left       + plot.data$z.cond * plot.data$x.mrg
  plot.data$bottom   <- plot.data$y.cond.cum + plot.data$y.cnt * plot.margin.y
  plot.data$top      <- plot.data$bottom     + plot.data$y.cond

  p <- ggplot(plot.data, aes_(xmin=~left, xmax=~right, ymin=~bottom, ymax=~top)) +
    # make outline color same as fill color - black outline will hide colors
    # for very narrow stripes
    # size=1: make sure empty rects are still visible
    geom_rect(aes_string(color=z, fill=z), size=1, alpha=0.6, na.rm=TRUE, ...) +
    scale_y_continuous(breaks=y.tab.df$y.center, labels=levels(data[[y]])) +
    scale_x_continuous(breaks=x.tab.df$x.center, labels=levels(data[[x]]))

  p <- add.color.fill(p, data, z,
                      palette.brewer.seq=palette.brewer.seq,
                      palette.brewer.qual=palette.brewer.qual)
  p <- add.color.legend(p, data, z, 'fill')

  p <- p + xlab(x) + ylab(y) +
    theme_panel_fac_y +
    theme(axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())

  # write labels slanted to mitigate overlap
  if (estimate.label.overlap(labels=levels(data[[x]]), breaks=x.tab.df$x.center)) {
    p <- p + theme_slanted_text_x
  }
  list(plot=p, data=data)
}


# blank plot with optional message text
gplt.blank <- function(text=NULL, ...) {
  p <- ggplot(data.frame(x=c(0,1), y=c(0,1)), aes_(x=~x, y=~y)) +
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
  list(plot=p, data=NULL)
}


#' Create options structure for \code{plotluck}
#'
#' @param opts An (optional) named list to start with. Anything not specified in ... will be inherited from opts.
#' @param ... Parameters to override default settings
#' @return A named list of options, usable as argument to function \code{\link{plotluck}}.
#'
#'  \code{plotluck} accepts a list of options to modify its behavior. Calling
#'  \code{plotluck.options} without arguments produces a list with the default
#'  values. Specifying any number of attribute/value pairs overrides them
#'  selectively.
#'
#'\tabular{lll}{
#' \strong{Option}\tab\strong{Default}\tab\strong{Comment}\cr
#' \code{na.rm} \tab \code{FALSE} \tab Do not show missing factor values as separate level.\cr
#' \code{geom} \tab \code{"auto"} \tab Override type of plot; the available types for a given formula and variables can be inspected with \code{verbose=TRUE}.\cr
#' \code{sample.max.rows} \tab \code{100000} \tab If the data set has more rows than that, sample it down.\cr
#' \code{trans.log.thresh} \tab \code{2} \tab Threshold for logarithmic axis scaling. Visible magnification factor of the central region of the distribution.\cr
#' \code{n.breaks.histogram} \tab \code{NA} \tab Override the number of histogram breaks.\cr
#' \code{min.points.hex} \tab \code{5000} \tab Minimum data points required to display a hexbin plot.\cr
#' \code{min.points.density} \tab \code{20} \tab Minimum data points required to display a density or histogram plot.\cr
#' \code{min.points.violin} \tab \code{20} \tab Minimum data points required to display a violin or box plot.\cr
#' \code{resolution.heat} \tab \code{30} \tab Grid spacing for heat maps.\cr
#' \code{dedupe.scatter} \tab \code{'area'} \tab To represent multiple instances of the same coordinates in scatter plot: scale the point size, or use jitter?\cr
#' \code{min.points.jitter} \tab \code{3} \tab Minimum number of coordinate duplicates to start jittering points.\cr
#' \code{jitter.x} \tab \code{0.4} \tab Amount of jitter to apply in horizontal direction, as a fraction of resolution.\cr
#' \code{jitter.y} \tab \code{0.4} \tab  Amount of jitter to apply in vertical direction, as a fraction of resolution.\cr
#' \code{few.unique.as.factor} \tab \code{5} \tab If a numeric variable has less than that many unique values, make it an ordered factor.\cr
#' \code{max.factor.levels} \tab \code{30} \tab For factors with more than that many levels, least frequent ones will be pruned into "other".\cr
#' \code{max.factor.levels.color} \tab \code{3} \tab Maximum number of factor levels that can be represented as colors in the same plot.\cr
#' \code{max.factor.levels.violin} \tab \code{20} \tab Maximum number of levels to plot violins; rather switch to box plot. \cr
#' \code{max.factor.levels.spine.x} \tab \code{20} \tab Maximum number of levels to plot in x-direction in a spine plot. \cr
#' \code{max.factor.levels.spine.y} \tab \code{10} \tab Maximum number of levels to plot in y-direction in a spine plot. \cr
#' \code{max.factor.levels.spine.z} \tab \code{5} \tab Maximum number of levels to represent as colors in a spine plot. \cr
#' \code{spine.plot.margin.x} \tab \code{0.05} \tab Horizontal gap between rectangles in a spine plot. \cr
#' \code{spine.plot.margin.y} \tab \code{0.02} \tab Vertical gap between rectangles in a spine plot. \cr
#' \code{facet.max.cols} \tab \code{10} \tab Maximum number of facet columns for conditional variables. \cr
#' \code{facet.max.rows} \tab \code{10} \tab Maximum number of facet rows for conditional variables. \cr
#' \code{facet.num.wrap} \tab \code{6} \tab Preferred number of facets for single conditional variable. \cr
#' \code{facet.num.grid} \tab \code{3} \tab Preferred number of facets for each of two conditional variables. \cr
#' \code{prefer.factors.vert} \tab \code{TRUE} \tab In mixed numeric/factor plots, use vertical axis for the factor. \cr
#' \code{fill.default} \tab \code{"deepskyblue"} \tab Default fill color for density and histogram plots. \cr
#' \code{palette.brewer.seq} \tab \code{"YlGn"} \tab Sequential brewer palette name. \cr
#' \code{palette.brewer.qual} \tab \code{"Set1"} \tab Qualitative brewer palette name. \cr
#' \code{multi.entropy.order} \tab \code{TRUE} \tab Use estimated conditional entropy to order multi-plots. \cr
#' \code{multi.max.rows} \tab \code{6} \tab Maximum number of rows for multi-plots. \cr
#' \code{multi.max.cols} \tab \code{6} \tab Maximum number of columns for multi-plots. \cr
#' \code{multi.in.grid} \tab \code{TRUE} \tab In multi-plots, make a page with multiple plots, or generate each one separately. \cr
#' \code{verbose} \tab \code{FALSE}\tab Print information about plot types, ordering, scaling, etc. \cr}
#'
#' @note \code{plotluck}'s aim is to provide a function that is usable
#'  "out-of-the-box", with no or very little manual tweaking. If you find
#'  yourself needing to change option values repeatedly or find the presets to
#'  be suboptimal, please contact the author.
#'
#' @seealso \code{\link{plotluck}}
#' @export
#'
#' @examples
#' # list all default options
#' plotluck.options()
#'
#' data(iris)
#' # default with violin plot
#' plotluck(iris, Petal.Length~Species)
#'
#' # use box-and-whiskers plot instead
#' plotluck(iris, Petal.Length~Species, opts=plotluck.options(geom='box'))
#'
#' @export
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
      min.points.violin=20,
      min.points.jitter=3,
      jitter.x=0.4,
      jitter.y=0.4,
      resolution.heat=30,
      dedupe.scatter='area',
      few.unique.as.factor=5,
      max.factor.levels=30,
      max.factor.levels.color=3,
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
      palette.brewer.seq='YlGn',
      palette.brewer.qual='Set1',
      multi.entropy.order=TRUE,
      multi.max.rows=6,
      multi.max.cols=6,
      multi.in.grid=TRUE,
      verbose=FALSE)
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
      stop('Invalid formula: at most two dependent or conditional variables allowed')
    }
    if (class(form[[3]]) == 'name') {
      res <- c(res, as.character(form[[3]]))
    } else if (class(form[[3]]) != 'numeric') {
      stop('Invalid formula: at most two dependent or conditional variables allowed')
    }
    return(res)
  }
  stop('Invalid formula')
}

# expected input: conditional formula with up to three variables
# output: list of lists (response, independent variable, conditionals)
parse.formula <- function(form) {
  if (as.character(form[[1]]) != '~') {
    stop('Invalid formula')
  }
  node <- form[[2]]
  if (class(node) != 'name') {
    stop('Invalid formula: only one dependent variable allowed')
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
    stop('Invalid formula: variable cannot be both dependent and independent')
  }
  if (length(intersect(indep,cond))>0) {
    stop('Variables can only be used once')
  }
  if (length(indep)+length(cond) > 2) {
    stop('Invalid formula: at most 3 variables allowed')
  }
  return(list(resp, indep, cond))
}

# try to match user input to column names, if not exact
correct.varnames <- function(x, data) {
  if (length(x) == 0 || length(x) == 1 && x == '.') {
    return(x)
  }
  for (i in seq(length(x))) {
    idx <- pmatch(tolower(x[i]), tolower(names(data)))
    if (!is.na(idx)) {
      x[i] <- names(data)[idx]
    } else {
      stop(sprintf('No such variable: %s', x[i]))
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
    message('Weight is NA for %d instances, deleting', length(which(weight.na)))
    data <- data[!weight.na,]
  }
  # if weights are integer type, Hmisc::wtd.quantile() can lead to NA due to overflow
  data[[weights]] <- as.double(data[[weights]])
  data
}

info.options <- function(chosen, choices, verbose=FALSE) {
  if (verbose) {
    choices <- setdiff(choices, 'auto')
    cat(sprintf("Choosing geom='%s' out of possible options: '%s'\n", chosen, paste0(choices,collapse=', ')))
  }
}

info.threshold <- function(cond, msg, threshold, ...) {
  if (cond) {
    cat(sprintf(msg, ...), sprintf('(%s = %s)\n', deparse(match.call()$threshold), threshold))
  }
}

#' "I'm feeling lucky" for ggplot
#'
#' The purpose of \code{plotluck} is to let the user focus on \emph{what} to plot,
#' and automate the \emph{how}. Given a dependency formula with up to three
#' variables, it tries to choose the most suitable type of plot. It also automates
#' sampling large datasets, correct handling of observation weights, logarithmic
#' axis scaling, ordering and pruning of factor levels, and overlaying smoothing
#' curves or median lines.
#'
#' @param data a data frame.
#' @param formula an object of class \code{\link[stats]{formula}}: a symbolic description
#'  of the relationship of up to three variables.
#' \tabular{lll}{
#' \strong{Formula}\tab\strong{Meaning}\tab\strong{Plot types}\cr
#' \code{y~1}\tab Distribution of single variable\tab Density, histogram, scatter, dot, bar\cr
#' \code{y~x}\tab One explanatory variable\tab Scatter, hex, violin, box, spine, heat\cr
#' \code{y~x+z}\tab Two explanatory variables\tab heat, spine\cr
#' \code{y~1|z} or \code{y~x|z}\tab One conditional variable\tab Represented through coloring or facetting\cr
#' \code{y~1|x+z}\tab Two conditional variables\tab Represented through facetting\cr}
#' In addition to these base plot types, the dot symbol \code{"."} can also be used,
#' and denotes all variables in the data frame. This gives rise to a lattice or
#' series of plots (use with caution, can be slow).
#' \tabular{lll}{
#' \strong{Formula}\tab\strong{Meaning}\cr
#' \code{.~1}\tab Distribution of each variable in the data frame, separately\cr
#' \code{y~.}\tab Plot \code{y} against each variable in the data frame\cr
#' \code{.~x}\tab Plot each variable in the data frame against \code{x}\cr
#' \code{.~.}\tab Plot each variable in the data frame against each other.\cr}
#'  See also section "Generating multiple plots at once" below.
#' @param weights observation weights or frequencies (optional).
#' @param opts a named list of options (optional); See also \code{\link{plotluck.options}}.
#' @param ... additional parameters to be passed to the respective ggplot2 geom objects.
#' @return a ggplot object, or a plotluck.multi object if the dot symbol was used.
#' @export
#' @keywords hplot
#' @keywords aplot
#' @keywords dplot
#' @concept automation
#' @concept visualization
#' @concept plotting
#' @concept exploratory data analysis
#' @concept ggplot
#' @concept ggplot2
#' @concept heat map
#' @concept density plot
#' @concept violin plot
#' @concept hexbin
#' @concept histogram
#' @concept bar plot
#' @concept box plot
#' @concept spine plot
#' @concept mosaic plot
#' @concept scatter plot
#' @concept heat map
#'
#' @seealso \code{\link{plotluck.options}}, \code{\link{sample.plotluck}}, \code{\link[ggplot2]{ggplot}}
#'
#' @section Determining the type of plot: Besides the shape of the formula, the
#'   algorithm takes into account the type of variables as either numeric, ordered,
#'   or unordered factors. Often, it makes sense to treat ordered factors similarly
#'   as numeric types.
#'
#'   One-variable numeric (resp. factor) distributions are usually represented by
#'   density (resp. Cleveland dot) charts, but can be overridden to histograms or
#'   bar plots using the \code{geom} option. Density plots come with an overlaid
#'   vertical median line.
#'
#'   For two numerical variables, by default a scatter plot is produced, but for
#'   high numbers of points a hexbin is preferred (option \code{min.points.hex}).
#'   These plots come with a smoothing line and standard deviation.
#'
#'   The relation between two factor variables can be depicted best by spine
#'   (a.k.a., mosaic) plots, unless they have too many levels (options
#'   \code{max.factor.levels.spine.x}, \code{max.factor.levels.spine.y},
#'   \code{max.factor.levels.spine.z}). Otherwise, a heat map is produced.
#'
#'   For a mixed-type (factor/numeric) pair of variables, violin (overridable
#'   to box) plots are generated. However, if the resulting graph would contain
#'   too many (more than \code{max.factor.levels.violin}) violin plots in a row,
#'   the algorithm switches automatically. The number of bins of a histogram can
#'   be customized with \code{n.breaks.histogram}. The default setting, \code{NA},
#'   applies a heuristic estimate.
#'
#'   The case of a response two dependent variables (`y~x+z`) is covered by
#'   either a spine plot (if all are factors) or a heat map.
#'
#'   In many cases with few points for one of the aggregate plots, a scatter
#'   looks better (options \code{min.points.density}, \code{min.points.violin},
#'   \code{min.points.hex}).
#'
#'   If each factor combination occurs only once in the data set, we resort to
#'   bar plots.
#
#' @section Conditional variables: Conditional variables are represented by either
#'   trying to fit into the same graph using coloring (\code{max.factor.levels.color}),
#'   or by facetting (preferred dimensions \code{facet.num.wrap} (resp.
#'   \code{facet.num.grid}) for one resp. two variables). Numeric vectors are
#'   discretized accordingly. Facets are laid out horizontally or vertically
#'   according to the plot type, up to maximum dimensions of \code{facet.max.rows}
#'   and \code{facet.max.cols}.
#'
#' @section Reordering of factor levels: To better illustrate the relation
#'  between an independent factor variable and a dependent numerical variable
#'  (or an ordered factor), levels are reordered according to the value of the
#'  dependent variable. If no other numeric or ordered variable exists, we
#'  sort by frequency.
#'
#' @section Instance weights: Argument \code{weights} allows to specify weights
#'  or frequency counts for each row of data. All plots and summary statistics
#'  take weights into account when supplied. In scatter and heat maps, weights
#'  are indicated either by a shaded disk with proportional area (default) or by
#'  jittering (option \code{dedupe.scatter}), if the number of duplicated points
#'  exceeds \code{min.points.jitter}. The amount of jittering can be controlled
#'  with \code{jitter.x} and \code{jitter.y}.
#'
#' @section Axis scaling: \code{plotluck} supports logarithmic and log-modulus
#'  axis scaling. log-modulus is considered if values are both positive and
#'  negative; in this case, the transform function is \code{f(x) = sign(x) *
#'  log(1+abs(x))}.
#'
#'  The heuristic to apply scaling is based on the proportion of total display
#'  range that is occupied by the 'core' region of the distribution between the
#'  lower and upper quartiles; namely, the fact whether the transform could
#'  magnify this region by a factor of at least \code{trans.log.thresh}.
#'
#' @section Missing values: By default, missing (\code{NA} or \code{NaN}) values
#'  in factors are are shown as a special factor level code{"?"}. They can be
#'  removed by setting \code{na.rm=TRUE}. Conventionally, missing numeric values
#'  are not shown.
#'
#' @section Sampling: For very large data sets, plots can take a very long time
#'  (or even crash R). \code{plotluck} has a built-in stop-gap: If the data
#'  comprises more than \code{sample.max.rows}, it will be sampled down to that
#'  size (taking into account \code{weights}, if supplied).
#'
#' @section Factor preprocessing:  Character (resp. logical) vectors are converted to
#'  unordered (resp. ordered) factors.
#'
#'  Frequently, when numeric variables have very few values despite sufficient
#'  data size, it helps to treat these values as the levels of a factor; this is
#'  governed by option \code{few.unique.as.factor}.
#'
#'  If an unordered factor has too many levels, plots can get messy. In this
#'  case, only the \code{max.factor.levels} most frequent ones are retained,
#'  while the rest are merged into a default level \code{".other."}.
#'
#' @section Coloring: If \code{color} or \code{fill} aesthetics are used to
#'  distinguish different levels or ranges of a variable, the color scheme adjusts
#'  to the type. Preferably, a sequential (resp. qualitative) palette is chosen
#'  for a numeric/ordered (unordered) factor (\code{palette.brewer.seq},
#'  \code{palette.brewer.qual}); see also \link[RColorBrewer]{RColorBrewer}.
#'
#' @section Generating multiple plots at once: If \code{formula} contains a dot
#'  (\code{"."}) symbol, the function creates a number of 1D or 2D plots by calling
#'  \code{plotluck} repeatedly. As described above, this allows either single
#'  distribution, one-vs-all and all-vs-all variable plots. To save space,
#'  rendering is minimal without axis labels.
#'
#'  In the all-vs-all case, the diagonal contains 1D distribution plots, analogous
#'  to the behavior of the default plot method for data frames, see
#'  \code{\link[graphics]{plot.data.frame}}.
#'
#'  With setting \code{in.grid=FALSE}, plots are produced in a sequence, otherwise
#'  together on one or multiple pages, if necessary (default). Page size is
#'  controlled by \code{multi.max.rows} and  \code{multi.max.cols}.
#'
#'  With \code{entropy.order=TRUE}, plots are sorted by an estimate of
#'  empirical conditional entropy, with the goal of prioritizing the more
#'  predictive variables. Set \code{verbose=TRUE} if you want to see the actual
#'  values. For large data sets the calculation can be time consuming; entropy
#'  calculation can be suppressed by setting \code{multi.entropy.order=FALSE}.
#'
#'  @note The return value is an object of class \code{plotluck_multi}. This
#'   class does not have any functionality; its sole purpose is to make this
#'   function work in the same way as \code{ggplot} and \code{plotluck}, namely,
#'   do the actual drawing if and only if the return value is not assigned.
#'
#' @section Debugging: With the option \code{verbose=TRUE} turned on, the function
#'  will print out information about the chosen and applicable plot types, ordering,
#'  log scaling, etc.
#'
#' @section Column name matching: Variable names can be abbreviated if they match
#'  a column name uniquely by prefix.
#'
#' @section Remarks on supported plot types: By default, \code{plotluck}
#'  uses violin and density plots in place of the more traditional box-and-whisker
#'  plots and histograms; these modern graph types convey the shape of a
#'  distribution better. In the former case, summary statistics like mean and
#'  quantiles are less useful if the distribution is not unimodal; a wrong
#'  choice of the number of bins of a histogram can create misleading artifacts.
#'
#'  Following Cleveland's advice, factors are plotted on the y-axis to make labels
#'  most readable and compact at the same time. This direction can be controlled
#'  using option \code{prefer.factors.vert}.
#'
#'  Due to their well-documented problematic aspects, pie charts and stacked bar
#'  graphs are not supported.
#'
#'  With real-world data (as opposed to smooth mathematical functions),
#'  three-dimensional scatter, surface, or contour plots can often be hard to
#'  read if the shape of the distribution is not suitable, data coverage is
#'  uneven, or if the perspective is not carefully chosen depending on the data.
#'  Since they usually require manual tweaking, we have refrained from
#'  incorporating them.
#'
#' @section Remarks on the use of options: For completeness, we have included the
#'  description of option parameters in the current help page. However, the
#'  tenet of this function is to be usable "out-of-the-box", with no or very
#'  little manual tweaking required. If you find yourself needing to change
#'  option values repeatedly or find the presets to be suboptimal, please
#'  contact the author.
#
#' @section Limitations: \code{plotluck} is designed for generic out-of-the-box
#'  plotting, and not suitable to produce more specialized types of plots that
#'  arise in specific application domains (e.g., association, stem-and-leaf,
#'  star plots, geographic maps, etc). It is restricted to at most three variables.
#'  Parallel plots with variables on different scales (such as time
#'  series of multiple related signals) are not supported.
#'
#' @examples
#' # Single-variable density
#' data(diamonds, package='ggplot2')
#' plotluck(diamonds, price~1)
#' invisible(readline(prompt="Press [enter] to continue"))
#'
#' # Violin plot
#' data(iris)
#' plotluck(iris, Species~Petal.Length)
#' invisible(readline(prompt="Press [enter] to continue"))
#'
#' # Scatter plot
#' data(mpg, package='ggplot2')
#' plotluck(mpg, cty~model)
#' invisible(readline(prompt="Press [enter] to continue"))
#'
#' # Spine plot
#' data(Titanic)
#' plotluck(as.data.frame(Titanic), Survived~Class+Sex, weights=Freq)
#' invisible(readline(prompt="Press [enter] to continue"))
#'
#' # Facetting
#' data(msleep, package='ggplot2')
#' plotluck(msleep, sleep_total~bodywt|vore)
#' invisible(readline(prompt="Press [enter] to continue"))
#'
#' # Heat map
#' plotluck(diamonds, price~cut+color)
#'
#'\donttest{
#' # Multi plots
#
#' # All 1D distributions
#' plotluck(iris, .~1)
#'
#' # 2D dependencies with one fixed variable on vertical axis
#' plotluck(iris, Species~.)
#'}
#' # See also tests/testthat/test_plotluck.R for more examples!
#'
plotluck <- function(data, formula, weights,
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

  # validation and preprocessing
  if (!('.' %in% vars)) {
     data <- data[,vars, drop=FALSE]
  }
  data <- sample.data(data, weights, opts$sample.max.rows)
  data <- process.weights(data, weights)
  data <- preprocess.factors(data, setdiff(vars, c(weights, '.')), opts$na.rm)

  info.threshold(opts$verbose && nrow(data) > opts$sample.max.rows,
                 'Data set has %s rows, sampling down', opts$sample.max.rows,
                 nrow(data))

  if ('.' %in% vars) {
    # trellis of plots
    return(plotluck.multi(response, indep, data, w=weights,
                          in.grid=opts$multi.in.grid,
                          max.rows=opts$multi.max.rows,
                          max.cols=opts$multi.max.cols,
                          opts=opts,
                          ...))
  }

  for (v in c(response, indep)) {
    if (is.numeric(data[[v]])) {
      data <- discretize.few.unique(data, v, opts$few.unique.as.factor, opts$verbose)
    }
  }

  # conditionals: discretize if numeric, and order
  # note: to have a single data frame containing all plot and facet variables,
  # we have to preprocess the conditionals beforehand.

  if (length(cond) > 0) {

    # if we have to discretize, we want that many facets
    intervals <- NULL

    if (length(cond) == 1) {
      intervals <- opts$facet.num.wrap
      info.threshold(opts$verbose && length(unique(data[[cond]])) > opts$facet.num.wrap, 'Discretizing %s into intervals', opts$facet.num.wrap, cond)
    } else {
      intervals <- opts$facet.num.grid
      info.threshold(opts$verbose && length(unique(data[[cond[1]]])) > opts$facet.num.grid, 'Discretizing %s into intervals', opts$facet.num.wrap, cond[1])
      info.threshold(opts$verbose && length(unique(data[[cond[2]]])) > opts$facet.num.grid, 'Discretizing %s into intervals', opts$facet.num.wrap, cond[2])
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
        data[[cond[i]]] <- ordered(data[[cond[i]]])
      } else {
        data <- limit.factor.levels(data, cond[i], w=weights,
                                    max.levels=min(intervals, opts$max.factor.levels))
        if (!is.ordered(data[[cond[i]]])) {
          data <- order.factor.by.value(data, cond[i], order.by, weights, verbose=opts$verbose)
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
    max.points.per.group <- max(t)
  } else {
    num.groups <- 1
    med.points.per.group <- nrow(data)
    max.points.per.group <- nrow(data)
  }

  # precompute medians
  if (length(vars.numeric) == 1) {
     grp.med <- group.central(data, vars.numeric, vars.non.numeric, w=weights, method='median')
     if (length(vars.non.numeric) == 0) {
        # single value
        data$.center. <- grp.med$.center.
     } else {
      data <- merge(data, grp.med)
     }
  }

  # note: factor levels cannot be limited before ordering is known!

  # determine type of plot
  if (length(indep) == 0) {
    # distribution plot
    res <- determine.plot.type.0(data, response, weights,
                                 med.points.per.group, num.groups, opts, ...)
  } else if (length(indep) == 1) {
    res <- determine.plot.type.1(data, response, indep, cond, weights,
                                 med.points.per.group, max.points.per.group,
                                 num.groups, opts, ...)
  } else { # length(indep) == 2
    res <- determine.plot.type.2(data, response, indep[[1]], indep[[2]],
                                 weights, opts, ...)
  }

  p <- res$plot
  data <- res$data
  type.plot <- res$type.plot

  if (length(indep) < 2) {
     # data is possibly modified, we need to pass it in
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
  res <- NULL
  type.plot <- NULL
  opts.geom <- NULL
  choices.geom <- NULL
  if (is.numeric(data[[response]])) {
    choices.geom <- c('density', 'histogram', 'scatter', 'auto')
    opts.geom <- match.arg(opts$geom, choices.geom)
    if (opts.geom == 'auto') {
      # enough points for density plot?
      opts.geom <- ifelse(med.points.per.group > opts$min.points.density, 'density', 'scatter')
    }
    if (opts.geom == 'density') {
      type.plot <- 'density'
      res <- gplt.density(data, response, w=w,
                        trans.log.thresh=opts$trans.log.thresh,
                        verbose=opts$verbose,
                        ...)
    } else if (opts.geom == 'histogram') {
      type.plot <- 'histogram'
      res <- gplt.histogram(data, response, w=w,
                          n.breaks=opts$n.breaks.histogram,
                          trans.log.thresh=opts$trans.log.thresh,
                          verbose=opts$verbose,
                          ...)
    } else { # opts.geom == 'scatter'

      # to plot a distribtion as a line plots,
      # make up a constant y coordinate
      indep <- '.y_const'
      data[[indep]] <- 1
      data[[indep]] <- as.factor(data[[indep]])
      type.plot <- 'scatter.num.1'
      if (!opts$prefer.factors.vert) {
         res <- gplt.scatter(data, indep, response, w=w,
                             dedupe.scatter=opts$dedupe.scatter,
                             min.points.jitter=opts$min.points.jitter,
                             trans.log.thresh=opts$trans.log.thresh,
                             max.factor.levels=opts$max.factor.levels,
                             verbose=opts$verbose,
                             ...)
         res$plot <- res$plot +
            theme(axis.ticks.x=element_blank(),
                  axis.text.x=element_blank()) +
            xlab('')
         # known issue: the conditioning layer can add theme_slanted_text_x,
         # in which case spurious '1's appear on the x-axis.
      } else {
         res <- gplt.scatter(data, response, indep, w=w,
                             dedupe.scatter=opts$dedupe.scatter,
                             min.points.jitter=opts$min.points.jitter,
                             trans.log.thresh=opts$trans.log.thresh,
                             max.factor.levels=opts$max.factor.levels,
                             verbose=opts$verbose,
                             ...)
         res$plot <- res$plot +
            theme(axis.ticks.y=element_blank(),
                  axis.text.y=element_blank()) +
            xlab('')
      }
    }
  } else { # !is.numeric(data[[response]])
    choices.geom <- c('dot', 'bar', 'auto')
    opts.geom <- match.arg(opts$geom, choices.geom)
    if (opts.geom == 'auto') {
      opts.geom <- ifelse(num.groups < 6, 'bar', 'dot')
    }
    type.plot <- opts.geom
    res <- gplt.dot(data, response, w=w, vertical=opts$prefer.factors.vert,
                  max.factor.levels=opts$max.factor.levels, geom=opts.geom,
                  verbose=opts$verbose, ...)
  }
  info.options(opts.geom, choices.geom, opts$verbose)
  res$type.plot=type.plot
  res
}

determine.plot.type.1 <- function(data, response, indep, cond,
                                  w, med.points.per.group,
                                  max.points.per.group,
                                  num.groups, opts,
                                  ...) {
  res <- NULL
  type.plot <- NULL
  opts.geom <- NULL
  choices.geom <- NULL
  if (is.numeric(data[[indep]]) && is.numeric(data[[response]])) {
    choices.geom <- c('hex', 'scatter', 'auto')
    opts.geom <- match.arg(opts$geom, choices.geom)
    if (opts.geom == 'auto') {
      # if a lot of points, choose hex
      opts.geom <- ifelse(nrow(data) >= opts$min.points.hex, 'hex', 'scatter')
    }
    # HACK: the smoothing line can only be colored with a default color if it
    # wouldn't be used to distinguish multiple groups.
    fill.smooth <- ifelse(length(cond) == 0, opts$fill.default, 'NULL')
    if (opts.geom == 'hex') {
      type.plot <- 'hex'
      res <- gplt.hex(data, indep, response, w=w, trans.log.thresh=opts$trans.log.thresh,
                    # hex plot already has color, make the smoothing line neutral
                    fill.smooth='grey', ...)
    } else {
      type.plot <- 'scatter.num'
      res <- gplt.scatter(data, indep, response, w,
                        dedupe.scatter=opts$dedupe.scatter,
                        min.points.jitter=opts$min.points.jitter,
                        jitter.x=opts$jitter.x,
                        jitter.y=opts$jitter.y,
                        trans.log.thresh=opts$trans.log.thresh,
                        max.factor.levels=opts$max.factor.levels,
                        fill.smooth=fill.smooth,
                        verbose=opts$verbose,
                        ...)
    }
  } else if (!is.numeric(data[[indep]]) && !is.numeric(data[[response]])) {
    choices.geom <-  c('spine', 'heat', 'auto')
    opts.geom <- match.arg(opts$geom, choices.geom)
    if (opts.geom == 'auto') {
      opts.geom <- ifelse(length(cond) == 0 && # HACK: spine plot cannot be faceted, due to difficulty of implementation
                            length(levels(data[[response]])) <= opts$max.factor.levels.spine.y &&
                            length(levels(data[[indep]])) <= opts$max.factor.levels.spine.x, 'spine', 'heat')
    }

    if (opts.geom == 'spine') {
       type.plot <- 'spine'
       res <- gplt.spine(data, indep, response, w=w,
                       plot.margin.x=opts$spine.plot.margin.x,
                       max.factor.levels=opts$max.factor.levels.spine.x,
                       palette.brewer.seq=opts$palette.brewer.seq,
                       palette.brewer.qual=opts$palette.brewer.qual,
                       verbose=opts$verbose,
                       ...)
    } else {
       type.plot <- 'heat'
       res <- gplt.heat(data, indep, response, z='NULL', w=w,
                        trans.log.thresh=opts$trans.log.thresh,
                        max.factor.levels=opts$max.factor.levels,
                        resolution.heat=opts$resolution.heat,
                        palette.brewer.seq=opts$palette.brewer.seq,
                        palette.brewer.qual=opts$palette.brewer.qual,
                        verbose=opts$verbose,
                        ...)
    }
  } else {
     # mixed types: one factor, one numeric

     if (max.points.per.group <= 1) {
        # special case: single point each - use bar plot
        choices.geom <- c('dot', 'bar', 'auto')
        opts.geom <- match.arg(opts$geom, choices.geom)
        if (opts.geom == 'auto') {
           opts.geom <- ifelse(num.groups < 6, 'bar', 'dot')
        }
        type.plot <- opts.geom
        var.fac <- response
        var.num <- indep
        if (is.numeric(data[[response]])) {
           var.fac <- indep
           var.num <- response
        }
        res <- gplt.dot(data, var.fac, w=var.num, vertical=opts$prefer.factors.vert,
                      max.factor.levels=opts$max.factor.levels, geom=opts.geom,
                      verbose=opts$verbose, ...)
     } else { #

        choices.geom <- c('violin', 'box', 'scatter', 'auto')
        opts.geom <- match.arg(opts$geom, choices.geom)
        var.fac <- response
        var.num <- indep
        # assignment of variables to axes depends only on opts$prefer.factors.vert
        switch <- is.numeric(data[[response]]) != opts$prefer.factors.vert
        if (switch) {
           var.num <- response
           var.fac <- indep
        }
        if (opts.geom == 'auto') {
           # if there are too many violin plots in a horizontal row, they
           # just look like black lines
           opts.geom <- ifelse(med.points.per.group > opts$min.points.violin,
                               ifelse(length(levels(var.fac)) <= opts$max.factor.levels.violin, 'violin', 'box'),
                               'scatter')
        }

        if (opts.geom == 'violin') {
           type.plot <- 'violin'
           res <- gplt.violin(data, var.fac, var.num, w=w,
                            trans.log.thresh=opts$trans.log.thresh,
                            max.factor.levels=opts$max.factor.levels,
                            verbose=opts$verbose,
                            ...)
        } else if (opts.geom == 'box') {
           type.plot <- 'box'
           res <- gplt.box(data, var.fac, var.num, w=w,
                         trans.log.thresh=opts$trans.log.thresh,
                         max.factor.levels=opts$max.factor.levels,
                         verbose=opts$verbose,
                         ...)
        } else { # opts.geom == 'scatter'
           # not enough points for violin or box plot
           type.plot <- 'scatter.mixed'
           res <- gplt.scatter(data, var.fac, var.num, w=w,
                             dedupe.scatter=opts$dedupe.scatter,
                             min.points.jitter=opts$min.points.jitter,
                             trans.log.thresh=opts$trans.log.thresh,
                             max.factor.levels=opts$max.factor.levels,
                             verbose=opts$verbose,
                             ...)
        }
     }
  }
  info.options(opts.geom, choices.geom, opts$verbose)
  res$type.plot=type.plot
  res
}

determine.plot.type.2 <- function(data, response, indep1, indep2, w, opts,
                                  ...) {
  res <- NULL
  type.plot <- NULL
  choices.geom <- c('spine', 'heat', 'auto')
  opts.geom <- match.arg(opts$geom, choices.geom)
  if (opts.geom == 'auto') {

    opts.geom <- ifelse((!is.numeric(data[[response]])) && (!is.numeric(data[[indep1]])) &&
                          (!is.numeric(data[[indep2]])) &&
                          length(levels(data[[response]])) <= opts$max.factor.levels.spine.z &&
                          length(levels(data[[indep2]])) <= opts$max.factor.levels.spine.y &&
                          length(levels(data[[indep1]])) <= opts$max.factor.levels.spine.x, 'spine', 'heat')
  }

  if (opts.geom == 'spine') {
    type.plot <- 'spine3'
    res <- gplt.spine3(data, indep1, indep2, response, w=w,
                     plot.margin.x=opts$spine.plot.margin.x,
                     plot.margin.y=opts$spine.plot.margin.y,
                     max.factor.levels=opts$max.factor.levels.spine.x,
                     palette.brewer.seq=opts$palette.brewer.seq,
                     palette.brewer.qual=opts$palette.brewer.qual,
                     verbose=opts$verbose,
                     ...)
  } else {
    type.plot <- 'heat3'
    res <- gplt.heat(data, indep1, indep2, response, w=w,
                     trans.log.thresh=opts$trans.log.thresh,
                     max.factor.levels=opts$max.factor.levels,
                     resolution.heat=opts$resolution.heat,
                     palette.brewer.seq=opts$palette.brewer.seq,
                     palette.brewer.qual=opts$palette.brewer.qual,
                     verbose=opts$verbose,
                     ...)
  }
  info.options(opts.geom, choices.geom, opts$verbose)
  res$type.plot=type.plot
  res
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
  show.var <- 'first'
  if (preferred.order == 'row') {
    ncol <- opts$facet.max.cols
  } else if (preferred.order == 'col') {
    nrow <- opts$facet.max.rows
    show.var <- 'last'
  }
  facet.labels <- format.facets(data, cond, show.var)

  p <- p + theme_slanted_text_x +  # axis text might overlap in small diagrams
    facet_wrap(cond,
               labeller=as_labeller(facet.labels),
               nrow=nrow, ncol=ncol) +
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
  } else if (!(type.plot %in% c('hex', 'scatter.num', 'heat'))) { # pure numeric -> no preference
    # mixed num/factor -> opposite to graph layout
    preferred.order <- ifelse(opts$prefer.factors.vert, 'row', 'col')
  }

  if (length(cond) == 1) {
    u <- length(unique(data[[cond]]))

    if (u <= opts$max.factor.levels.color && !(type.plot %in% c('spine', 'heat', 'hex'))) {
      p <- add.color.fill(p, data, cond,
                          palette.brewer.seq=opts$palette.brewer.seq,
                          palette.brewer.qual=opts$palette.brewer.qual)

      aesth.legend <- ifelse(type.plot %in% c('density', 'histogram', 'bar', 'violin'), 'fill', 'color')
      p <- add.color.legend(p, data, cond, aesth.legend)

    } else {
      p <- redundant.factor.color(p, data, response, indep, type.plot, opts)
      add.facet.wrap(p, data, cond, preferred.order, opts)
    }
  } else { #length(cond) == 2
    #u1 <- length(unique(data[[cond[1]]]))
    #u2 <- length(unique(data[[cond[2]]]))
    # if (min(u1,u2) <= opts$max.factor.levels.color) {
    #    if (u1 < u2) {
    #       p + aes_string(fill=cond[1]) + aes_string(color=cond[1]) +
    #          add.facet.wrap(p, data, cond[2], preferred.order, opts)
    #    } else {
    #       p + aes_string(fill=cond[2]) + aes_string(color=cond[2]) +
    #          add.facet.wrap(p, data, cond[1], preferred.order, opts)
    #    }
    # } else {
    # example: plotluck(diamonds, price~1|cut+clarity)
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

plotluck.multi <- function(response, indep, data, w='NULL',
                           in.grid=TRUE,
                           max.rows=10, max.cols=10,
                           opts=plotluck.options(),
                           ...) {

  main <- deparse(substitute(data))

  vars.x <- NULL
  vrs.y <- NULL
  if (response == '.') {
    vars.y <- names(data)
    if (length(indep) == 0) {
      vars.x <- '1'
      if (opts$verbose) {
        cat('Plotting distributions for all variables\n')
      }
    } else if (indep == '.') {
      vars.x <- names(data)
      if (opts$verbose) {
        cat('Plotting each variable against each other\n')
      }
    }  else {
      vars.x <- indep
      if (opts$verbose) {
        cat(sprintf('Plotting each variable against %s\n', vars.x))
      }
    }
  } else { # response != '.'
    vars.y <- response
    if(indep != '.') {
      stop('Error: invalid formula')
    }
    if (opts$verbose) {
      cat(sprintf('Plotting %s against each variable\n', vars.y))
    }
    vars.x <- names(data)
  }

  # order variables by conditional entropy
  cond.ent <- NULL
  if (opts$multi.entropy.order) {
    if (length(vars.y) > 1 && length(vars.x) == 1 && vars.x != '1') {
      cond.ent <- suppressWarnings(cond.entropy.data(data, given=vars.x, w=w))
      cond.ent <- cond.ent[,vars.x]
      vars.y <- vars.y[order(cond.ent)]
    } else if (length(vars.y) == 1 && length(vars.x) > 1) {
      cond.ent <- suppressWarnings(cond.entropy.data(data, target=vars.y, w=w))
      cond.ent <- cond.ent[vars.y,]
      vars.x <- vars.x[order(cond.ent)]
    } else if (length(vars.y) > 1 && length(vars.x) > 1) {
      cond.ent <- suppressWarnings(cond.entropy.data(data, w=w))
      cond.ent <- (apply(cond.ent, 1, mean) + apply(cond.ent, 1, mean))/2
      vars.x <- vars.x[order(cond.ent)]
      vars.y <- vars.x
    }
    if (opts$verbose && !is.null(cond.ent)) {
      cat('Ordering variables according to conditional entropy:\n')
      cond.ent<-sort(cond.ent)
      out<-data.frame(var=names(cond.ent),cond.ent=cond.ent)
      print(out, row.names=FALSE)
    }
  }

  combi <- expand.grid(seq(length(vars.x)),
                       seq(length(vars.y)), stringsAsFactors=FALSE)
  names(combi) = c('grid.x', 'grid.y')
  combi$x <- vars.x[combi$grid.x]
  combi$y <- vars.y[combi$grid.y]
  combi$fac.x <- sapply(combi$x, function(v) !is.numeric(data[[v]]))
  combi$fac.y <- sapply(combi$y, function(v) !is.numeric(data[[v]]))

  # try to make a square layout
  cols <- ceiling(sqrt(nrow(combi)))
  rows <- ceiling(nrow(combi)/cols)

  # figure out which plots are not at the margins, and
  # don't show axis labels for them if redundant
  suppress.xlab <- FALSE
  suppress.ylab <- FALSE

  theme.multi <- ''

  is.square <- TRUE

  if (!in.grid && opts$verbose) {
    cat('Not plotting all graphs on one page because multi.in.grid=FALSE\n')
  }

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
    opts$dedupe.scatter <- 'jitter'

    # use the limited space to make subgraphs
    # relatively rectangular
    opts$facet.max.rows <- NULL
    opts$facet.max.cols <- NULL

    # do all plots fit on a single page?
    if (cols > max.cols || cols > max.rows) {
      is.square <- FALSE
      cols <- max.cols
      rows <- min(max.rows, ceiling(nrow(combi)/cols))
      if (opts$verbose) {
        cat(sprintf('Cannot fit all %d variables on one page; multi.max.rows = %d, multi.max.cols= %d\n',
                    max.rows, max.cols))
      }
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
        main <- vars.y
        suppress.ylab <- 'margin'
      } else {
        suppress.ylab <- 'all'
      }
    }

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
  }

  combi$opts <- sprintf('prefer.factors.vert=%s', opts$prefer.factors.vert)
  combi$labs <- ''
  if (length(vars.x) > 1 && length(vars.y) > 1 &&
      ((!in.grid) || (!is.square))) {
    # plots with the axis transposed are redundant;
    # if not part of grid, no need to show them
    combi <- combi[combi$grid.x >= combi$grid.y,]

  } else {

     # For mixed factor/num combinations, plotluck() assigns the axes only
     # depending on opts$prefer.factors.vert. In the output lattice for
     # plotluck.multi, we need to control the direction (and show opposites
     # at mirror position)
     combi$opts <- 'prefer.factors.vert=TRUE'
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
  }

  # plot single-variable distributions for x~x
  combi$x[combi$x == combi$y] <- '1'

  if (w != 'NULL')
  {
    call.strs <- sprintf('plotluck(data, %s~%s, w=%s, opts=plotluck.options(opts,%s), ...)%s%s',
                         combi$y, combi$x, w, combi$opts, combi$labs, theme.multi)
  } else {
    call.strs <- sprintf('plotluck(data, %s~%s, opts=plotluck.options(opts,%s), ...)%s%s',
                         combi$y, combi$x, combi$opts, combi$labs, theme.multi)
  }

  call.str <- sprintf('list(%s)', paste(call.strs, collapse=','))

  # don't show verbose messages for all individual plots
  opts$verbose <- FALSE

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
#' @export
print.plotluck_multi <- function(x, ...) {
  if (x$in.grid) {
    mplot(x$plots, rows=x$rows, cols=x$cols,
          suppress.xlab=x$suppress.xlab, suppress.ylab=x$suppress.ylab)
  } else {
    lapply(x$plots, print)
  }
}


#' Run \code{plotluck} for a randomly generated formula.
#'
#' \code{sample.plotluck} samples a formula as follows:
#' \itemize{
#' \item Uniformly draw the number of variables (1-3).
#' \item For each variable, uniformly choose one of the existing variable types from the data set (numeric, ordered or unordered factor).
#' \item Uniformly select one of the data frame columns of that type.
#'}
#'
#' @param data a data frame
#' @param ... additional parameters to be passed to \code{plotluck}, such as
#'       \code{weights} and \code{opts}.
#' @return a ggplot2 object.
#'
#' @seealso \code{\link{plotluck}}
#' @export
#' @examples
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
  if (length(arg) == 1) {
     s <- sprintf('plotluck(%s, %s)\n', as.character(match.call()[2]), form)
  } else {
     s <- list()
     for (i in seq(2, length(arg))) {
        s[[length(s)+1]] <- sprintf('%s = %s', names(arg)[i], arg[i])
     }
     s <- do.call(paste, c(s, sep=', '))
     s <- sprintf('plotluck(%s, %s, %s)\n', as.character(match.call()[2]), form, s)
  }
  cat(s)
  plotluck(data, as.formula(form), ...) + labs(title=s) + theme(plot.title=element_text(size=10))
}

# same as sample.plotluck, but can be called with a list of options
# only used for testing/debugging!

# e.g.
# opts.list<-list()
# opts.list[[1]]<-plotluck.options(verbose=T)
# opts.list[[2]]<-plotluck.options(verbose=T,prefer.factors.vert=F)
# opts.list[[3]]<-plotluck.options(verbose=T,prefer.factors.vert=F,max.factor.levels.color=100)
# opts.list[[4]]<-plotluck.options(verbose=T,prefer.factors.vert=T,max.factor.levels.color=100,dedupe.scatter='jitter',min.points.hex=10000,min.points.density=1E20,min.points.violin=1E20)
# opts.list[[5]]<-plotluck.options(verbose=T,prefer.factors.vert=F,max.factor.levels=3)


sample.plotluck.testopts <- function(data, opts.list, ...) {
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
   if (length(arg) <= 2) {
      s <- sprintf('plotluck(%s, %s)', as.character(match.call()[2]), form)
   } else {
      s <- list()
      for (i in seq(3, length(arg))) {
         s[[length(s)+1]] <- sprintf('%s = %s', names(arg)[i], arg[i])
      }
      s <- do.call(paste, c(s, sep=', '))
      s <- sprintf('plotluck(%s, %s, %s)', as.character(match.call()[2]), form, s)
   }

   for (i in seq(length(opts.list))) {
      s.opt <- paste(s, sprintf(' [OPTION %d]', i), sep='')
      cat(s.opt)
      print(plotluck(data, as.formula(form), opts=opts.list[[i]],...) + labs(title=s.opt) + theme(plot.title=element_text(size=10)))
   }
}
