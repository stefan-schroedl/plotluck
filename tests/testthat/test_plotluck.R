
# convenience functions

.file.cnt <- 0
.prefix <- 'test'

# set .test.mode to 'demo' to visually inspect all tests
# set .test.mode to 'regression' to ensure plots stay identical

if (exists('.test.mode') && .test.mode == 'regression') {
   context("Regression Tests")

   test_that_ref <- function(prefix, desc, code) {
      .prefix <<- prefix
      .file.cnt <<- 0
      test_that(desc, code)
   }
   eq_ref <- function(x) {
      filename <- sprintf('%s_%d.rds', .prefix, .file.cnt)
      .file.cnt <<- .file.cnt + 1
      expect_equal_to_reference(x, filename)
   }

} else {
   test_that_ref <- function(prefix, desc, code) {
      cat(sprintf('*** %s ***\n', desc))
      code
   }

   if (exists('.test.mode') && .test.mode == 'demo') {
      eq_ref <- function(x) {
         print(x)
         l <- readline(prompt="Hit <RETURN> to continue, anything else to quit: ")
         if (nchar(l) > 0) {
            stop('', call.=FALSE, domain=NA)
         }
      }
   } else {
      eq_ref <- function(x)  { print(x); expect_true(TRUE) }
   }
}


# load data
data(iris)
if (!requireNamespace('ggplot2movies', quietly=TRUE)) {
   install.packages('ggplot2movies')
}
data(movies, package='ggplot2movies')

data(Titanic)
data(occupationalStatus)
data(diamonds, package='ggplot2')

# 1D

test_that_ref("1d_density", "1D density plot", {
   eq_ref(plotluck(iris, Petal.Length~1))
})

test_that_ref("1d_scatter", "1D scatter num/fact", {
   # area
   eq_ref(plotluck(iris, Petal.Length~1, opts=plotluck.options(min.points.density=1E20)))
   # jitter
   eq_ref(plotluck(iris, Petal.Length~1, opts=plotluck.options(min.points.density=1E20, dedupe.scatter='jitter')))
})


test_that_ref("1d_bar", "1D bar", {
   eq_ref(plotluck(movies, mpaa~1))
})


test_that_ref("1d_scaling", "log scaling", {
   set.seed(0)
   n <- 1000
   m <- 100
   df <- data.frame(a=rnorm(n, mean=m, sd=5),
                    b=c(10*m, rnorm(n-1, mean=m, sd=5)),
                    c=m*c(-2*m, rexp(n-1, 1)),
                    d=-m*rexp(n, 1))
   eq_ref(plotluck(df, a~1))
   eq_ref(plotluck(df, b~1))
   eq_ref(plotluck(df, c~1))
   eq_ref(plotluck(df, d~1))
})


test_that_ref("2d_scatter", "2D scatter", {
   eq_ref(plotluck(iris, Petal.Length~Petal.Width))

   # hex
   # note: when running devtools::check(), weird dependency problem with mgce,
   #  nlme; comment out the following three tests
   eq_ref(plotluck(movies, votes~rating))

   # scatter
   eq_ref(plotluck(movies, rating~votes,
                  opts=plotluck.options(min.points.hex=1E20)))

   # scatter, no log scaling
   eq_ref(plotluck(movies, rating~votes,
                   opts=plotluck.options(min.points.hex=1E20, trans.log.thresh=100)))
})


test_that_ref("2d_density", "2D density", {
   # using facets
   eq_ref(plotluck(movies, rating~1|mpaa))

   # using colors
   eq_ref(plotluck(iris, Petal.Length~1|Species))

   # scatter
   eq_ref(plotluck(iris, Petal.Length~Species,
                   opts=plotluck.options(min.points.violin=1E20)))
   # ordering
   i2 <- iris
   i2$Species <- as.ordered(i2$Species)
   eq_ref(plotluck(i2, Petal.Length~1|Species))   # note color scheme
   eq_ref(plotluck(diamonds, price~cut))

})


test_that_ref("2d_box", "2D box/violin plot", {
   eq_ref(plotluck(movies, rating~mpaa))

   # box instead of violin
   eq_ref(plotluck(movies, rating~mpaa, opts=plotluck.options(geom='box')))

   # as scatter plot
   eq_ref(plotluck(movies, rating~mpaa, opts=plotluck.options(min.points.violin=1E20)))

   # with jittering
   eq_ref(plotluck(movies, rating~mpaa, opts=plotluck.options(min.points.violin=1E20, dedupe.scatter='jitter')))

   # implicit conversion of binary variables
   eq_ref(plotluck(movies, budget~Documentary))

})


test_that_ref("2d_id", "2D identity bar", {
   # identity bar
   df <- data.frame(f=factor(c('aaaaaaa','bbbbbbbbb','ccccccccc','dddddddd')), val=c(5,6,2,8))
   eq_ref(plotluck(df, val~f))
})


test_that_ref("2d_spine", "2D spine", {

   # spine plot
   eq_ref(plotluck(as.data.frame(Titanic), Survived~Class, weights=Freq))

   df <- as.data.frame(occupationalStatus)
   df$origin <- ordered(df$origin)
   df$destination <- ordered(df$destination)
   eq_ref(plotluck(df, destination~origin, weights=Freq))
})

test_that_ref("3d_identity", "3D identity", {
   df <- data.frame(f=factor(c('aaaaaaa','bbbbbbbbb','ccccccccc','dddddddd')),
                    f2=factor(c(1,1,2,2)), val=c(5,6,2,8))
   eq_ref(plotluck(df, val~f|f2))

   eq_ref(plotluck(df, val~f|f2, opts=plotluck.options(max.factor.levels.color=0)))
})

test_that_ref("3d_heat", "3D heat map", {
   eq_ref(plotluck(diamonds, price~cut+color))
})


test_that_ref("3d_spine", "3D spine", {
   eq_ref(plotluck(as.data.frame(Titanic), Survived~Class+Sex, weights=Freq))
})


test_that_ref("3d_scatter", "3D scatter", {
   # using hex, facets
   eq_ref(plotluck(movies, rating~length|mpaa))

   # using scatter, facets
   eq_ref(plotluck(movies, rating~length|mpaa,
                   opts=plotluck.options(min.points.hex=1E20)))

   # using colors
   eq_ref(plotluck(iris, Petal.Width~Petal.Length|Species))
})



test_that_ref("3d_density", "3D density", {
   # density, facets
   eq_ref(plotluck(diamonds, price~1|cut+color))

   # scatter, facets
   eq_ref(plotluck(diamonds, price~1|cut+color,
                   opts=plotluck.options(min.points.density=1E20)))

   # colors
   eq_ref(plotluck(diamonds, price~1|cut,
                   opts=plotluck.options(max.factor.levels.color=1E20)))
})



test_that_ref("3d_violin", "3D violin", {
   # color
   eq_ref(plotluck(movies, rating~mpaa|Action))

   # facets
   eq_ref(plotluck(movies, rating~mpaa|Action,
                   opts=plotluck.options(max.factor.levels.color=0)))

   # scatter, color
   eq_ref(plotluck(movies, rating~mpaa|Action,
                   opts=plotluck.options(min.points.violin=1E20)))

   # scatter, color, jitter
   eq_ref(plotluck(movies, rating~mpaa|Action,
                   opts=plotluck.options(min.points.violin=1E20, dedupe.scatter='jitter')))

   # scatter, facets
   eq_ref(plotluck(movies, rating~mpaa|Action,
                   opts=plotluck.options(min.points.violin=1E20, max.factor.levels.color=0)))
})


test_that_ref("3d_spine", "3D spine", {
   eq_ref(plotluck(as.data.frame(Titanic), Survived~Class+Sex, weights=Freq))
})


test_that_ref("missing", "missing values", {
   set.seed(0)
   df<-data.frame(f=factor(sample(c(letters[1:3], NA), 20, replace=TRUE)),
                  f2=factor(sample(c(1,2,3,NA), 20, replace=TRUE)),
                  v=runif(20),
                  v2=runif(20))
   df$v2[c(1,5,6)] <- NA
   eq_ref(plotluck(df, v~f))
   eq_ref(plotluck(df, v~f, opts=plotluck.options(na.rm=TRUE)))
   eq_ref(plotluck(df, v~f, opts=plotluck.options(min.points.violin=0)))
   eq_ref(plotluck(df, v~f, opts=plotluck.options(min.points.violin=0, na.rm=TRUE)))
   eq_ref(plotluck(df, f~v))
   eq_ref(plotluck(df, v2~v|f))
   eq_ref(plotluck(df, f2~f))
   eq_ref(plotluck(df, f2~f, opts=plotluck.options(na.rm=TRUE)))
})

test_that_ref("2d_weight", "instance weights", {
   set.seed(0)
   df<-data.frame(f=factor(sample(letters[1:3], 20, replace=TRUE), exclude=FALSE),
                  v1=runif(20),
                  v2=runif(20))
   df$w<-runif(20) + ifelse(df$v1>0.6, 3 * runif(20), ifelse(df$v1>0.3, 2*runif(20), 0))

   # num/num scatter
   eq_ref(plotluck(df, v1~v2))
   eq_ref(plotluck(df, v1~v2, weights=w))

   # violin
   eq_ref(plotluck(df, v1~f, opts=plotluck.options(geom='violin')))
   eq_ref(plotluck(df, v1~f, weights=w, opts=plotluck.options(geom='violin')))

   # box
   eq_ref(plotluck(df, v1~f, opts=plotluck.options(geom='box')))
   eq_ref(plotluck(df, v1~f, weights=w, opts=plotluck.options(geom='box')))

   # density
   eq_ref(plotluck(df, v1~1, opts=plotluck.options(min.points.density=0)))
   eq_ref(plotluck(df, v1~1, weights=w, opts=plotluck.options(min.points.density=0)))

   # histogram
   eq_ref(plotluck(df, v1~1,
                   opts=plotluck.options(geom='histogram')))
   eq_ref(plotluck(df, v1~1, weights=w,
                   opts=plotluck.options(geom='histogram')))

   # heat map/spine
   df<-expand.grid(1:5, 1:5)
   df <- rbind(df, df, df)
   df$v <- runif(75)
   df$w<-runif(75) + ifelse(df$v>0.6, 10 * runif(75), ifelse(df$v>0.3, 5*runif(75), 0))
   eq_ref(plotluck(df, v~Var1+Var2, opts=plotluck.options(geom='spine')))
   eq_ref(plotluck(df, v~Var1+Var2, weights=w, opts=plotluck.options(geom='spine')))

})

test_that_ref("multi", "multiple plots", {

   #testthat::skip_on_cran()
   if (identical(Sys.getenv("NOT_CRAN"), "true")) {
      eq_ref(plotluck(diamonds, .~1))
      eq_ref(plotluck(diamonds, price~.))
      eq_ref(plotluck(diamonds, .~price))
   }

})
