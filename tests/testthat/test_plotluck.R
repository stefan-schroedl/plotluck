
# convenience functions

.file.cnt <- 0
.prefix <- 'test'

# set .debug to TRUE to visually inspect all tests

if (!exists('.debug') || !.debug) {
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

   eq_ref <- function(x) {
      print(x)
      l <- readline(prompt="Hit <RETURN> to continue, anything else to quit: ")
      if (nchar(l) > 0) {
         stop('', call.=FALSE, domain=NA)
      }
   }
}


# load data
data(iris)
data(movies, package='ggplot2')
data(Titanic)
data(occupationalStatus)
data(diamonds, package='ggplot2')

# 1D

test_that_ref("1d_density", "1D density plot", {
   eq_ref(plotluck(iris, Petal.Length))
})

test_that_ref("1d_scatter", "1D scatter num/fact", {
   eq_ref(plotluck(iris, Petal.Length, opts=plotluck.options(min.points.density=1E20)))
   # with jitter
   eq_ref(plotluck(iris, Petal.Length, opts=plotluck.options(min.points.density=1E20, convert.duplicates.to.weights=FALSE)))
})


test_that_ref("1d_bar", "1D bar", {
   eq_ref(plotluck(movies, mpaa))
})


test_that_ref("1d_scaling", "scaling", {
   set.seed(0)
   n <- 1000
   m <- 100
   df <- data.frame(a=rnorm(n, mean = m, sd = 5),
                    b=c(10*m, rnorm(n-1, mean = m, sd = 5)),
                    c=m*c(-2*m, rexp(n-1, 1)),
                    d=-m*rexp(n, 1))
   eq_ref(plotluck(df, a))
   eq_ref(plotluck(df, b))
   eq_ref(plotluck(df, c))
   eq_ref(plotluck(df, d))
})


test_that_ref("2d_scatter", "2D scatter", {
   eq_ref(plotluck(iris, Petal.Length, Petal.Width))

   # hex
   eq_ref(plotluck(movies, votes, rating))

   # scatter
   eq_ref(plotluck(movies, votes, rating,
                   opts=plotluck.options(min.points.hex = 1E20)))

   # scatter, no log scaling
   eq_ref(plotluck(movies, votes, rating,
                   opts=plotluck.options(min.points.hex = 1E20, trans.log.thresh = 0)))
})


test_that_ref("2d_density", "2D density", {
   # using facets
   eq_ref(plotluck(movies, rating, mpaa))

   # using colors
   eq_ref(plotluck(iris, Petal.Length, Species))

   # scatter
   eq_ref(plotluck(iris, Petal.Length, Species,
                   opts=plotluck.options(min.points.density = 1E20)))
   # ordering
   i2 <- iris
   i2$Species <- as.ordered(i2$Species)
   eq_ref(plotluck(i2, Petal.Length, Species))   # scatter
   eq_ref(plotluck(diamonds, price, cut))

})


test_that_ref("2d_box", "2D box/violin plot", {
   eq_ref(plotluck(movies, mpaa, rating))

   # box instead of violin
   eq_ref(plotluck(movies, mpaa, rating, opts=plotluck.options(use.geom.violin=FALSE)))

   # as scatter plot
   eq_ref(plotluck(movies, mpaa, rating, opts=plotluck.options(min.points.box = 1E20)))

   # with jittering
   eq_ref(plotluck(movies, mpaa, rating, opts=plotluck.options(min.points.box = 1E20, convert.duplicates.to.weights=FALSE)))

   # implicit conversion of binary variables
   eq_ref(plotluck(movies, Documentary, budget))

})


test_that_ref("2d_id", "2D identity bar", {
   # identity bar
   df <- data.frame(f=factor(c('aaaaaaa','bbbbbbbbb','ccccccccc','dddddddd')), val=c(5,6,2,8))
   eq_ref(plotluck(df, f, val))
})


test_that_ref("2d_spine", "2D spine", {

   # spine plot
   eq_ref(plotluck(as.data.frame(Titanic), Class, Survived, w=Freq))

   df <- as.data.frame(occupationalStatus)
   df$origin <- ordered(df$origin)
   df$destination <- ordered(df$destination)
   eq_ref(plotluck(df, origin, destination, w=Freq))
})

test_that_ref("3d_identity", "3D identity", {
   df <- data.frame(f=factor(c('aaaaaaa','bbbbbbbbb','ccccccccc','dddddddd')),
                    f2=factor(c(1,1,2,2)), val=c(5,6,2,8))
   eq_ref(plotluck(df, f, val, f2))

   eq_ref(plotluck(df, f, val, f2, opts=plotluck.options(max.colors.bar=0)))
})

test_that_ref("3d_heat", "3D heat", {
   eq_ref(plotluck(diamonds, cut, color, price))
})


test_that_ref("3d_spine", "3D spine", {
   eq_ref(plotluck(as.data.frame(Titanic), Class, Sex, Survived, w=Freq))
})


test_that_ref("3d_scatter", "3D scatter", {
   # using facets, hex
   eq_ref(plotluck(movies, length, rating, mpaa))

   # using facets
   eq_ref(plotluck(movies, length, rating, mpaa,
                   opts=plotluck.options(min.points.hex = 1E20)))

   # using colors
   eq_ref(plotluck(iris, Petal.Length, Petal.Width, Species))
})



test_that_ref("3d_density", "3D density", {
   eq_ref(plotluck(diamonds, price, cut, color)) # facets

   # facets, scatter
   eq_ref(plotluck(diamonds, price, cut, color,
                   opts=plotluck.options(min.points.density = 1E20)))

   # colors
   eq_ref(plotluck(diamonds, price, cut, color,
                   opts=plotluck.options(max.colors.density = 1E20)))
})



test_that_ref("3d_box", "3D box", {
   # color
   eq_ref(plotluck(movies, mpaa, rating, Action))

   # facets
   eq_ref(plotluck(movies, mpaa, rating, Action,
                   opts=plotluck.options(max.colors.box = 0)))

   # scatter, color
   eq_ref(plotluck(movies, mpaa, rating, Action,
                   opts=plotluck.options(min.points.box = 1E20)))

   # scatter, color, jitter
   eq_ref(plotluck(movies, mpaa, rating, Action,
                   opts=plotluck.options(min.points.box = 1E20, convert.duplicates.to.weights = FALSE)))

   # scatter, facets
   eq_ref(plotluck(movies, mpaa, rating, Action,
                   opts=plotluck.options(min.points.box = 1E20, max.colors.scatter = 0)))
})

test_that_ref("3d_discrete", "3D discretization", {

   eq_ref(plotluck(diamonds, cut, price, carat))
   eq_ref(plotluck(diamonds, cut, price, carat, opts=plotluck.options(min.points.box = 1E20)))
   eq_ref(plotluck(diamonds, cut, price, carat, opts=plotluck.options(min.points.box = 1E20, max.colors.scatter = 1E20)))
   eq_ref(plotluck(diamonds, cut, color, carat, opts=plotluck.options(min.size.grid.heat = 100)))
})


test_that_ref("3d_spine", "3D spine", {
   eq_ref(plotluck(as.data.frame(Titanic), Class, Sex, Survived, w=Freq))
})


test_that_ref("missing", "missing values", {
   set.seed(0)
   df<-data.frame(f=factor(sample(c(letters[1:3], NA), 20, replace=TRUE), exclude=FALSE),
                  f2=factor(sample(c(1,2,3,NA), 20, replace=TRUE), exclude=FALSE),
                  v=runif(20),
                  v2=runif(20))
   df$v2[c(1,5,6)] <- NA
   eq_ref(plotluck(df, f, v))
   eq_ref(plotluck(df, f, v, opts=plotluck.options(min.points.box = 1E20)))
   eq_ref(plotluck(df, v, f))
   eq_ref(plotluck(df, v, v2, f))
   eq_ref(plotluck(df, f, f2))
   eq_ref(plotluck(df, f, f2, opts=plotluck.options(exclude.factor=NA)))
})

test_that_ref("2d_weight", "instance weights", {
   set.seed(0)
   df<-data.frame(f=factor(sample(letters[1:3], 20, replace=TRUE), exclude=FALSE),
                  v1=runif(20),
                  v2=runif(20))
   df$w<-runif(20) + ifelse(df$v1>0.6, 3 * runif(20), ifelse(df$v1>0.3, 2*runif(20), 0))

   #violin
   eq_ref(plotluck(df, f, v1))
   eq_ref(plotluck(df, f, v1, w=w))

   # box
   eq_ref(plotluck(df, f, v1, opts=plotluck.options(use.geom.violin=FALSE)))
   eq_ref(plotluck(df, f, v1, w=w, opts=plotluck.options(use.geom.violin=FALSE)))

   # factor/num scatter
   eq_ref(plotluck(df, f, v1, opts=plotluck.options(min.points.box=1E20)))
   eq_ref(plotluck(df, f, v1, w=w, opts=plotluck.options(min.points.box=1E20)))

   # num/num scatter
   eq_ref(plotluck(df, v2, v1))
   eq_ref(plotluck(df, v2, v1, w=w))

   # density
   eq_ref(plotluck(df, v1, f, opts=plotluck.options(min.points.density=0)))
   eq_ref(plotluck(df, v1, f, w=w, opts=plotluck.options(min.points.density=0)))

   # histogram
   eq_ref(plotluck(df, v1, f,
                   opts=plotluck.options(min.points.density=0, use.geom.density=FALSE)))
   eq_ref(plotluck(df, v1, f, w=w,
                   opts=plotluck.options(min.points.density=0, use.geom.density=FALSE)))

   # heat map
   df<-expand.grid(1:5, 1:5)
   df <- rbind(df, df, df)
   df$v <- runif(75)
   df$w<-runif(75) + ifelse(df$v>0.6, 10 * runif(75), ifelse(df$v>0.3, 5*runif(75), 0))
   eq_ref(plotluck(df, Var1, Var2, v))
   eq_ref(plotluck(df, Var1, Var2, v, w=w))

})

test_that_ref("multi", "multiple plots", {
   eq_ref(plotluck.multi(diamonds))

   eq_ref(plotluck.multi(diamonds, price))

   eq_ref(plotluck.multi(diamonds, y=price))

})
