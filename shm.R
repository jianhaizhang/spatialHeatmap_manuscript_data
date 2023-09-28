
library(spatialHeatmap) 

# Multiple input aSVGs, two developmental stages.
set.seed(10)

df.random <- data.frame(matrix(runif(50, min=0, max=10), nrow=10))
colnames(df.random) <- c('shootTotalA__treatment1', 'shootTotalA__treatment2', 'shootTotalB__treatment1', 'shootTotalB__treatment2', 'notMapped') # Assign column names
rownames(df.random) <- paste0('gene', 1:10) # Assign row names 
df.random[1:2, ]


svg.sh1 <- "img/arabidopsis.thaliana_shoot_shm1.svg"
svg.sh2 <- "img/arabidopsis.thaliana_shoot_shm2.svg"
svg.sh.mul <- read_svg(c(svg.sh1, svg.sh2))


dat.mul.svg <- SPHM(svg=svg.sh.mul, bulk=df.random)
shm(data=dat.mul.svg, ID=c('gene2'), legend.r=0.25, legend.width=1, preserve.scale=TRUE, bar.width=0.09, line.color='grey50', sub.title.size=30, legend.text.size=25, legend.plot.title.size=30, legend.nrow=2)


