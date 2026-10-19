# GOAL: analyze differences between observation protocols
# NOTE: the figure used in the publication is in the "richness_comparisons.Rmd" file.

# ==============================================================================================
# load preamble (required functions & data sets) 
# ==============================================================================================
source(here::here("scripts", "preamble.R"), echo = TRUE)

# ==============================================================================================
# define relevant parameters
# ==============================================================================================
## Naming conventions:
### ao: ARU only
### ap: ARU plus PC
### pc: Point count (only)
pc.dat = all_data_long$PC
ao.dat = all_data_long$ARU
n.dat = length(ao.dat)

# define theoretical lambda values with the same mean as the data
lambda_ap = mean(ap.dat)		# 6.70; ARU + PC
lambda_ao = mean(ao.dat)		# 6.30; ARU only		
lambda_pc	= mean(pc.dat)		# 3.13; PC only	

# differences between testing parameters
ao.pc.dif	= ao.dat - pc.dat

# object for parameters in the graph
dmax = 15
ymax = 0.356 # sets the max y-value
xmax = 16 # sets the max x-value
theo.w = 12 # sets width of the theoretical bars

# custom margins for the tick marks of the smaller graphs
small_y_ticks <- seq(0, ymax, length.out = 4)
small_y_labels <- sprintf("%.1f", small_y_ticks)

# ==============================================================================================
# Build function for a hanging label
# ==============================================================================================

# add a function so you can label each plot within the figure later on
add_hanging_label <- function(label, 
                              corner = c("topleft", "topright"),
                              cex = 1.3, 
                              bg = "white", 
                              border = "black") {
  
  corner <- match.arg(corner)
  
  u <- par("usr")
  
  pad.x <- diff(u[1:2]) * 0.03
  pad.y <- diff(u[3:4]) * 0.035
  
  label.width  <- strwidth(label, cex = cex)
  label.height <- strheight(label, cex = cex)
  
  box.width  <- label.width + 2 * pad.x
  box.height <- label.height + 2 * pad.y
  
  y.top    <- u[4]
  y.bottom <- y.top - box.height
  
  if (corner == "topleft") {
    x.left  <- u[1]
    x.right <- u[1] + box.width
    x.text  <- u[1] + box.width / 2
  }
  
  if (corner == "topright") {
    x.right <- u[2]
    x.left  <- u[2] - box.width
    x.text  <- u[2] - box.width / 2
  }
  
  rect(
    xleft   = x.left,
    ybottom = y.bottom,
    xright  = x.right,
    ytop    = y.top,
    col = bg,
    border = border,
    xpd = NA
  )
  
  text(
    x = x.text,
    y = y.top - box.height / 2,
    labels = label,
    cex = cex,
    xpd = NA
  )
}
# ==============================================================================================
# BUILD GRAPH
# ==============================================================================================

# Generate the PNG
png(file = "distributionsUPDATE.png", width = 7.5, height = 5, units = "in", res = 300)


# Establish layout. This code makes the left column twice as big as the right.
## Can tweak as needed.
layout(
  matrix(c(2, 1,
           3, 1),
         nrow = 2, byrow = TRUE),
  widths = c(1.3, 3),
  heights = c(1, 1)
)

# define outer margins
par(oma = c(4, 1, 2, 1))

#### ARU vs. PC: the main plot ==============================================================================================
# margins for main plot
par(mar = c(0, 0, 0, 4))

# add theoretical distribution
plot(	x = seq(-dmax,dmax),
      y = dskellam(seq(-dmax,dmax), lambda1 = lambda_pc, lambda2 = lambda_pc),
      axes = FALSE,
      type = "h", lwd = theo.w, col = "lightsalmon", lend = 1,
      ylab = "",
      xlab = "",
      xlim = c(-10,10), ylim = c(0,ymax)); box()

# add marks for the axes
axis(4); axis(1, labels = TRUE)

# add observed distribution
lines( 	x = sort(unique(ao.pc.dif)),
        y = table(ao.pc.dif)/n.dat,
        type = "h", lwd = 6.0, col = "black", lend = 1)

# add line at zero for visual aid
abline(v = 0, lty = 3, lwd = 2, col = "gray40")

# add label
add_hanging_label(expression(italic("ARU v. PC")), cex = 1.4, corner = "topleft")

# put the legend in the white space
legend(x = -10, y = 0.25, # position using graph's coordinates
       bty = "n", # specify that I don't want a border around the legend
        c("Theoretical Dist.", "Observed Dist."), # the two things the legend should say
        lty = c(1,1), # tells R to make both things a solid line
       lwd = c(13,13), # sets line width of each thing (width is along the y-axis, weirdly...)
       seg.len = c(1,1), # x-axis length
       cex = 1.3, # font size
        col = c("lightsalmon", "black")) # line colors

# big y-axis label
mtext("Probability Mass", 4, line = 3)
# x-axis label
mtext("Difference in Species Richess", 1, line = 2.5)



# top plot: PC distribution ==============================================================================================

# margins 
par(mar = c(0, 2, 0, 0.5))

# plot 
plot(	x = seq(0,30),
      y = dpois(x = seq(0,30), lambda = lambda_pc), 
      axes = FALSE,
      type = "h", lwd = 8, col = "lightsalmon", lend = 1,
      ylab = "",
      xlab = "",
      xlim = c(0,xmax), ylim = c(0,ymax)); box()
axis(side = 2,
       at = small_y_ticks,
       labels = small_y_labels,
       las = 1); axis(1, labels = FALSE)

# label the plot
add_hanging_label(expression(italic("PC")), cex = 1.4, corner = "topright")

# top plot: PC only, observed dist. 
lines( 	x = sort(unique(pc.dat)),
        y = table(pc.dat)/n.dat,
        type = "h", lwd = 4, col = "black", lend = 1)

# bottom left plot: ARU dist. ==============================================================================================

# margins
par(mar = c(0, 2, 0, 0.5))

# build plot
plot(	x = seq(0,30),
      y = dpois(x = seq(0,30), lambda = lambda_ao), 
      axes = FALSE,
      type = "h", lwd = 8, col = "lightsalmon", lend = 1,
      ylab = "", xlab = "",
      xlim = c(0,xmax), ylim = c(0,ymax)); box()
axis(  side = 2,
       at = small_y_ticks,
       labels = small_y_labels,
       las = 1); axis(1, labels = TRUE)

# add observed dist.
lines( 	x = sort(unique(ao.dat)),
        y = table(ao.dat)/n.dat,
        type = "h", lwd = 4, col = "black", lend = 1)

# label the plot
add_hanging_label(expression(italic("ARU")), cex = 1.4, corner = "topright")
# label the x-axes
mtext("Species Richess", 1, line = 2.5)


# That's the plot! Save it:
dev.off()
 # 

