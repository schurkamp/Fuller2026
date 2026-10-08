# R Script: comparing bioacoustic and point count richness for Fuller 2026
## Finalized by Sam Schurkamp on July 1, 2026

# ==============================================================================================
# load preamble (required functions & data sets) 
# ==============================================================================================
source(here::here("scripts", "preamble.R"), echo = TRUE)

# ==============================================================================================
# statistical tests for richness comparisons
# ==============================================================================================

## Naming conventions:
### ao: ARU only
### ap: ARU plus PC
### pc: Point count only

pc.dat = all_data_long$PC
ao.dat = all_data_long$ARU
ap.dat = all_data_long$`ARU + PC`
n.dat = length(ao.dat)

# Dispersion tests. Is the data poisson-distributed?
## null is equidispersion, aternative is over- or underdispersion

## Build and run models
ap.model = glm(ap.dat ~ 1, family = poisson)
dispersiontest(ap.model, alternative = "less", trafo = 1)

ao.model = glm(ao.dat ~ 1, family = poisson)
dispersiontest(ao.model, alternative = "less", trafo = 1)

pc.model = glm(pc.dat ~ 1, family = poisson)
dispersiontest(pc.model, alternative = "less", trafo = 1)

## Alternative: goodness-of-fit test
### build models
ap.goodfit = goodfit(ap.dat, type = "poisson", method = "ML")
ao.goodfit = goodfit(ao.dat, type = "poisson", method = "ML")
pc.goodfit = goodfit(pc.dat, type = "poisson", method = "ML")

### run models
summary(ap.goodfit)
summary(ao.goodfit)
summary(pc.goodfit)

# Correlation, scatter plot of (1) ARU and ARU/PC vs PC and (2) ARU vs ARU/PC
## correlation
cor(pc.dat, ap.dat)
cor(pc.dat, ao.dat)
cor(ap.dat, ao.dat)

# scatter plot
par(mfrow = c(1,2))
par(mar = c(4,4,0.5,0.5))
plot(	jitter(pc.dat), jitter(ap.dat), 
      pch = 20, xlim = c(0,12), ylim = c(0,12),
      mgp = c(2.5,1,0),
      xlab = "Species richess (PC)",
      ylab = "Species richess (ARU and ARU+PC)")
points(jitter(pc.dat), jitter(ao.dat), pch = 20, col = "red")
abline(a = 0, b = 1, col = 8)
plot(	jitter(ao.dat), jitter(ap.dat), 
      pch = 20, xlim = c(0,12), ylim = c(0,12),
      mgp = c(2.5,1,0),
      xlab = "Species richess (ARU only)",
      ylab = "Species richess (ARU+PC)")
abline(a = 0, b = 1, col = 8)
points(jitter(pc.dat), jitter(ao.dat), pch = 20, col = "red")


# ==============================================================================================
# FIGURE: combine observed/theoretical graphs and their differences
# ==============================================================================================

# create objects that represent differences between techniques
ap.pc.dif	= ap.dat - pc.dat
ao.pc.dif	= ao.dat - pc.dat
ap.ao.dif 	= ap.dat - ao.dat

## Adjust for correlation
lambda_ap_adjust_pc = lambda_ap - cor(ap.dat, pc.dat) * sqrt(lambda_ap * lambda_pc)
lambda_pc_adjust_ap = lambda_pc - cor(ap.dat, pc.dat) * sqrt(lambda_ap * lambda_pc)
lambda_ao_adjust_pc = lambda_ao - cor(ao.dat, pc.dat) * sqrt(lambda_ao * lambda_pc)
lambda_pc_adjust_ao = lambda_pc - cor(ao.dat, pc.dat) * sqrt(lambda_ao * lambda_pc)
lambda_ao_adjust_ap = lambda_ao - cor(ao.dat, ap.dat) * sqrt(lambda_ao * lambda_ap)
lambda_ap_adjust_ao = lambda_ap - cor(ao.dat, ap.dat) * sqrt(lambda_ao * lambda_ap)

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

# COLUMN 1: theoretical vs. obserced distributions to evaluate dispersal and Poisson status ==============================================================================================

### Visualize the data; does it look Poisson or overdispersed?

# define theoretical lambda values with the same mean as the data
lambda_ap = mean(ap.dat)		# 6.70; ARU + PC
lambda_ao = mean(ao.dat)		# 6.30; ARU only		
lambda_pc	= mean(pc.dat)		# 3.13; PC only	

# object for parameters in the graph
dmax = 15
ymax = 0.356 # sets the max y-value
xmax = 16 # sets the max x-value
theo.w = 12 # sets width of the theoretical bars

layout(
  matrix(c(1, 4,
           2, 5,
           3, 6),
         nrow = 3, byrow = TRUE),
  heights = c(1, 1, 1)
)

par(oma = c(4, 0.5, 0.5, 0.5))

# margins for first column: y-axis on left
par(mar = c(1, 4.5, 0, 1))

# top plot: PC only, theoretical dist.
plot(	x = seq(0,30),
      y = dpois(x = seq(0,30), lambda = lambda_pc), 
      axes = FALSE,
      type = "h", lwd = theo.w, col = "lightsalmon", lend = 1,
      ylab = " ",
      xlim = c(0,xmax), ylim = c(0,ymax)); box()
axis(2); axis(1, labels = FALSE)

add_hanging_label(expression(italic("PC only")), cex = 1.4, corner = "topright")

# top plot: PC only, observed dist.
lines( 	x = sort(unique(pc.dat)),
        y = table(pc.dat)/n.dat,
        type = "h", lwd = 5.0, col = "black", lend = 1)

# put the legend in the PC plot (most white space)
legend(	x = 7.5, y = 0.2, bty = "n",
        c("Theoretical Dist.", "Observed Dist."),
        lty = c(1,1), lwd = c(5,3), cex = 1.3, 
        col = c("lightsalmon", "black"))

# mid plot: ARU only, theoretical dist.
plot(	x = seq(0,30),
      y = dpois(x = seq(0,30), lambda = lambda_ao), 
      axes = FALSE,
      type = "h", lwd = theo.w, col = "lightsalmon", lend = 1,
      ylab = " ", xlab = " ",
      xlim = c(0,xmax), ylim = c(0,ymax)); box()
axis(2); axis(1, labels = FALSE)
mtext("Probability Mass", 2, line = 3.2)
add_hanging_label(expression(italic("ARU only")), cex = 1.4, corner = "topright")

# mid plot: ARU only, observed dist.
lines( 	x = sort(unique(ao.dat)),
        y = table(ao.dat)/n.dat,
        type = "h", lwd = 5.0, col = "black", lend = 1)

# bottom plot: ARU and PC, theoretical dist.
plot(	x = seq(0,30),
      y = dpois(x = seq(0,30), lambda = lambda_ap), 
      axes = FALSE,
      type = "h", lwd = theo.w, col = "lightsalmon", lend = 1,
      ylab = "", 
      xlim = c(0,xmax), ylim = c(0,ymax)); box()
axis(2); axis(1)

# bottom plot: ARU and PC, observed dist.
lines( 	x = sort(unique(ap.dat)),
        y = table(ap.dat)/n.dat,
        type = "h", lwd = 5.0, col = "black", lend = 1)

add_hanging_label(expression(italic("ARU + PC")), cex = 1.4, corner = "topright")

mtext("Species Richess", 1, line = 2.95)


# COLUMN 2: actual distributions vs. a theoretical distribution with no differences ==============================================================================================

# margins for second column: y-axis on right
par(mar = c(1, 1.5, 0, 4.5))

## Assess distribution of pairwise differences 
## Theoretical distribution of the difference between two paired Poisson samples is Skellam 
## (also sometimes called a Bessel function -- same thing)
## Play around with whether the correlation affects this....

plot.a = table(ap.pc.dif); plot.a
plot.b = table(ao.pc.dif); plot.b
plot.c = table(ap.ao.dif); plot.c

ymax = 0.30

#### Plot theoretical skellam distribution WITH correlation vs observed differences
plot(	x = seq(-dmax,dmax),
      y = dskellam(seq(-dmax,dmax), lambda1 = lambda_pc, lambda2 = lambda_pc),
      axes = FALSE,
      type = "h", lwd = theo.w, col = "lightsalmon", lend = 1,
      ylab = "", 
      xlim = c(-10,10), ylim = c(0,ymax)); box()
axis(4); axis(1, labels = FALSE)

lines( 	x = sort(unique(ao.pc.dif)),
        y = table(ao.pc.dif)/n.dat,
        type = "h", lwd = 5.0, col = "black", lend = 1)
abline(v = 0, lty = 3, lwd = 2, col = "gray40")

add_hanging_label(expression(italic("ARU v. PC")), cex = 1.4, corner = "topleft")

plot(	x = seq(-dmax,dmax),
      y = dskellam(seq(-dmax,dmax), lambda1 = lambda_pc, lambda2 = lambda_pc),
      axes = FALSE,
      type = "h", lwd = theo.w, col = "lightsalmon", lend = 1,
      ylab = "", 
      xlim = c(-10,10), ylim = c(0,ymax)); box()
axis(4); axis(1, labels = FALSE)

lines( 	x = sort(unique(ap.pc.dif)),
        y = table(ap.pc.dif)/n.dat,
        type = "h", lwd = 5.0, col = "black", lend = 1)
abline(v = 0, lty = 3, lwd = 2, col = "gray40")

add_hanging_label(expression(italic("ARU + PC vs. PC")), cex = 1.4, corner = "topleft")

plot(	x = seq(-dmax,dmax),
      y = dskellam(seq(-dmax,dmax), lambda1 = lambda_ao, lambda2 = lambda_ao),
      axes = FALSE,
      type = "h", lwd = theo.w, col = "lightsalmon", lend = 1,
      ylab = "", 
      xlim = c(-10,10), ylim = c(0,0.75)); box()		# use different y limits in this panel only
axis(4); axis(1)
mtext("Difference in Species Richess", 1, line = 2.95)

lines( 	x = sort(unique(ap.ao.dif)),
        y = table(ap.ao.dif)/n.dat,
        type = "h", lwd = 5.0, col = "black", lend = 1)
abline(v = 0, lty = 3, lwd = 2, col = "gray40")

add_hanging_label(expression(italic("ARU + PC v. ARU")), cex = 1.4, corner = "topleft")



