# GOAL: analyze differences between observation protocols
# NOTE: the figure used for the publication is at the last section of this code.

# ==============================================================================================
# Load required libraries
# ==============================================================================================

## organizational functions
library(tidyr)
library(dplyr)
library(stringr)
library(lubridate)
library(kableExtra)
library(data.table)
library(plyr)                                                 
library(ggh4x)
library(reshape)

# used during the second lollipop graph prep
library(forcats)

# improve the import functions for data within repositories
library(readr)
library(here)

# graphing functions
library(ggplot2)
library(patchwork) # required for "plot layout" when merging lollipop graphs

# packages used in summary data sets
library(flextable)
library(officer)
library(purrr)

# package for running statistical tests
library(car)

# for the Bayesian summaries
library(AER)
library(vcd)
library(skellam)



################################################    
# Step 1: import and prepare a combined dataset
################################################ 

# read the appropriate CSV from GitHub into R
all_pc_aru_IMPORT <- readr::read_csv(here::here("CSVs/analysis", "all_occupancy_analysis.csv"))
# 'problem' upon reading is due to columns we don't use. Okay to ignore

# calculate richness for ARU and PC separately
pc_aru_richness <- all_pc_aru_IMPORT %>%
  filter(!CommonName %in% c("----", "Trumpeter Swan", "Sandhill Crane", "Snowy Egret", "Forster's Tern", "Sedge Wren", "Common Tern", "Wilson's Snipe")) %>% # had to add the sedge wren and common tern
  dplyr::filter(IndivCount >=1) %>%
  dplyr::group_by(Protocol, Route, Year, Point) %>%
  dplyr::summarize(Richness = length(unique(BirdCd)))

# calculate richness when ARU + PC are added together
pc_aru_richness_added_together <- all_pc_aru_IMPORT %>%
  filter(!CommonName %in% c("----", "Trumpeter Swan", "Sandhill Crane", "Snowy Egret", "Forster's Tern", "Sedge Wren", "Common Tern", "Wilson's Snipe")) %>%
  dplyr::filter(IndivCount >=1) %>%
  dplyr::group_by(Route, Point, Year) %>% # exclude Protocol to get combined richness
  dplyr::summarize(Richness = length(unique(BirdCd))) %>%
  dplyr::mutate(Protocol = "ARU + PC")

# combine the resulting DFs into one DF that has 3 protocol (ARU alone, PC alone, ARU + PC combined)
all_richness_no0 <- rbind(pc_aru_richness, pc_aru_richness_added_together)

# zero-fill richness summary  
all_data_plot <- all_richness_no0 %>%
  dplyr::select(Point, Protocol, Richness, Route, Year) %>%
  pivot_wider(names_from = Protocol, values_from = Richness) %>%
  pivot_longer(cols = c(ARU, PC, `ARU + PC`), names_to = "Protocol", values_to = "Richness") %>%
  dplyr::mutate(Richness = ifelse(is.na(Richness), 0, Richness)) 

# convert to wide-form data
all_data_long <- cast(all_data_plot, Point+Route+Year~Protocol, value = "Richness")





## Naming conventions:
### ao: ARU only
### ap: ARU plus PC
### pc: Point count (only)

pc.dat = all_data_long$PC
ao.dat = all_data_long$ARU
ap.dat = all_data_long$`ARU + PC`
n.dat = length(ao.dat)

########################################################################
########################################################################
## 2. Look at the data; 
#### a. see if it appears Poisson and overlay theoretical Poisson distribution
#### b. see if there is correlation (will effect how we model differences)
########################################################################
########################################################################

### Visualize the data; does it look Poisson or overdispersed?

# define theoretical lambda values with the same mean as the data
lambda_ap = mean(ap.dat)		# 6.70; ARU + PC
lambda_ao = mean(ao.dat)		# 6.30; ARU only		
lambda_pc	= mean(pc.dat)		# 3.13; PC only	

ymax = 0.35 # sets the max y-value
xmax = 16 # sets the max x-value
theo.w = 15 # sets width of the theoretical bars

par(mfrow = c(3,1))
par(mar = c(1, 4.5, 0, 0.5))
par(oma = c(4, 0, 2, 0))

# top plot: ARU and PC, theoretical dist.
plot(	x = seq(0,30),
      y = dpois(x = seq(0,30), lambda = lambda_ap), 
      axes = FALSE,
      type = "h", lwd = theo.w, col = "lightsalmon", lend = 1,
      ylab = "", 
      xlim = c(0,xmax), ylim = c(0,ymax)); box()
axis(2); axis(1, labels = FALSE)

# top plot: ARU and PC, observed dist.
lines( 	x = sort(unique(ap.dat)),
        y = table(ap.dat)/n.dat,
        type = "h", lwd = 5.0, col = "black", lend = 1)

# legend, which goes in top plot
legend(	x = 9, y = 0.4, bty = "n",
        c("Theoretical distribution", "Observed distribution"),
        lty = c(1,1), lwd = c(5,3), cex = 1.3, 
        col = c("lightsalmon", "black"))

# mid plot: ARU only, theoretical dist.
plot(	x = seq(0,30),
      y = dpois(x = seq(0,30), lambda = lambda_ao), 
      axes = FALSE,
      type = "h", lwd = theo.w, col = "lightsalmon", lend = 1,
      ylab = " ", 
      xlim = c(0,xmax), ylim = c(0,ymax)); box()
axis(2); axis(1, labels = FALSE)
mtext("Probability Mass", 2, line = 3.2)

# mid plot: ARU only, observed dist.
lines( 	x = sort(unique(ao.dat)),
        y = table(ao.dat)/n.dat,
        type = "h", lwd = 5.0, col = "black", lend = 1)

# bottom plot: PC only, theoretical dist.
plot(	x = seq(0,30),
      y = dpois(x = seq(0,30), lambda = lambda_pc), 
      axes = FALSE,
      type = "h", lwd = theo.w, col = "lightsalmon", lend = 1,
      ylab = " ", xlab = " ",
      xlim = c(0,xmax), ylim = c(0,ymax)); box()
axis(2); axis(1)
mtext("Species Richess", 1, line = 2.95)

# bottom plot: PC only, observed dist.
lines( 	x = sort(unique(pc.dat)),
        y = table(pc.dat)/n.dat,
        type = "h", lwd = 5.0, col = "black", lend = 1)


########################################################################
## b. Correlation and scatter plot of (1) ARU and ARU/PC vs PC and (2) ARU vs ARU/PC
########################################################################

cor(pc.dat, ap.dat)
cor(pc.dat, ao.dat)
cor(ap.dat, ao.dat)

par(mfrow = c(1,2))
par(mar = c(4,4,0.5,0.5))
#par(oma = c(4,4,0.2,0.2))

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


########################################################################
########################################################################
## 3. Formally test whether data is underdispersed
########################################################################
########################################################################

### For dispersion test (AER package), null is equidispersion, laternative is over or underdispersion
ap.model = glm(ap.dat ~ 1, family = poisson)
dispersiontest(ap.model, alternative = "less", trafo = 1)

ao.model = glm(ao.dat ~ 1, family = poisson)
dispersiontest(ao.model, alternative = "less", trafo = 1)

pc.model = glm(pc.dat ~ 1, family = poisson)
dispersiontest(pc.model, alternative = "less", trafo = 1)

### Also could use goodness-of-fit test

ap.goodfit = goodfit(ap.dat, type = "poisson", method = "ML")
ao.goodfit = goodfit(ao.dat, type = "poisson", method = "ML")
pc.goodfit = goodfit(pc.dat, type = "poisson", method = "ML")

summary(ap.goodfit)
summary(ao.goodfit)
summary(pc.goodfit)

########################################################################
########################################################################
## 4. Test for differences between testing techniques
########################################################################
########################################################################

ap.pc.dif	= ap.dat - pc.dat
ao.pc.dif	= ao.dat - pc.dat
ap.ao.dif 	= ap.dat - ao.dat

########################################################################
### 4a. Visual appraisal
########################################################################

## Assess distribution of pairwise differences 
## Theoretical distribution of the difference between two paired Poisson samples is Skellam 
## (also sometimes called a Bessel function -- same thing)
## Play around with whether the correlation affects this....

plot.a = table(ap.pc.dif); plot.a
plot.b = table(ao.pc.dif); plot.b
plot.c = table(ap.ao.dif); plot.c

dmax = 15
ymax = 0.30

## Adjust for correlation;

lambda_ap_adjust_pc = lambda_ap - cor(ap.dat, pc.dat) * sqrt(lambda_ap * lambda_pc)
lambda_pc_adjust_ap = lambda_pc - cor(ap.dat, pc.dat) * sqrt(lambda_ap * lambda_pc)

lambda_ao_adjust_pc = lambda_ao - cor(ao.dat, pc.dat) * sqrt(lambda_ao * lambda_pc)
lambda_pc_adjust_ao = lambda_pc - cor(ao.dat, pc.dat) * sqrt(lambda_ao * lambda_pc)

lambda_ao_adjust_ap = lambda_ao - cor(ao.dat, ap.dat) * sqrt(lambda_ao * lambda_ap)
lambda_ap_adjust_ao = lambda_ap - cor(ao.dat, ap.dat) * sqrt(lambda_ao * lambda_ap)

#### Plot theoretical skellam distribution WITH correlation vs observed differences

par(mfrow = c(3,1))
par(mar = c(1, 4.5, 0, 0.5))
par(oma = c(4, 0, 2, 0))

plot(	x = seq(-dmax,dmax),
      # y = dskellam(seq(-dmax,dmax), lambda1 = lambda_ap, lambda2 = lambda_pc),
      y = dskellam(seq(-dmax,dmax), lambda1 = lambda_ap_adjust_pc, lambda2 = lambda_pc_adjust_ap),
      axes = FALSE,
      type = "h", lwd = theo.w, col = "lightgray", lend = 1,
      ylab = "", 
      xlim = c(-dmax,dmax), ylim = c(0,ymax)); box()
axis(2); axis(1, labels = FALSE)

lines( 	x = sort(unique(ap.pc.dif)),
        y = table(ap.pc.dif)/n.dat,
        type = "h", lwd = 2.0, col = "black", lend = 1)

legend(	x = -15, y = 0.3, bty = "n",
        c("Theoretical distribution", "Observed distribution"),
        lty = c(1,1), lwd = c(5,3), cex = 1.3, 
        col = c("Light grey", "black"))

plot(	x = seq(-dmax,dmax),
      # y = dskellam(seq(-dmax,dmax), lambda1 = lambda_ao, lambda2 = lambda_pc),
      y = dskellam(seq(-dmax,dmax), lambda1 = lambda_ao_adjust_pc, lambda2 = lambda_pc_adjust_ao),
      axes = FALSE,
      type = "h", lwd = theo.w, col = "lightgrey", lend = 1,
      ylab = "", 
      xlim = c(-dmax,dmax), ylim = c(0,ymax)); box()
axis(2); axis(1, labels = FALSE)
mtext("Probability mass (or relative frequency)", 2, line = 3.2)

lines( 	x = sort(unique(ao.pc.dif)),
        y = table(ao.pc.dif)/n.dat,
        type = "h", lwd = 2.0, col = "black", lend = 1)

plot(	x = seq(-dmax,dmax),
      # y = dskellam(seq(-dmax,dmax), lambda1 = lambda_ap, lambda2 = lambda_ao),
      y = dskellam(seq(-dmax,dmax), lambda1 = lambda_ap_adjust_ao, lambda2 = lambda_ao_adjust_ap),
      axes = FALSE,
      type = "h", lwd = theo.w, col = "light grey", lend = 1,
      ylab = "", 
      xlim = c(-dmax,dmax), ylim = c(0,0.75)); box()		# use different y limits in this panel only
axis(2); axis(1)
mtext("Difference in species richess", 1, line = 2.95)

lines( 	x = sort(unique(ap.ao.dif)),
        y = table(ap.ao.dif)/n.dat,
        type = "h", lwd = 2.0, col = "black", lend = 1)


########################################################################
### 4b. Formal test for statistical significance of difference
########################################################################

## Best approach, after a bit of reading and research...
## .... appears to be the Conditional method (or Binomial test)
## NOTE: Some of this code was generated with guidance from Google AI

# For this test, only paired counts that are DIFFERENT have "value"
# (need to read more about this, but let's do it for now)

sum(ap.pc.dif > 0); sum(!ap.pc.dif == 0)
ap.pc.test = binom.test(x = sum(ap.pc.dif > 0), n = sum(!ap.pc.dif == 0), alternative = "greater") 
print(ap.pc.test)

sum(ao.pc.dif > 0); sum(!ao.pc.dif == 0)
ao.pc.test = binom.test(x = sum(ao.pc.dif > 0), n = sum(!ao.pc.dif == 0), alternative = "greater") 
print(ao.pc.test)

sum(ap.ao.dif > 0); sum(!ap.ao.dif == 0)
ap.ao.test = binom.test(x = sum(ap.ao.dif > 0), n = sum(!ap.ao.dif == 0), alternative = "greater") 
print(ap.ao.test)



#############################################################################
### END OF OFFICIAL CODE FOR MANUSCRIPT
#############################################################################





#################################################################
#################################################################
#################################################################
## Extra: there's also a poisson ratio test (if the two posison variables are no different, 
## the distribution of their ratios will be close to 1).
#################################################################
#################################################################
#################################################################

#ratios = rpois(10000,lambda_aru)/rpois(10000,lambda_pc)
#clean.r = ratios[!is.infinite(ratios)]

#hist(clean.r, freq = F)




########################################################################
## ...looking at ratios instead of (or in addition to) differences (let's see what reviewers ask for)

#par(mfrow = c(2,3))
#hist(ao.pc.dif, xlim = c(-10,10)); hist(ap.pc.dif, xlim = c(-10,10)); hist(ap.ao.dif, xlim = c(-10,10))
#rat.xlim = c(0.1,10)
#hist(ao.pc.rat, log = "x", xlim = rat.xlim); hist(ap.pc.rat, log = "x", xlim = rat.xlim); hist(ap.ao.rat, log = "x", xlim = rat.xlim)














# ==============================================================================================
# PUB FIGURE: combine observed/theoretical graphs and their differences
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
# COLUMN 1
# ==============================================================================================
### Visualize the data; does it look Poisson or overdispersed?

# define theoretical lambda values with the same mean as the data
lambda_ap = mean(ap.dat)		# 6.70; ARU + PC
lambda_ao = mean(ao.dat)		# 6.30; ARU only		
lambda_pc	= mean(pc.dat)		# 3.13; PC only	

ymax = 0.356 # sets the max y-value
xmax = 16 # sets the max x-value
theo.w = 12 # sets width of the theoretical bars

# Generate the PNG
png(file = "distributions.png", width = 8, height = 10, units = "in", res = 300)
#layout(
 # matrix(c(1, 1, 
  #         2, 5,
   #        3, 6,
    #       4, 7),
     #    nrow = 4, byrow = TRUE),
#  heights = c(0.35, 1, 1, 1)
#)


layout(
  matrix(c(1, 4,
           2, 5,
           3, 6),
         nrow = 3, byrow = TRUE),
  heights = c(1, 1, 1)
)

par(oma = c(4, 0.5, 0.5, 0.5))





#par(oma = c(4, 0.5, 2, 0))

# shared legend above all six graphs
#par(mar = c(0, 0, 0, 0))
#plot.new()

# put legend above all plots
#legend(
#  "center",
#  legend = c("Theoretical Distribution", "Observed Distribution"),
 # lty = c(1, 1),
#  lwd = c(6, 6),
#  cex = 1.6,
#  col = c("lightsalmon", "black"),
#  bty = "n",
#  horiz = TRUE
#)

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

# ==============================================================================================
# COLUMN 2
# ==============================================================================================

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

# Generate the values for the theoretical distributions
# NOTE! We decided to display the theoretical distributions assuming no 
# correlation and no significant differences. 
# All the following commented out code is our original code,
# Which was a theoretical distribution of the correlation and relationship.
# That code's N_o was simply that the dist. is Poisson... less interesting!
## lambda adjustments for plot 1
#lambda_ap_adjust_pc = lambda_ap - cor(ap.dat, pc.dat) * sqrt(lambda_ap * lambda_pc)
#lambda_pc_adjust_ap = lambda_pc - cor(ap.dat, pc.dat) * sqrt(lambda_ap * lambda_pc)

## adjustments for plot 2
#lambda_ao_adjust_pc = lambda_ao - cor(ao.dat, pc.dat) * sqrt(lambda_ao * lambda_pc)
#lambda_pc_adjust_ao = lambda_pc - cor(ao.dat, pc.dat) * sqrt(lambda_ao * lambda_pc)

## adjustments for plot 3
#lambda_ao_adjust_ap = lambda_ao - cor(ao.dat, ap.dat) * sqrt(lambda_ao * lambda_ap)
#lambda_ap_adjust_ao = lambda_ap - cor(ao.dat, ap.dat) * sqrt(lambda_ao * lambda_ap)

# INSTEAD: skellam lambdas w/o correlation, sig diff
## note 2: we also determined that we should always use PC, except in ARUvsARUPC, where we use ARU.

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

# That's the plot! Save it:
dev.off()












