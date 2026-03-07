# Use this with SOL Super Computer when loading graphs/plots
###############################################################
# if there is a Error in dev.control(displaylist = "enable") : 
# dev.control() called without an open graphics device

# Restart R in Session > Restart R

dev.list()
# If this returns NULL

dev.new()
# This may also return NULL
# This might also return "In (function ()  : Only one RStudio graphics device is permitted"

options(device = "RStudioGD")

dev.list()
# Should return a list with RStudioGD in it

# Try running the plot again, should work now

