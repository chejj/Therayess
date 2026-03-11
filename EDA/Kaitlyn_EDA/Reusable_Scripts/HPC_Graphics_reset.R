# Use this with SOL Super Computer when loading graphs/plots
###############################################################


# if there is a Error in dev.control(displaylist = "enable") : 
# dev.control() called without an open graphics device

# --- Restart R in Session > Restart R ---
#  ↓↓↓↓

dev.list() # If this returns NULL
dev.off()
dev.new() # This may also return NULL # This might also return "In (function ()  : Only one RStudio graphics device is permitted"
options(device = "RStudioGD") # may need to run this a bunch of times
dev.list() # Should return a list with RStudioGD in it
plot(1:10) # test plot

# Try running the plot again, should work now