
# global helper objects -- available to all modules

def restrict_bounds(x, lower, upper):
	''' Return x if in bounds, otherwise return the bound it crossed. '''
	return max(lower, min(x, upper))
	