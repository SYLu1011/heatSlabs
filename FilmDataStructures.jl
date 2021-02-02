__precompile__()
# See readme for module description. 
module FilmDataStructures
export lyrDsc, dptDsc
# Data structure definition for calculating heat transfer between two layered infinite half spaces.  
"""
Data structure for permittivity profile of layered structure.
Each provided permittivity in the permittivity list is presumed to be constant within an interval beginning on the previous interface depth, and terminating on the next interface depth (infinity when no other boundary exists). 

# Arguments
.bdrLoc: relative positions of interfaces, units of microns. 
.tmpLst: assumed (constant) temperature of each layer in Kelvin.
.rspPrf: permittivity values for each depth.
.tfrFac: layer transfer factor, imaginary part of electric susceptibility for heat transfer. 
"""
struct lyrDsc

	bdrLoc::Array{Float64,1}
	tmpLst::Array{Float64,1}
	rspPrf::Any
	tfrFac::Any
end
"""
Data structure describing properties of supposed dopants.

# Arguments
.enr: offset of impurity energy from band edge in electron volts.
.ccn: impurity concentration in electron volts. 
"""
# Antimonide (n-type)
# enr = 0.039
# Phosphorous, most common n-type
# enr = 0.044
# Boron, most common p-type
# enr = 0.045
struct dptDsc

	enr::Float64
	ccn::Float64
end
end