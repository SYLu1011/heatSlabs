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