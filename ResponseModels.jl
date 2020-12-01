module ResponseModels
export hStep, cstRsp, sicRsp, siDsc, dptDsc, siRsp, prmMSi
# Collection of material permittivity models for use in heatLayers code.
# Conversion factors.
# Microns to electron volts.
const muEv = 1.239
# Radial frequency to electron volts.
const hBEv = 6.582e-16
# Temperature to electron volts.
const blzK = 8.6173e-5
# Electron volts to Joules.
const evJ = 1.602176565e-19
## Constant value for testing.
function cstRsp(cst::ComplexF64, enr::Float64)::ComplexF64

	return cst
end
## Heaviside function for testing
function hStep(arg::Float64)::Float64

	if arg >= 0.0

		return 1.0
	else
		return 0.0
	end
end
### Silicon carbide.
"""

	sicRsp(enr::Float64)::ComplexF64

Simplified permittivity response model for silicon carbide. 
Approximately valid for energies between 0.1 and 10 electron volts. 
"""
# Reproduction of the model reported in J. Appl. Phys. 106, 044306 (2009).
function sicRsp(enr::Float64)::ComplexF64

	# Energy values (in electron volts) for the longitudinal and transverse optical phonons, 
	# and decay rate. 
	enrL = 0.12014342321784009
	enrT = 0.09831711568536007
	enrD = 0.0005901528146448003

	return 6.7 * (1.0 + /(enrL^2.0 - enrT^2.0, enrT^2.0 - enr^2.0 - im * enrD * enr))
end
### Silicon 
# Sample usage
# Boron p-type doping
# brnDpt = dptDsc(0.044, 1.0e18)
# Empty doping
# nllDpt = dptDsc(0.0, 0.0)
# Calculate silicon model parameters
# siModE = prmMSi(tmpLst[1], nllDpt, nllDpt)
# siRspE(enr) = siRsp(enr, siModE)
"""
Data structure holding computed model parameters for the permittivity response of silicon for 
given dopant properties at a fixed temperature. Call prmMSi(tmp::Float64, dnr::dptDsc, 
acp::dptDsc)::siDsc for construction. 

# Arguments
.tmp: temperature at which parameters where calculated. 
.dnrFac: strength of donor impurity Lorentzian.
.acpFac: strength of acceptor impurity Lorentzian.
.dnrDcy: decay rate (gamma) for donor impurity Lorentzian.
.acpDcy: decay rate (gamma) for acceptor impurity Lorentzian.
"""
struct siDsc

	tmp::Float64
	dnrFac::Float64
	acpFac::Float64
	dnrDcy::Float64
	acpDcy::Float64
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
# Include model code
include("siMod.jl")
"""

	siRsp(enr::Float64, tmp::Float64, siVal::siDsc)::ComplexF64

Permittivity response model for doped silicon. See source file for additional details. 
"""
function siRsp(enr::Float64, siVal::siDsc)::ComplexF64

	return siRspBck(enr, siVal.tmp) - /(siVal.dnrFac, enr^2 + im * enr * siVal.dnrDcy) - /(siVal.acpFac, enr^2 + im * enr * siVal.acpDcy)
end
end