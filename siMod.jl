# Model permittivity of doped silicon, reproducing the model used in Fu, C. J., 
# & Zhang, Z. M. (2006). Nanoscale radiation heat transfer for silicon at different doping levels. 
# (International Journal of Heat and Mass Transfer, 49(9-10)), and 1703-1718 American Journal of 
# Physics 44, 353 (1976); doi: 10.1119/1.10194.
## Dependencies
using Roots
# Conversion factors are loaded from call through filmUtilities.jl
"""

	blkFD(arg::Float64)::Float64

Blakemore approximation of the Fermi-Dirac integral, accurate to with two percent of true values. 
"""
function blkFD(arg::Float64)::Float64

	if arg < -5.0 

		return exp(arg)
	elseif arg < 0.7

		return exp(arg) - 0.3 * exp(2.0 * arg) + 0.06 * exp(3.0 * arg)
	elseif arg < 1.1

		return exp(arg) * /(1.0, 1.0 + 0.27 * exp(arg))
	elseif arg < 1.3

		return /(4.0 * ^(arg, 1.5), 3.0 * sqrt(pi)) * (1.0 + /(1.15, arg^2))
	else

		return /(4.0 * ^(arg, 1.5), 3.0 * sqrt(pi)) * (1.0 + /(pi^2, 8.0 * arg^2))
	end
end
precompile(blkFD, (Float64, ))
"""
	
	siRspBck(enr::Float64, tmp::Float64)::ComplexF64

Temperature dependent background permittivity of intrinsic silicon. Loss is assumed to be 
negligible in comparison to free carrier contributions. 
"""
# Reproduction of the model reported in Lee, B. J., Zhang, Z. M., Early, E. A., DeWitt, D. P., & 
# Tsai, B. K. (2005). Modeling radiative properties of silicon with coatings and comparison with 
# reflectance measurements. Journal of thermophysics and heat transfer, 19(4), 558-565. 
@inline function siRspBck(enr::Float64, tmp::Float64)::ComplexF64

	return ^(sqrt(4.565 + /(97.3, ^(3.648,2) - ^(/(enr * 1.24, muEv), 2))) + (-1.864 * ^(10, -4) + 5.394 * /(^(10, -3), ^(3.648,2) - ^(/(enr * 1.24, muEv), 2))) * tmp, 2) - /(^(hBEv * 1.0e13, 2), enr * (enr + im * 1.0e10 * hBEv))
end
precompile(siRspBck, (Float64, Float64))
"""

	siEnrGap(tmp::Float64)::Float64 

Temperature dependent energy gap for intrinsic silicon in electron volts.
"""
@inline function siEnrGap(tmp::Float64)::Float64 

	return 1.17 - /(0.000473 * tmp^2.0, tmp + 636.0)
end
precompile(siEnrGap, (Float64, ))
"""

	itrnV(tmp::Float64)::Float64

Number of intrinsic valence band carriers in units of cm^{-3} as a function of temperature in 
Kelvin. 
"""
@inline function itrnV(tmp::Float64)::Float64

	return 2.66 * ^(10.0, 19.0) * ^(/(tmp, 300), 1.5)
end
precompile(itrnV, (Float64, ))
"""

	itrnC(tmp::Float64)::Float64

Number of intrinsic conduction band carriers in units of cm^{-3} as a function of temperature in 
Kelvin. 
"""
@inline function itrnC(tmp::Float64)::Float64 

	return 2.86 * ^(10.0, 19.0) * ^(/(tmp, 300), 1.5)
end
precompile(itrnC, (Float64, ))
"""

	function dnrCcn(enrFrm::Float64, tmp::Float64, acp::dptDsc)::Float64	

Acceptor concentration function. 
"""

@inline function acpCcn(enrFrm::Float64, tmp::Float64, acp::dptDsc)::Float64

	return acp.ccn * (1.0 - /(4.0, 4.0  + exp(/(enrFrm - acp.enr, blzK * tmp))))
end
precompile(acpCcn, (Float64, Float64, dptDsc))
"""

	function dnrCcn(enrFrm::Float64, tmp::Float64, dnr::dptDsc)::Float64

Donor concentration function. 
"""
@inline function dnrCcn(enrFrm::Float64, tmp::Float64, dnr::dptDsc)::Float64

	return dnr.ccn * (1.0 - /(2.0, 2.0  + exp(/(siEnrGap(tmp) - dnr.enr - enrFrm, blzK * tmp))))
end
precompile(dnrCcn, (Float64, Float64, dptDsc))
"""
	
	frmZero(enrFrm::Float64, tmp::Float64, dnr::dptVal, acp::dptVal)::Float64

Net charge concentration as a function of supposed Fermi energy and temperature. 
Used as an implicit equation for determining the correct Fermi energy. 
"""
@inline function frmZero(enrFrm::Float64, tmp::Float64, dnr::dptDsc, acp::dptDsc)::Float64

	return itrnV(tmp) * blkFD(-/(enrFrm, blzK * tmp)) + /(dnr.ccn, 1.0 + 2.0 * exp(/(enrFrm - siEnrGap(tmp) + dnr.enr, blzK * tmp))) - itrnC(tmp) * blkFD(/(enrFrm - siEnrGap(tmp), blzK * tmp)) - /(acp.ccn, 1.0 + 4.0 * exp(/(acp.enr - enrFrm, blzK * tmp)))
end
precompile(frmZero, (Float64, Float64, dptDsc, dptDsc))
"""

	prmMSi(tmp::Float64, dnr::dptDsc, acp::dptDsc)::siDsc

Computation of model parameters for doped silicon, see preamble of source file for details. 
"""
function prmMSi(tmp::Float64, dnr::dptDsc, acp::dptDsc)::siDsc

	# Self-consistent Fermi level.
	frmFnc(flvl) = frmZero(flvl, tmp, dnr, acp)
	enrFrm = find_zero(frmFnc, /(siEnrGap(tmp), 2.0))

	# Number of ionized acceptor and donor dopants (carriers), units cm^-3.
	cncDnr = dnrCcn(enrFrm, tmp, dnr)
	cncAcp = acpCcn(enrFrm, tmp, acp)

	# Room temperature total electron scattering time.
	stETR = (141.0 + /(19.5, 1.0 + ^(/(dnr.ccn, 1.3 * ^(10.0, 17.0)), 0.91))) * 1.0e-15
	# Room temperature total hole scattering time.
	stHTR = (10.0 + /(94.0, 1.0 + ^(/(acp.ccn, 1.9 * ^(10.0, 17.0)), 0.76))) * 1.0e-15

	# Room temperature electron-lattice scattering time.
	stELR = 2.23 * 1.0e-13
	# Room temperature hole-lattice scattering time.
	stHLR = 1.06 * 1.0e-13
	# Room temperature electron-impurity scattering time. 
	stEIR = /(stETR * stELR, stELR - stETR) 
	# Room temperature hole-impurity scattering time. 
	stHIR = /(stHTR * stHLR, stHLR - stHTR) 

	# Include temperature dependence scattering times.
	# Electrons 
	stEI = stEIR * ^(/(tmp, 300), 1.5)
	stEL = stELR * ^(/(tmp, 300), -3.8)
	# Holes
	stHI = stHIR * ^(/(tmp, 300), 1.5)
	stHL = stHLR * ^(/(tmp, 300), -3.6)

	# Total donor and acceptor decay rates.  
	dnrDcy = hBEv * (/(1.0, stEI) + /(1.0, stEL))
	acpDcy = hBEv * (/(1.0, stHI) + /(1.0, stHL))

	# Electron and hole effective masses.
	mEffE = 0.27 
	mEffH = 0.37

	# Conversion factors for donor and acceptor contributions.
	dnrFac = /(cncDnr * hBEv^2 * ^(10.0, 9.0), mEffE * 0.5199895 * 55.26349406 * 0.01112265)
	acpFac = /(cncAcp * hBEv^2 * ^(10.0, 9.0), mEffH * 0.5199895 * 55.26349406 * 0.01112265)

	return siDsc(tmp, dnrFac, acpFac, dnrDcy, acpDcy)
end
precompile(prmMSi, (Float64, dptDsc, dptDsc))