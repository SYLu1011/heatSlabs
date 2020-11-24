# module FilmUtilities
using Cubature
export lyrDsc, femFlx, plnFnc, flxPfc, flxFncPR
## Constants
# Order parameter for cubature.
const qudOrd = 32
# Relative tolerance for cubature.
const relTol = 1.0e-3
# Cutoff for unbounded integrals tan and atan variable transforms.
# tan(pi/2 - intDltInf) scaling. 
# Use number smaller than 0.5 * ^(10.0, -4) for sub nanometer gaps. 
const intDltInf = 1.0e-5
# Initial guess for width of ``intermediate wave'' behavior. 
const lcWdt = 1.0e-3
# Desired relative magnitude for defining width of ``intermediate wave'' behavior. 
const lcRel = 1.0e-2
# Ratio of /(real(rsp(ener)), imag(rsp(ener))) before result is treated as intermediate wave. 
const iTrnV = 1.0e3
# Near light line normalization, effectively number of wavelengths before decay.
const mxSz = 50.0
# Near light regularization, imposed as imaginary part of perpendicular wave vector. 
const lcReg = 1.0e-2
# Conversion factors.
# Microns to electron volts.
const muEv = 1.239
# Radial frequency to electron volts.
const hBEv = 6.582e-16
# Temperature to electron volts.
const blzK = 8.6173e-5
# Electron volts to Joules.
const evJ = 1.602176565e-19

## Load data structures. 
include("filmDataStructures.jl")
# Repository functions no longer being used, but possibly helpful for testing. 
# include("filmUtilitiesRepo.jl")
# Functions definitions for calculating heat transfer between two layered infinite half spaces. 
# See notes for theoretical details.
"""

	fndMax(f::Any, lb::Float64, ub::Float64, wdt::Float64, prc::Float64)::Float64

# Naively locate the maximum of a function within a specified interval.
"""
function fndMax(f::Any, lb::Float64, ub::Float64, wdt::Float64, prc::Float64)::Tuple{Float64, Float64}

	lWdt = wdt
	r = lb : wdt : ub
	(val, ind) = findmax(f.(r))

	while lWdt > prc
		
		lWdt *= 0.1
		
		if ind > 1 && ind < length(r)

			r = r[ind - 1] : lWdt : r[ind + 1]
		elseif ind == 1

			r = r[1] : lWdt : r[2]
		else

			r = r[end - 1] : lWdt : r[end]
		end

		(val, ind) = findmax(f.(r))
	end

	return (val, r[ind])
end
"""

	fndWdt(f::Any, sLoc::Float64, sVal::Float64, initWdt::Float64, relV::Float64)::Float64

# Naively determine width of nearly singular behavior.
"""
function fndWdt(f::Any, sLoc::Float64, sVal::Float64, initWdt::Float64, relV::Float64)::Float64

	lrelV = 1.0
	lWdt = 2.0 * initWdt

	while lrelV < relV
		
		lWdt *= 0.5	
		lrelV = /(f(sLoc + lWdt), sVal)
	end

	return lWdt
end
"""

	geoSrs(intTrm::ComplexF64, geoFct::ComplexF64)::ComplexF64

Analytic extension of the geometric series for a given initial term and geometric factor.
"""
function geoSrs(intTrm::ComplexF64, geoFct::ComplexF64)::ComplexF64

	return /(intTrm, 1.0 - geoFct)
end
"""

	zSqrt(val::ComplexF64)::ComplexF64

Square root defined with cut line along real axis. 
"""
function zSqrt(val::ComplexF64)::ComplexF64

	if angle(val) >= 0.0

		return sqrt(val)
	else

		return sqrt(val) * exp(im * pi)
	end
end
"""

	wvcPrp(rsp::ComplexF64, lTck::Float64, wvc::Float64)::ComplexF64

Perpendicular wave, incorporating effective normalization.
"""
function wvcPrp(rsp::ComplexF64, lTck::Float64, wvc::Float64)::ComplexF64

	return zSqrt(rsp - wvc^2) + im * lcReg * ^(/(sqrt(abs(zSqrt(rsp - wvc^2))^2 + wvc^2) * lTck, mxSz), 4)
end
"""

	plcFnc(enr::Float64, tmp::Float64)::Float64

Planck function for energy in electron volts and temperature in Kelvin.
"""
function plcFnc(enr::Float64, tmp::Float64)::Float64

	return /(enr, exp(/(enr, blzK * tmp)) - 1.0)
end
"""
	flxPfc(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, enr::Float64)::Float

Energy prefactor in Watts per electron volts cm^-2 for fluctuating electric field intensity 
integrand. 
"""
function flxPfc(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, enr::Float64)::Float64

	return /(evJ * enr^2 * pi^3, 2.0^2 * muEv^2 * hBEv * ^(10.0, -8)) * (plcFnc(enr, lVar.tmpLst[lPair[1]]) - plcFnc(enr, lVar.tmpLst[lPair[2]]))
end
"""

	prpLyr(dst::Float64, rsp::ComplexF64, wvc::Float64)::ComplexF64

Calculation of phase accumulation under propagation, distance assumed to be in units of wavelength.
"""
function prpLyr(dst::Float64, rsp::ComplexF64, wvc::Float64)::ComplexF64

	return exp(2.0 * pi * im * wvcPrp(rsp, dst, wvc) * dst)
end
"""
	lyrTck(lVar::lyrDsc, lNum::Int)::Float64

Layer thickness in microns. 
"""
function lyrTck(lVar::lyrDsc, lNum::Int)::Float64

	if lNum == 1 || lNum == length(lVar.tmpLst)

		return 0.0
	else

		return lVar.bdrLoc[lNum] - lVar.bdrLoc[lNum - 1]
	end
end
"""

	lyrTckRel(lVar::lyrDsc, lNum::Int, enr::Float64)::Float64

Compute layer thickness in units of relative wavelength.
"""
function lyrTckRel(lVar::lyrDsc, lNum::Int, enr::Float64)::Float64

	return /(enr * lyrTck(lVar, lNum), muEv)
end
"""

	strTck(lVar::lyrDsc, enr::Float64)::Float64

Compute total structure thickness in units of relative wavelength.
"""
function strTck(lVar::lyrDsc, enr::Float64)::Float64

	return /(enr * (lVar.bdrLoc[end] - lVar.bdrLoc[1]), muEv)
end
"""

	xfrBdr(kPrp::NTuple{2,ComplexF64}, rsp::NTuple{2,ComplexF64})::Array{ComplexF64,1}

Isotopic electric reflection and transmission coefficients for 
a interface between two half spaces. 
First entry in the rsp permittivity tuple is consider to be 
incident. 
The first two returned entries are for 'p' polarized waves, 
the next two for 's' polarized waves. 
Odd numbers are transmission, even modes are reflection.
"""
function xfrBdr(wvc::Float64, kPrp::NTuple{2,ComplexF64}, rsp::NTuple{2,ComplexF64})::Array{ComplexF64,1}

	# Fresnel coefficients.
	return [/(2.0 * zSqrt(rsp[1] * rsp[2]) * kPrp[1], rsp[1] * kPrp[2] + rsp[2] * kPrp[1]), /(rsp[1] * kPrp[2] - rsp[2] * kPrp[1], rsp[1] * kPrp[2] + rsp[2] * kPrp[1]), /(2.0 * kPrp[1], kPrp[2] + kPrp[1]), /(kPrp[2] - kPrp[1], kPrp[2] + kPrp[1])]
end
"""

	bdrFill!(lVar::lyrDsc, enr::Float64, wvc::Float64, mode::Int, bdrArr::Array{ComplexF64,2})::Nothing

Fills bdrArr all individual boundary transfer coefficients, 
i.e. as though every boundary were between half-spaces. 
"""
# If mode == 1 smaller indices are treated as initial. 
# If mode == 2 larger indices are treated as initial.
# Convention for bdrArr as described for xfrBdr. 
function bdrFill!(lVar::lyrDsc, enr::Float64, wvc::Float64, mode::Int, bdrArr::Array{ComplexF64,2})::Nothing

	if mode == 1

		for bdr = 1:length(lVar.bdrLoc)

			bdrArr[:,bdr] .= xfrBdr(wvc, (wvcPrp(lVar.rspPrf[bdr](enr), lyrTckRel(lVar, bdr, enr), wvc), wvcPrp(lVar.rspPrf[bdr + 1](enr), lyrTckRel(lVar, bdr + 1, enr), wvc)), (lVar.rspPrf[bdr](enr), lVar.rspPrf[bdr + 1](enr)))
		end
	elseif mode == 2

		for bdr = 1:length(lVar.bdrLoc)
			
			bdrArr[:,bdr] .= xfrBdr(wvc, (wvcPrp(lVar.rspPrf[bdr + 1](enr), lyrTckRel(lVar, bdr + 1, enr), wvc), wvcPrp(lVar.rspPrf[bdr](enr), lyrTckRel(lVar, bdr, enr), wvc)), (lVar.rspPrf[bdr + 1](enr), lVar.rspPrf[bdr](enr)))
		end
	else

	 error("Unrecognized fill mode.")	
	end

	return nothing 
end
"""
	
	xfrFlmExp!(lVar::lyrDsc, enr::Float64, wvc::Float64, lNum::Int, mode::Int, refArrL::SubArray{ComplexF64,2}, refArrR::SubArray{ComplexF64,2}, xfrArr::SubArray{ComplexF64,1})::Nothing

Updates transfer coefficients under expansion of nested reflection coefficient. 
See associated notes for additional details. 
"""
# If mode == 1, the increment is ``added on the right'', 
# i.e. building from the left. 
# If mode == 2, the increment is ``added on the left'', i.e. 
# building from the right.
# xfrArr holds the expanded and nested coefficients of the 
# numerator as the first entries of a given column. 
# The next two entries are the expanded and nested 
# coefficients of the denominator. 
# See notes for additional details. 
function xfrFlmInc!(lVar::lyrDsc, enr::Float64, wvc::Float64, lNum::Int, mode::Int, refArrL::SubArray{ComplexF64,2}, refArrR::SubArray{ComplexF64,2}, xfrArr::SubArray{ComplexF64,1})::Nothing

	if mode == 1
	
		copyto!(xfrArr, [xfrArr[1] + refArrR[lNum - 1] * xfrArr[2], prpLyr(2.0 * lyrTckRel(lVar, lNum - 1, enr), lVar.rspPrf[lNum - 1](enr), wvc) * (xfrArr[2] - refArrL[lNum - 1] * xfrArr[1]), xfrArr[3] + refArrR[lNum - 1] * xfrArr[4], prpLyr(2.0 * lyrTckRel(lVar, lNum - 1, enr), lVar.rspPrf[lNum - 1](enr), wvc) * (xfrArr[4] - refArrL[lNum - 1] * xfrArr[3])])

	elseif mode == 2
		
		copyto!(xfrArr, [xfrArr[1] + refArrL[lNum] * xfrArr[2], prpLyr(2.0 * lyrTckRel(lVar, lNum - 1, enr), lVar.rspPrf[lNum - 1](enr), wvc) * (xfrArr[2] - refArrR[lNum] * xfrArr[1]), xfrArr[3] + refArrL[lNum] * xfrArr[4], prpLyr(2.0 * lyrTckRel(lVar, lNum - 1, enr), lVar.rspPrf[lNum - 1](enr), wvc) * (xfrArr[4] - refArrR[lNum] * xfrArr[3])])
	else 

		error("Unrecognized increment mode.")	
	end

	return nothing
end
"""
	
	xfrFlmNst!(lVar::lyrDsc, enr::Float64, wvc::Float64, lNum::Int, mode::Int, refArrL::SubArray{ComplexF64,2}, refArrR::SubArray{ComplexF64,2}, xfrArr::SubArray{ComplexF64,1})::Nothing

Updates transfer coefficients by nesting current values, i.e. inclusion of an additional layer. 
See associated notes for additional details. 
"""
# If mode == 1, the increment is ``added on the right'', 
# i.e. building from the left. 
# If mode == 2, the increment is ``added on the left'', i.e. 
# building from the right.
# xfrArr holds the expanded and nested coefficients of the 
# numerator as the first entries of a given column. 
# The next two entries are the expanded and nested 
# coefficients of the denominator. 
# See notes for additional details. 
function xfrFlmInc!(lVar::lyrDsc, enr::Float64, wvc::Float64, lNum::Int, mode::Int, refArrL::SubArray{ComplexF64,2}, refArrR::SubArray{ComplexF64,2}, xfrArr::SubArray{ComplexF64,1})::Nothing

	if mode == 1
	
		copyto!(xfrArr, [xfrArr[1] + refArrR[lNum - 1] * xfrArr[2], prpLyr(2.0 * lyrTckRel(lVar, lNum - 1, enr), lVar.rspPrf[lNum - 1](enr), wvc) * (xfrArr[2] - refArrL[lNum - 1] * xfrArr[1]), xfrArr[3] + refArrR[lNum - 1] * xfrArr[4], prpLyr(2.0 * lyrTckRel(lVar, lNum - 1, enr), lVar.rspPrf[lNum - 1](enr), wvc) * (xfrArr[4] - refArrL[lNum - 1] * xfrArr[3])])

	elseif mode == 2
		
		copyto!(xfrArr, [xfrArr[1] + refArrL[lNum] * xfrArr[2], prpLyr(2.0 * lyrTckRel(lVar, lNum - 1, enr), lVar.rspPrf[lNum - 1](enr), wvc) * (xfrArr[2] - refArrR[lNum] * xfrArr[1]), xfrArr[3] + refArrL[lNum] * xfrArr[4], prpLyr(2.0 * lyrTckRel(lVar, lNum - 1, enr), lVar.rspPrf[lNum - 1](enr), wvc) * (xfrArr[4] - refArrR[lNum] * xfrArr[3])])
	else 

		error("Unrecognized increment mode.")	
	end

	return nothing
end
"""

	relFlm(lVar::lyrDsc, enr::Float64, wvc::Float64, lPair::Tuple{Int64,Int64}, trfCff::Array{ComplexF64,1})::Array{ComplexF64,1}

Relative field amplitude coefficients from transfer 
coefficients. 
See tfrFlm for structure for transfer coefficient array. 
"""
function relFlm(lVar::lyrDsc, enr::Float64, wvc::Float64, lPair::Tuple{Int64,Int64}, trfCff::Array{ComplexF64,1})::Array{ComplexF64,1}
	
	# Layer phase propagation 
	if lPair[1] > 1

		lLyrPhz = prpLyr(lyrTckRel(lVar, lPair[1], enr), lVar.rspPrf[lPair[1]](enr), wvc)
	else

		lLyrPhz = 0.0
	end

	if lPair[2] < length(lVar.tmpLst)

		rLyrPhz = prpLyr(lyrTckRel(lVar, lPair[2], enr), lVar.rspPrf[lPair[2]](enr), wvc)
	else

		rLyrPhz = 0.0
	end
	# Field dressing, absolute value used for stability, consistent with later expressions
	pDrs = abs(geoSrs(1.0 + im * 0.0, trfCff[3] * lLyrPhz * trfCff[5] * lLyrPhz)) * abs(geoSrs(1.0 + im * 0.0, trfCff[7] * rLyrPhz * trfCff[9] * rLyrPhz))
	sDrs = abs(geoSrs(1.0 + im * 0.0, trfCff[4] * lLyrPhz * trfCff[6] * lLyrPhz)) * abs(geoSrs(1.0 + im * 0.0, trfCff[8] * rLyrPhz * trfCff[10] * rLyrPhz))
	# Propagation term that should be included in left moving waves is shifted to the 
	# flux integrand to avoid overflow. 
	return [pDrs * trfCff[1], trfCff[9] * trfCff[1], sDrs * trfCff[2], sDrs * trfCff[10] * trfCff[2]]
end
"""
	
	tfrFlm!(lVar::lyrDsc, enr::Float64, wvc::Float64, lPairs::Array{Int64,2}, bdrArrL::Array{ComplexF64,2}, bdrArrR::Array{ComplexF64,2}, trfCff::Array{ComplexF64,2})::Nothing

Calculates field transfer coefficients for a list of ordered 
pairs of layers. 
"""
# bdrArrL should contain interface transfer coefficients with 
# smaller layer numbers treated as initial. 
# bdrArrR should contain interface transfer coefficients with 
# larger layer numbers treated as initial. 
# trfCff[:, lyr] = [trnCf; refLL; refLR; refRL; refRR]
# p pol followed by s pol
function tfrFlm!(lVar::lyrDsc, enr::Float64, wvc::Float64, lPairs::Array{Int64,2}, bdrArrL::Array{ComplexF64,2}, bdrArrR::Array{ComplexF64,2}, trfCff::Array{ComplexF64,2})::Nothing

	# Double check that pairs in list are correctly ordered. 
	for ind = 1:size(lPairs)[2]

		if lPairs[ind, 2] <= lPair[ind, 1]
			
			error("Program convention requires ordered layers.")
		return nothing
	end

	srcLyr = sort!(unique(lPairs[1,:]))
	trgLyr = sort!(unique(lPairs[2,:]))
	# Number of layers
	lyrN = length(lVarL.tmpLst)

	# Preallocate
	# Coefficient for calculating reflection. 
	# Asymmetry comes from assumption of ordered layers. 
	xfrCffLP = fill(0.0 + 0.0im, 4, length(trgLyr))
	xfrCffLS = fill(0.0 + 0.0im, 4, length(trgLyr))
	xfrCffRP = fill(0.0 + 0.0im, 4)
	xfrCffRS = fill(0.0 + 0.0im, 4)
	

	### Calculate layer reflection and transmission coefficients
	## Right incident reflection coefficients. 
	# Treat possibility that source layer is first layer.
	if srcLyr[1] == 1

		for lyr in findall(x -> x == 1, lPairs[1,:])

			trfCff[3:4,lyr] .= [bdrArrR[2,1], bdrArrR[4,1]]
		end
	end 

	# Left boundary seed.
	xfrCffLP 

	ind = 2
	
	while ind < rLyr 

		if lLyr == ind
			# Save left reflection coefficients for left cell.
			refLL .= [xfrRef[2], xfrRef[4]]
		end	
		# Increment reflection field coefficients.
		ind = ind + 1
		xfrFlmInc!(lVarL, enr, wvc, ind, 1, xfrRef)
	end
	# Save reflection coefficient for right cell
	refRL .= [xfrRef[2], xfrRef[4]]

	## Right coefficients coefficients
	# Right boundary 
	ind = lyrN - 1
	# Reflection seed
	xfrRef .= xfrBdr(wvc, (wvcPrp(lVarL.rspPrf[lyrN - 1](enr), lyrTckRel(lVar, lyrN - 1, enr), wvc), wvcPrp(lVarL.rspPrf[lyrN](enr), lyrTckRel(lVar, lyrN, enr), wvc)), (lVarL.rspPrf[lyrN - 1](enr), lVarL.rspPrf[lyrN](enr)))
	# Transmission seed
	xfrTrn = xfrBdr(wvc, (wvcPrp(lVarL.rspPrf[rLyr - 1](enr), lyrTckRel(lVar, rLyr - 1, enr), wvc), wvcPrp(lVarL.rspPrf[rLyr](enr), lyrTckRel(lVar, rLyr, enr), wvc)), (lVarL.rspPrf[rLyr - 1](enr), lVarL.rspPrf[rLyr](enr)))

	while ind > lLyr

		if rLyr == ind 
			# Save right reflection coefficients for right cell.
			refRR .= [xfrRef[2], xfrRef[4]]
		end	
		# Increment reflection field coefficients.
		ind = ind - 1
		xfrFlmInc!(lVarL, enr, wvc, ind, 2, xfrRef)
		# Increment transmission field coefficients. 
		if ind < rLyr - 1
			xfrFlmInc!(lVarL, enr, wvc, ind, 2, xfrTrn)
		end
	end

	# Transmission coefficients between cell pair.
	trnCf .= [xfrTrn[1], xfrTrn[3]]
	# Reflection coefficients for right edge of left cell. 
	refLR .= [xfrRef[2], xfrRef[4]]
	
	return relFlm(lVarL, enr, wvc, (lLyr, rLyr), [trnCf; refLL; refLR; refRL; refRR])
end
"""
	flxIntFunc(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, enr::Float64, wvc::Float64, trfCff::Array{ComplexF64,1})::Float64

Inner function computing wave vector integrand for heat flux between a pair of layers. 
"""
function flxIntFunc(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, enr::Float64, wvc::Float64, trfCff::Array{ComplexF64,1})::Float64
	# Permittivity response. 
	sRsp = lVar.rspPrf[lPair[1]](enr)
	tRsp = lVar.rspPrf[lPair[2]](enr)
	# Perpendicular wave vectors.
	swv = wvcPrp(sRsp, lyrTckRel(lVar, lPair[1], enr), wvc)
	twv = wvcPrp(tRsp, lyrTckRel(lVar, lPair[2], enr), wvc)
	# Layer thickness factors.
	# Account for possible finite thickness of source layer.
	if (lPair[1] != 1) && (lPair[1] != length(lVar.tmpLst))

		sLyrFac = (1.0 - exp(-4.0 * pi * imag(swv) * lyrTckRel(lVar, lPair[1], enr)))		
	else

		sLyrFac = 1.0
	end
	# Account for possible finite thickness of target layer. 
	# If lPair[1] > lPair[2], field transfer coefficients are computed for a mirror system. 
	# Hence, only ``right'' moving coefficients, trfCff[1] and trfCff[3] need to be considered.
	if (lPair[2] != 1) && (lPair[2] != length(lVar.tmpLst))

		lTck = lyrTckRel(lVar, lPair[2], enr)
		# Single layer factor comes from shifting factors to avoid overflow.
		tLyrFacR = 1.0 - exp(-4.0 * pi * imag(twv) * lTck)
		tLyrFacL = exp(-4.0 * pi * imag(twv) * lTck) * (1.0 - exp(-4.0 * pi * imag(twv) * lTck))
		tLyrFacM = exp(-4.0 * pi * lTck * (im * real(twv) + imag(twv))) * (1.0 - exp(4.0 * pi * im * real(twv) * lTck))
	else
		
		tLyrFacR = 1.0 
		tLyrFacL = 0.0
		tLyrFacM = 0.0
	end
	
	# Integrand components.
	intFac = /(imag(sRsp) * imag(tRsp), 4.0 * abs(swv)^2)
	srcFac = /(sLyrFac, imag(swv))
	# Magnitude factor for p-polarized right moving waves.
	pFacR = /(abs(wvc^2 + conj(swv) * twv)^2, abs(tRsp) * abs(sRsp))
	# Magnitude factor for p-polarized left moving waves.
	pFacL = /(abs(wvc^2 - conj(swv) * twv)^2, abs(tRsp) * abs(sRsp))
	# Magnitude factor for p-polarized wave mixing.
	pFacM = /((wvc^2 + conj(swv) * twv) * (wvc^2 - swv * conj(twv)), abs(tRsp) * abs(sRsp))
	## Polarization contributions
	pPol = /(abs(trfCff[1])^2 * tLyrFacR * pFacR + abs(trfCff[2])^2 * tLyrFacL * pFacL, imag(twv)) + /(2.0 * imag(trfCff[1] * conj(trfCff[2]) * tLyrFacM * pFacM), real(twv)) 	
	
	sPol = /(abs(trfCff[3])^2 * tLyrFacR + abs(trfCff[4])^2 * tLyrFacL, imag(twv)) + /(2.0 * imag(trfCff[3] * conj(trfCff[4]) * tLyrFacM), real(twv))

	# Correct possible NaN underflow.
	if isnan(intFac * srcFac * (pPol + sPol))

		return 0.0
	else
		# Impose theoretical transfer cut off in case of numerical inaccuracy.
		return wvc * (min(1.0, intFac * srcFac * pPol) + min(1.0, intFac * srcFac * sPol))
	end
end
"""
	evaVarED

Variable transformation for decaying waves. 
"""
function evaVarED(dlt::Float64, uvc::Float64)::Float64

	return tan(/(pi, 2.0) * uvc - dlt)
end
"""
	
	flxFncED(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, enr::Float64, uvc::Float64)::Float64

Integrand for evanescent tail of heat flux. 
"""
function flxFncED(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, dlt::Float64, enr::Float64, uvc::Float64)::Float64

	# Evanescent integrand for quadrature integration. 
	return /(flxIntFunc(lVar, lPair, enr, evaVarED(intDltInf, uvc), tfrFlm(lVar, enr, evaVarED(intDltInf, uvc), lPair)), cos(/(pi, 2.0) * uvc - intDltInf)^2)
end
"""
	evaVarEILD(cen::Float64, dlt::Float64, uvc::Float64)::Float64

Double Lorentzian variable transformation for intermediate waves. 
"""
function evaVarEILD(cen::Float64, dlt::Float64, uvc::Float64)::Float64

	return sqrt(/(1.0, sqrt(uvc)) - ^(dlt, 2)) + cen
end
"""
	
	flxFncEILD(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, enr::Float64, shft::Float64, lwb::Float64, uvc::Float64)::Float64

Double Lorentzian scaled integrand for ``intermediate waves''. 
"""
function flxFncEILD(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, cen::Float64, wdt::Float64, dlt::Float64, enr::Float64, uvc::Float64)::Float64

	# Intermediate field integrand for quadrature integration. 
	return /(flxIntFunc(lVar, lPair, enr, evaVarEILD(cen, dlt, uvc), tfrFlm(lVar, enr, evaVarEILD(cen, dlt, uvc), lPair)) + flxIntFunc(lVar, lPair, enr, evaVarEILD(cen, dlt, ^(wdt^2 + dlt^2, -2)) - evaVarEILD(cen, dlt, uvc), tfrFlm(lVar, enr, evaVarEILD(cen, dlt, ^(wdt^2 + dlt^2, -2)) - evaVarEILD(cen, dlt, uvc), lPair)), 4.0 * uvc * sqrt(sqrt(uvc) - ^(dlt, 2)))
end
"""
	
	flxFncPR(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, enr::Float64, uvc::Float64)::Float64

Integrand for propagating portion of heat flux. 
"""
function flxFncPR(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, enr::Float64, uvc::Float64)::Float64

	# Propagating field integrand for quadrature integration.
	return flxIntFunc(lVar, lPair, enr, uvc, tfrFlm(lVar, enr, uvc, lPair))
end
"""

	function iwLocs(lVar::lyrDsc, enr::Float64)::Array{Float64,2}

Determine location and spacing for intermediate wave transformations. 
Note that output may have NaN entries, and function relies on global constants. 
"""
function iwLocs(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, enr::Float64)::Array{Float64,2}
	
	iwLcs = Float64[]
	
	# Find possible points needing more precise treatment.
	for rsp in unique(lVar.rspPrf)

		if /(real(rsp(enr)), imag(rsp(enr))) > iTrnV && imag(rsp(enr)) > 0

			push!(iwLcs, sqrt(real(rsp(enr)))) 	
		end
	end

	sort!(iwLcs)
	# Determine intermediate wave widths.
	iwWdt = Array{Float64}(undef, length(iwLcs))

	quadFlxPR(evlPt) = flxFncPR(lVar, lPair, enr, evlPt)

	for ind in 1:length(iwWdt)

		iwWdt[ind] = fndWdt(quadFlxPR, iwLcs[ind], quadFlxPR(iwLcs[ind]), lcWdt, lcRel)
	end
	# Remove locations that are too close together. 
	curr = 1
	
	for ind in 1:(length(iwWdt) - 1)
		
		if ind == 1

			if iwLcs[ind] - iwWdt[ind] < 0.0

				iwWdt[ind] = iwLcs[ind]
			end
		end

		if iwLcs[ind + 1] - iwLcs[curr] < /(iwWdt[curr], 10.0)

			iwLcs[ind + 1] = NaN
			iwWdt[ind + 1] = NaN
		else 

			if /(iwLcs[ind + 1] - iwLcs[curr], 2.0) < iwWdt[curr] 

				iwLcs[curr] = /(iwLcs[ind + 1] - iwLcs[curr], 3.0)
			end

			if /(iwLcs[ind + 1] - iwLcs[curr], 2.0) < iwWdt[ind + 1]

				iwLcs[curr] = /(iwLcs[ind + 1] - iwLcs[curr], 3.0)
			end

			curr = ind + 1
		end
	end
	# Order remaining locations
	iwDsc = sort([iwLcs iwWdt], dims = 1)
	return iwDsc
end
"""
	
	femFlx(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, enr::Float64)::Float64

Thermal flux between a pair of layers per ev cm^2. 
"""
# mode == 0 only includes evanescent contribution beyond k/ko > 5. 
# mode == 1 calculates full heat flux, without variable transformation for intermediate waves.
# mode == 2 calculates full heat flux, with variable transformation for intermediate waves. 
function femFlx(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, enr::Float64, mode::Int)::Float64

	# Define single variable functions for quadrature integration.
	# Propagating waves
	quadFlxPR(evlPt) = flxFncPR(lVar, lPair, enr, evlPt)
	# # Rapidly decaying waves, intDltInf sets effective upper bound as tan(pi/2 - inDlt)
	quadFlxED(evlPt) = flxFncED(lVar, lPair, intDltInf, enr, evlPt)
	# # Treat intermediate waves. 
	if mode == 2
		# Determine intermediate wave locations.
		iwDsc = iwLocs(lVar, lPair, enr)
		# Account for possible NaN entries.
		nans = findfirst(isnan, view(iwDsc, :, 1))

		if !isnothing(nans)

			numIW = nans - 1
		else

			numIW = length(iwDsc[:,1])
		end
		# Compute integral, filling in portions not treated as intermediate waves with 
		# linear quadrature.
		if iwDsc[1, 1] - iwDsc[1, 2] > 1.0e-6
		
			intVal = quadgk(quadFlxPR, 0.0, iwDsc[1, 1] - iwDsc[1, 2], rtol = relTol, order = 16)[1]
		else

			intVal = 0.0
		end

		for ind in 1:numIW

			# Width parameter for double Lorentzian.
			dltIW = min(^(quadFlxPR(iwDsc[ind, 1]), - 0.5), /(iwDsc[ind, 2], 10.0))
			# Transformed flux integrand
			quadFlxEILD(evlPt) = flxFncEILD(lVar, lPair, iwDsc[ind, 1], iwDsc[ind, 2], dltIW, enr, evlPt)
			# Perform integrations
			intVal += quadgk(quadFlxEILD, ^(^(iwDsc[ind, 2], 2) + ^(dltIW, 2), -2), ^(dltIW, -4) , rtol = relTol, order = qudOrd)[1]

			if ind < numIW
				
				intVal += quadgk(quadFlxPR, iwDsc[ind, 1] + iwDsc[ind, 2], iwDsc[ind + 1, 1] - iwDsc[ind + 1, 2], rtol = relTol, order = qudOrd)[1]
			end	
		end
		# # Transition to decaying waves  
		if iwDsc[numIW, 1] + iwDsc[numIW, 2] < 10.0
		
			intVal += quadgk(quadFlxPR, iwDsc[numIW, 1] + iwDsc[numIW, 2], 10.0, rtol = relTol, order = qudOrd)[1]
			
			evaDecTrn = /(2.0 * atan(10.0), pi)

		else

			evaDecTrn = /(2.0 * atan(iwDsc[numIW, 1] + iwDsc[numIW, 2]), pi)
		end
	
	elseif mode == 1

		 intVal = quadgk(quadFlxPR, 0.0, 1.0, rtol = relTol, order = qudOrd)[1]
		 evaDecTrn = /(2.0 * atan(1.0), pi)
	
	elseif mode == 0
		
		intVal = 0.0
		evaDecTrn = /(2.0 * atan(5.0), pi)	
	else

		error("Unrecognized flux calculation mode.")
		return nothing
	end
	# Compute decaying wave contribution 
	intVal += quadgk(quadFlxED, evaDecTrn, 1.0, rtol = relTol, order = qudOrd)[1]

	return flxPfc(lVar, lPair, enr) * intVal
end
# end