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
	testFlxIntFunc(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, enr::Float64, wvc::Float64)::Float64

Known wavevector integrand for heat transfer across a slab, three layer, structure. 
"""
function testFlxIntFunc(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, enr::Float64, wvc::Float64)::Float64

	# Permittivity response. 
	sRsp = lVar.rspPrf[lPair[1]](enr)
	tRsp = lVar.rspPrf[lPair[2]](enr)
	# Gap thickness
	gTck = lyrTckRel(lVar, lPair[2] - 1, enr)
	# Perpendicular wave vectors.
	sWV = wvcPrp(lVar.rspPrf[lPair[1]](enr), lyrTckRel(lVar, lPair[1], enr), wvc)
	tWV = wvcPrp(lVar.rspPrf[lPair[2]](enr), lyrTckRel(lVar, lPair[2], enr), wvc)
	gWV = wvcPrp(lVar.rspPrf[lPair[2] - 1](enr), gTck, wvc)
	# Gap propagation.
	gPrp = exp(4.0 * pi * im * gWV * gTck)
	# Fresnel coefficients. 
	xfrC = xfrBdr(wvc, (gWV, tWV), (lVar.rspPrf[lPair[2] - 1](enr), lVar.rspPrf[lPair[2]](enr)))

	if wvc > 1.0
		
		pT = /(4.0 * abs(gPrp) * imag(xfrC[2])^2, abs(1.0 - xfrC[2]^2 * gPrp)^2)
		sT = /(4.0 * abs(gPrp) * imag(xfrC[4])^2, abs(1.0 - xfrC[4]^2 * gPrp)^2)
	else
		
		pT = /((1.0 - abs(xfrC[2])^2)^2, abs(1.0 - xfrC[2]^2 * gPrp)^2)
		sT = /((1.0 - abs(xfrC[4])^2)^2, abs(1.0 - xfrC[4]^2 * gPrp)^2)
	end

	return wvc * (pT + sT)
end
"""
	evaVarEIAT

ArcTan variable transformation for intermediate waves. 
"""
function evaVarEIAT(shft::Float64, lwb::Float64, dlt::Float64, uvc::Float64)::Float64

	return /(shft * (atan(uvc) + /(pi, 2.0) - dlt), /(pi, 2.0) - dlt) + lwb
end
"""
	
	flxFncEIAT(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, shft::Float64, lwb::Float64, dlt::Float64, enr::Float64, uvc::Float64)::Float64

ArcTan scaled integrand for ``intermediate waves''. 
"""
function flxFncEIAT(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, shft::Float64, lwb::Float64, dlt::Float64, enr::Float64, uvc::Float64)::Float64

	# Intermediate field integrand for quadrature integration. 
	return /(shft * flxIntFunc(lVar, lPair, enr, evaVarEIAT(shft, lwb, dlt, uvc), tfrFlm(lVar, enr, evaVarEIAT(shft, lwb, dlt, uvc), lPair)), (/(pi, 2.0) - dlt) * (1.0 + ^(uvc, 2)))
end
"""
	prfIntFunc(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, enr::Float64, wvc::Float64, prfArr::Array{Float64,1}, trfCff::Array{ComplexF64,1})::Array{Float64,1}

Inner function computing wave vector integrand for heat flux profile between a pair of layers. 
"""
function prfIntFunc(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, enr::Float64, wvc::Float64, prfArr::Array{Float64,1}, trfCff::Array{ComplexF64,1})::Array{Float64,1}

	# Perpendicular wave vectors.
	swv = wvcPrp(lVar.rspPrf[lPair[1]](enr), lyrTckRel(lVar, lPair[1], enr), wvc)
	twv = wvcPrp(lVar.rspPrf[lPair[2]](enr), lyrTckRel(lVar, lPair[2], enr), wvc)

	# Layer thickness relative to wavelength.
	lTck = lyrTckRel(lVar, lPair[2], enr)
	# Account for possible finite thickness of source layer.
	if (lPair[1] != 1) && (lPair[1] != length(lVar.tmpLst))

		sLyrFac = (1.0 - exp(-4.0 * pi * imag(swv) * lyrTckRel(lVar, lPair[1], enr)))		
	else

		sLyrFac = 1.0
	end
	# Account for possible finite thickness of target layer. 
	# If lPair[1] > lPair[2], field transfer coefficients are computed for a mirror system. 
	# Hence only ``right'' moving coefficients, trfCff[1] and trfCff[3] need to be considered.
	if (lPair[2] != 1) && (lPair[2] != length(lVar.tmpLst))

		tPrf = lTck .* prfArr
		# Note that phase is shifted for left moving waves to avoid possible overflow. 
		tLyrFacL = @. exp(4.0 * pi * imag(twv) * (tPrf - 2.0 * lTck))
		tLyrFacM = @. exp(4.0 * pi * im * (real(twv) * tPrf + twv * lTck))
	# For half space layers, the relative thickness in prfArr is taken to refer to 
	# the separation between the two outermost boundaries.
	else
		
		tPrf = strTck(lVar, enr) .* prfArr
		tLyrFacL = zeros(length(prfArr))
		tLyrFacM = zeros(length(prfArr))
	end
	
	tLyrFacR = @. exp(-4.0 * pi * imag(twv) * tPrf)

	# Integrand components.
	intFac = /(imag(sRsp) * imag(tRsp) * wvc, 4.0 * abs(swv)^2)
	srcFac = /(sLyrFac, imag(swv))
	# Magnitude factor for p-polarized right moving waves.
	pRghFac = /(abs(wvc^2 + conj(swv) * twv)^2 , abs(tRsp) * abs(sRsp))
	# Magnitude factor for p-polarized left moving waves.
	pLftFac = /(abs(wvc^2 - conj(swv) * twv)^2 , abs(tRsp) * abs(sRsp))
	# Magnitude factor for p-polarized wave mixing.
	pMixFac = /((wvc^2 + conj(swv) * twv) * (wvc^2 - swv * conj(twv)), abs(tRsp) * abs(sRsp))

	pPol = @. abs(trfCff[1])^2 * tLyrFacR * pRghFac + abs(trfCff[2])^2 * tLyrFacL * pLftFac + 2.0 * real(trfCff[1] * conj(trfCff[2]) * tLyrFacM * pMixFac)

	sPol = @. abs(trfCff[3])^2 * tLyrFacR + abs(trfCff[4])^2 * tLyrFacL + 2.0 * real(trfCff[3] * conj(trfCff[4]) * tLyrFacM)

	intEPrf = @. intFac * srcFac * (pPol + sPol)

	# Correct possible not a number underflow.
	for ind = 1:length(intEPrf)

		if isnan(intEPrf[ind])

			intEPrf[ind] = 0.0
		end
	end

	return intEPrf
end
# NOT CURRENTLY CONSISTENT WITH REST OF CODE, NEEDS UPDATE
"""
	femPrf(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, enr::Float64,
	prfArr::Array{Float64})::Array{Float64}

Profile of absorption between a pair of layers per ev cm^2. 
prfArr should consist of numbers between 0 and 1, representing relative thicknesses 
within the layer. 
"""
function femPrf(lVar::lyrDsc, lPair::Tuple{Int64,Int64}, enr::Float64, prfArr::Array{Float64,1})::Array{Float64,1}

	quadPrfFnc(evlPt) = prfFnc(lVar, lPair, enr, evlPt, prfArr)
	
	# Specify possible poles. 
	if real(lVar.rspPrf[lPair[1]](enr)) > 0

		polLoc = /(2.0, pi) * atan(sqrt(abs(real(lVar.rspPrf[lPair[1]](enr)))))
	else 

		polLoc = 0.5
	end
	# Middle values is pole location, end values are interval.
	return flxPfc(lVar, lPair, enr) .* (quadgk(quadPrfFnc, 0.0, polLoc, 1.0 - intDltInf, rtol = relTol)[1])
end