using ProgressMeter, Plots, LaTeXStrings
using ResponseModels 
include("FilmUtilities.jl")
# Radial frequency to electron volts.
const hBEv = 6.582e-16
## Program settings
# Description of the layer structure.
# Boundary locations in units of microns.
bdrLoc = [0.000, 0.100]
# Layer temperatures in units of Kelvin.
tmpLst = [800.0, 0.0, 300.0] 

## Construct layer structure description container. 
# Permittivity model for gap.
gRsp(enr) = cstRsp(1.0 + im * ^(10.0, -7), enr)
inAsRsp(enr) = cstRsp(13.0 + im * ^(10.0, -7), enr) + im * 0.25 * hStep(enr - 0.43)
qtSCRsp(enr) = cstRsp(14.0 + im * ^(10.0, -7), enr) + im * enr * 2.2
# Permittivity model for silicon. 
# Boron p-type doping
brnDpt = dptDsc(0.044, 1.0e18)
# Empty doping
nllDpt = dptDsc(0.0, 0.0)
# Calculate silicon model parameters
siModE = prmMSi(tmpLst[1], nllDpt, nllDpt)
siRspE(enr) = siRsp(enr, siModE)

# Construct description of structure. 
lVar = lyrDsc(bdrLoc, tmpLst, [siRspE, gRsp, inAsRsp])

lPairs = Array{Int64,2}(undef, 2, 1)
lPairs[1,1] = 1
lPairs[2,1] = 3 

enrPts = 0.2:0.001:1.0
trgPts = Array{Float64,2}(undef, size(lPairs)[2], length(enrPts))

@time tfrIntEv!(lVar, lPairs, enrPts, trgPts)
@time tfrIntEv!(lVar, lPairs, enrPts, trgPts)
@time tfrIntEv!(lVar, lPairs, enrPts, trgPts)
# # List of energy values to compute in electron volts.
# enrLst = 0.01:0.005:1.0
# # Transfer pair.
# lPairB = (1,6)
# lPairC = (1,4)

# # Preallocate memory for output storage. 
# flxValsB = Array{Float64,1}(undef, length(enrLst))
# flxValsC = Array{Float64,1}(undef, length(enrLst))
# rspVals = Array{ComplexF64,1}(undef, length(enrLst))

# # Progress bar
# prog = Progress(length(enrLst), dt = 0.5, barglyphs = BarGlyphs("[=> ]"), barlen = 50, color = :green)

# # Last argument provided to femFlx is the calculation mode. 
# # mode == 0 only includes evanescent contribution beyond k/ko > 5. 
# # mode == 1 calculates full heat flux, without variable transformation for intermediate waves.
# # mode == 2 calculates full heat flux, with variable transformation for intermediate 
# # waves. 
# # Mode 2 should be most exact. 

# @threads for ind = 1:length(enrLst)
	
# 	flxValsB[ind] = femFlx(lVar, lPairB, enrLst[ind], 2)
# 	flxValsC[ind] = femFlx(lVar, lPairC, enrLst[ind], 2)
# 	next!(prog)
# end

# # Conversion for comparing with APS paper ./(hBEv .* ^(10.0, 4.0), 300) 
# pgfplotsx()

# # plt1 = plot(enrLst, [real.(rspVals), imag(rspVals)], leg = false, xaxis = "Energy (eV)", yaxis = L"$\epsilon$")

# plt1 = plot(enrLst, flxValsB, leg = false, xaxis = "Energy (eV)", yaxis = (L"Flux Log(W /(eV cm$\!^2$)", :log))

# plt2 = plot(enrLst, flxValsC, leg = false, xaxis = "Energy (eV)", yaxis = (L"Flux Log(W /(eV cm$\!^2$)", :log))


# plt = plot(plt1, plt2, layout = (1, 2))

# display(plt)