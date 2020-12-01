using Base.Threads, ProgressMeter, Plots, LaTeXStrings
using ResponseModels, FilmUtilities
# Radial frequency to electron volts.
const hBEv = 6.582e-16
## Program settings
# Description of the layer structure.
# Boundary locations in units of microns.
bdrLoc = [0.000, 0.001, 0.101, 0.102, 0.202, 0.203]
# Layer temperatures in units of Kelvin.
tmpLst = [300.0, 300.0, 0.0, 0.0, 0.0, 0.0, 0.0] 
## Construct layer structure description container. 
# Permittivity model for gap.
gRsp(enr) = cstRsp(1.0 + im * ^(10.0, -7), enr)
# inAsRsp(enr) = cstRsp(13.0 + im * ^(10.0, -7), enr) + im * 0.25 * hStep(enr - 0.43)
# qtSCRsp(enr) = cstRsp(14.0 + im * ^(10.0, -7), enr) + im * enr * 2.2
# Permittivity model for silicon. 
# Boron p-type doping
brnDpt = dptDsc(0.044, 1.0e18)
# Empty doping
nllDpt = dptDsc(0.0, 0.0)
# Calculate silicon model parameters
siModE = prmMSi(tmpLst[1], nllDpt, nllDpt)
siRspE(enr) = siRsp(enr, siModE)

# Construct description of structure. 
lVar = lyrDsc(bdrLoc, tmpLst, [gRsp, sicRsp, gRsp, sicRsp, gRsp, sicRsp, gRsp])

lPairs = Array{Int64,2}(undef, 2, 2)
lPairs[1,1] = 2
lPairs[2,1] = 4 
lPairs[1,2] = 2
lPairs[2,2] = 6 

enrPts = 0.06:0.0005:0.15
trgPts = Array{Float64,2}(undef, size(lPairs)[2], length(enrPts))
rspVals = Array{ComplexF64,1}(undef, length(enrPts))

tfrIntEv!(lVar, lPairs, enrPts, trgPts)

@threads for ind = 1:length(enrPts)

	rspVals[ind] = sicRsp(enrPts[ind])
end
# # Conversion for comparing with APS paper ./(hBEv .* ^(10.0, 4.0), 300) 
pgfplotsx()

plt1 = plot(enrPts, [real.(rspVals), imag(rspVals)], leg = false, xaxis = "Energy (eV)", yaxis = L"$\epsilon$")

plt2 = plot(enrPts, view(trgPts, 1, :), leg = false, xaxis = "Energy (eV)", yaxis = (L"Flux Log(W /(eV cm$\!^2$)", :log))

plt3 = plot(enrPts, view(trgPts, 2, :), leg = false, xaxis = "Energy (eV)", yaxis = (L"Flux Log(W /(eV cm$\!^2$)", :log))

plt = plot(plt1, plt2, plt3, layout = (1, 3))

display(plt)