# Construction of PV structure for University of Ottawa collaboration. 
# See slides for schematic, slab structure is Si, Gap, InAsSbP, InAs (N cell),
# InAsSbP (P cell), InAs. 
# All temperatures are in Kelvin, and all lengths are in microns. 
# Absorption functions are current setup to calculate the number of photons 
# absorbed as interband transitions. 
# Alter code under the ***Optical models and absorption functions*** heading 
# for alternative behavior. 
####
# Number of division of the N-type InAs layer. 
# Resolution of photon generation is thckInAs / divNCell.
# divNCell = 3
# Number of division of the P-type InAsSbP layer. 
# divPCell = 3
# Temperature of the emitter.
# tEmit = 800.0
# Background temperature.
# tBck = 300.0

const e = 1.602*^(10,-19)
const eV = 1.60218*^(10,-19)
const hbar = 1.05457*^(10,-34) #m2kg/s
const hbar_eV = 6.5821*^(10,-16) #eV s 
function uOttawaSlabs(tEmit::Float64, tBck::Float64, divNCell::Int, divPCell::Int)
	### Settings
	# Separation distances and thickness in uma
	distGap = 0.1
	# Total InAs slab.
	thckInAs = 0.5
	# Total InAsSbP slab.
	thckInAsSbP = 0.5
	# Quaternary thickness. thckQuat = 0.025
	# Protective layer thickness. 
	thckProt = 0.015
	### Derived parameters
	## Layer numbers. 
	resNCell = /(thckInAs, divNCell) #divNCell is the number
	resPCell = /(thckInAsSbP, divPCell)
	# Total number of layers.
	numLayers = 4 + divNCell + divPCell
	## Temperature list of the layers. 
	tmpLst = fill(tBck, numLayers)
	# Set temperature of the emitter.
	tmpLst[1] = tEmit
	## Boundary locations. 
	bdrLoc = Vector{Float64}(undef, numLayers - 1)
	bdrLoc[1] = 0.0
	bdrLoc[2] = distGap
	bdrLoc[3] = bdrLoc[2] + thckProt
	# Fill boundaries for N-type part of the cell.
	for ind = 1 : divNCell

		bdrLoc[3 + ind] = bdrLoc[2 + ind] + resNCell
	end
	# Fill boundaries for P-type part of the cell
	for ind = 1 : divPCell

		bdrLoc[3 + divNCell + ind] = bdrLoc[2 + divNCell + ind] + resPCell
	end
	## Optical models and absorption functions. 
	# Silicon emitter. 
	# Dopant densities, currently empty?  Francoeur set real part of p-doped Si to be same as real part of top layer. P7
	acpDpt = dptDsc(0.0, 0.0)
	dnrDpt = dptDsc(0.0, 0.0)
	# Calculate silicon model parameters, see siMod and ResponseModels.
	siModParamsE = prmMSi(tmpLst[1], dnrDpt, acpDpt)
	# Construct silicon response model. 
	siRspE(enr) = siRsp(enr, siModParamsE)
	siAbsE(enr) = imag(siRspE(enr))
	# Gap response.
	gRsp(enr) = cstRsp(1.0 + im * ^(10.0, -7), enr)
	gAbs(enr) = 0.0
	# InAs N-type (PV Cell) doping concentration and temperature
	InAs_param =eps_InAs_struct(10.0^23,300.0)
	eps_nInAs(enr) = eps_InAsntype(enr,InAs_param)
	eps_nInAs_IBimag(enr) = imag(epsIB(enr*eV/hbar,InAs_param.N0,InAs_param.T,InAs_param.E0_T_InAs_value))
	# InAsSbP
	InAsSbPstructure = InAsSbP_struct(0.25,0.25)
	eps_InAsSbP(enr) = eps_InAsSbP_xy(enr,InAsSbPstructure)
	eps_InAsSbP_imag(enr) = eps_InAsSbP_imag_xy(enr,InAsSbPstructure)
	eps_gold_imag(enr) = imag(epsgold(enr))
	## Generate lists of optical responses and transfer factors. 
	optRsp = []
	trfFacs = []  #absorption
	# Layers prior to N-type part of the cell.
	push!(optRsp, siRspE, gRsp,eps_InAsSbP ) #Si, air and protective InAsSbP
	push!(trfFacs, siAbsE, gAbs,eps_InAsSbP_imag) 
	# N-type part of the PV cell. n-InAs
	for ind = 1 : divNCell
		
		push!(optRsp, eps_nInAs)
		push!(trfFacs, eps_nInAs_IBimag)
	end
	# P-type part of the PV cell. InAsSbP
	for ind = 1 : divPCell
		
		push!(optRsp, eps_InAsSbP)
		push!(trfFacs,eps_InAsSbP_imag)
	end
	# Optically thick backing, gold/p-InAs
	push!(optRsp, epsgold)
	push!(trfFacs,eps_gold_imag)
	## Set which layers transmission should be calculated for. 
	lPairs = Array{Int64,2}(undef, 2, divNCell + divPCell)
	# N-type layer pairs.
	for ind = 1 : divNCell

		lPairs[1, ind] = 1
		lPairs[2, ind] = 3 + ind
	end
	# N-type layer pairs.
	for ind = 1:divPCell

		lPairs[1, ind + divNCell] = 1
		lPairs[2, ind + divNCell] = 3 + divNCell + ind
	end
	# Build layer description for heat transfer code. 
	return (lyrDsc(bdrLoc, tmpLst, optRsp, trfFacs),lPairs)
end
precompile(uOttawaSlabs, (Float64, Float64, Int, Int))