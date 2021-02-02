__precompile__()
module OttawaSlabCell
using FilmDataStructures, ResponseModels,eps_gold,eps_InAsSbP,eps_InAs_ntype
export uOttawaSlabs
# Code for initializing uOttawa slab structure 
include("uOttawaSlabStruct.jl")
# Creates uOttawaSlabs(tEmit::Float64, tBck::Float64, divNCell::Int, divPCell::Int)::Tuple{lyrDsc, Array{Int64,2}} 
precompile(uOttawaSlabs, (Float64, Float64, Int, Int))
end