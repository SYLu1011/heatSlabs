fdim = 2
rTol = 1.0e-8
lwBnd = [-1.0, -1.0]
upBnd = [1.0, 1.0]

function intFunc(dom::Array{Float64,2}, cod::Array{Float64,2})::Nothing

	@threads for ind = 1:size(cod)[2]

		cod[1,ind] = dom[1,ind]^2 * dom[2,ind]^2
		cod[2,ind] = dom[1,ind]^2 + dom[2,ind]^2
	end

	return nothing
end

function quad_tst1v(x::Array{Float64,2}, v::Array{Float64,1})::Nothing

    @threads for i = 1:length(v)
        v[i] = prod(cos.(x[:,i]))
    end

    return nothing
end