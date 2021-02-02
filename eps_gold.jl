__precompile__()
module eps_gold
export epsgold
@inline function epsgold(e_photon::Float64)::ComplexF64
   # eVtoOmega = eV/hbar_eV
    eps_0 = 9.84
    wp = 9.010 #eV
    gamma = 0.072 #eV
    eps_di = eps_0 - wp^2/(e_photon^2 + gamma^2) + im*(wp^2*gamma/(e_photon*(e_photon^2 + gamma^2)))
    wc =2.4#eV (the second energy gap of the interbad transitions).
    A = 5.6#eV
    delta = 0.17 #eV
    delta_eps = A/(1  + exp(-(e_photon-wc)/delta))
    
    if(e_photon>1.8)
        eps = eps_di + 1j*delta_eps
    else
        eps =eps_di
    end
   # if imag(eps)<0
   #     return real(eps)-im*imag(eps)
   # else 
    return eps
end
precompile(epsgold,(Float64, ))
end

