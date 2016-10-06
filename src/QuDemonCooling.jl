module QuDemonCooling
using QuBase
import QuBase.AbstractQuVector
export cooling!

# for general use
function cooling!(state::AbstractQuVector,U::AbstractMatrix,t::Real,gamma::Real)
    dice = rand()
    if dice <= 0.5*(1+sin(gamma))
        heat!(state,U,t,gamma)
        return 0.5*(1+sin(gamma))
    else
        cool!(state,U,t,gamma)
        return 0.5*(1-sin(gamma))
    end
end

function heat!(state::AbstractQuVector,U::AbstractMatrix,t::Real,gamma::Real)
    next_state = 0.5*( coeffs(state)+im*exp(im*gamma)*U*coeffs(state) ) |> normalize!
    CQST = QuBase.similar_type(state)
    return CQST(next_state,bases(state))
end

function cool!(state::AbstractQuVector,U::AbstractMatrix,t::Real,gamma::Real)
    next_state = 0.5*( coeffs(state)-im*exp(im*gamma)*U*coeffs(state) ) |> normalize!
    CQST = QuBase.similar_type(state)
    return CQST(next_state,bases(state))
end

end # module
