# measurements.jl
#
# Defines estimators and provides measurements
using Statistics


function sample(spin_left, operator_list)
    M = length(operator_list) ÷ 2
    spin_prop = copy(spin_left)

    for o in operator_list[1:M] #propagate half the list only (to the middle)
        if issiteoperator(o) && !isdiagonal(o)
            spin_prop[o.i] ⊻= 1 #spinflip
        end
    end

    return spin_prop
end


magnetization(spin_prop) = mean(x->2x-1, spin_prop)

#  (-2,i) is an off-diagonal site operator h(sigma^+_i + sigma^-_i)
#  (-1,i) is a diagonal site operator h
#  (0,0) is the identity operator I - NOT USED IN THE PROJECTOR CASE
#  (i,j) is a diagonal bond operator J(sigma^z_i sigma^z_j + 1)

function energy_abs_zero(h, J, spin_prop, operator_list)
    m_b = mean(isbondoperator, operator_list)
    # m_s = mean(issiteoperator, operator_list)

    E = (h/2 - J)*m_b
    return E
end



function num_single_site_diag(operator_list)
    return mean(x -> issiteoperator(x) && isdiagonal(x), operator_list)
end
