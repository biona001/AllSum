module AllSum

using LinearAlgebra

export solve

function solve(
    zscores::Vector{T},
    R::AbstractMatrix{T},
    lambda0::T,
    lambda1::T,
    lambda2::T,
    sample_sizes::Vector{T},
    max_iter::Int = 1000
    ) where T
    # check error
    p = length(zscores)
    p == size(R, 1) == size(R, 2) || error("Dimention of zscores does not match R")

    # initialize
    r = [tj / sqrt(nj - 2 + tj^2) for (nj, tj) in zip(sample_sizes, zscores)]
    beta = zeros(T, length(zscores))
    beta_tilde = zeros(T, length(zscores))

    for t in 1:max_iter
        # update beta_tilde
        mul!(beta_tilde, R, beta, -1.0, 0.0)
        for j in eachindex(beta_tilde)
            beta_tilde[j] += r[j] + R[j, j] * beta[j] 
        end
        # coordinate descent
        for j in eachindex(beta)
            numer = abs(beta_tilde[j] - lambda1)
            denom = 1 + 2lambda2
            if numer > sqrt(2lambda0 * denom)
                beta[j] = numer / denom * sign(beta_tilde[j])
            else
                beta[j] = 0
            end
        end
    end

    return beta
end

end # module AllSum
