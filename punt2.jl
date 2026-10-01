#1 sweep = N step normali di metropolis (faccio per eguagliare i tempi)
function metropolis_sweep!(J, β, σ)
    N = length(σ)
    for _ in 1:N
        i = rand(1:N) 
        
        
        sumsigma = 0
        for idx in J.colptr[i]:(J.colptr[i+1]-1) #posizione dei vicini dello spin i
            j = J.rowval[idx] #vicini di i 
            sumsigma += σ[j]
        end
        
        
        ΔE = 2 * σ[i] * sumsigma
        if ΔE <= 0 || rand() < exp(-β * ΔE)
            σ[i] *= -1 
        end
    end
end


function wolffsigma!(J, β, σ, nsteps; stats = (t, σ, Q)->nothing)
    n = size(J, 1)
    R = fill(false, n)
    Q = Int[]
    p = 1 - exp(-2β)

    for t in 1:nsteps
        # select initial spin
        push!(Q, rand(1:n))
        R[Q[1]] = true
        # grow Q
        k = 1
        while k <= length(Q)
            i = Q[k]
            for j in neighbors(J,i) 
                if !R[j] && σ[i] == σ[j] && rand() < p
                    R[j] = true
                    push!(Q, j)
                end
            end
            k += 1
        end
        # flip region Q
        σ[Q] .*= -1
        stats(t, σ, Q)
        # reset R and Q
        R[Q] .= false
        empty!(Q)
    end
end