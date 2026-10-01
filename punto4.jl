function wolff_hetero(J, β, nsteps; stats = (t, σ, Q)->nothing)
    n = size(J, 1)
    σ = rand((-1,1), n)
    R = fill(false, n)
    Q = Int[]

    for t in 1:nsteps
        push!(Q, rand(1:n))
        R[Q[1]] = true
        k = 1
        while k <= length(Q)
            i = Q[k]
            # leggo i vicini e l'intensità J_ij
            for idx in J.colptr[i]:(J.colptr[i+1]-1)
                j = J.rowval[idx]
                J_ij = J.nzval[idx]
                
                #p dipendente dal legame
                p_ij = 1 - exp(-2 * β * J_ij)
                
            
                if !R[j] && σ[i] == σ[j] && rand() < p_ij
                    R[j] = true
                    push!(Q, j)
                end
            end
            k += 1
        end
        σ[Q] .*= -1
        stats(t, σ, Q)
        R[Q] .= false
        empty!(Q)
    end
end

function heteroising(J, β, nsteps, n_therm)
    mags = Float64[]
    N = size(J, 1) 
    my_stats = (t, σ, Q) -> begin
        if t > n_therm
            push!(mags, abs(sum(σ)) / N)
        end
    end
    wolff_hetero(J, β, nsteps, stats=my_stats)
    m_avg = mean(mags)
    m_var = var(mags)
    chi = β * N * m_var
    
    return m_avg, chi
end
