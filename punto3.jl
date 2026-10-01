function magnetization_3d(J, β, nsteps, n_therm)
    mags = Float64[]
    N = size(J, 1) 
    
    my_stats = (t, σ, Q) -> begin
        if t > n_therm
            push!(mags, abs(sum(σ)) / N)
        end
    end
    
    wolff(J, β, nsteps, stats=my_stats)

    
    m_avg = mean(mags)
    
    #chi = β * N * Varianza(M)
    m_var = var(mags)
    chi = β * N * m_var 
    
    return m_avg, chi
end