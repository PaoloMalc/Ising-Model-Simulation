function magnetization(J, β, nsteps, n_therm)
    mags = Float64[]
    N = size(J, 1) #NUmero spin
    
    
    my_stats = (t, σ, Q) -> begin
        if t > n_therm 
            m = abs(sum(σ)) / N
            push!(mags, m)
        end
    end
    
   
    wolff(J, β, nsteps, stats=my_stats)
    
    
    return mean(mags)
end