function h_lattice(L, J_val, h_val)
    J_base = clattice(L, L) * J_val
    N = L * L
    
    I_1, J_1, V_1 = findnz(J_base)
    
    I_new = copy(I_1)  #indice riga 
    J_new = copy(J_1)  #indice colonna 
    V_new = copy(V_1)  #valore in quella posizione
    
    #es. I[5] J[5] V[5] riga 5 colonna 5, collegamento di valore V[5]
    idx_field = N + 1
    
    #collego ogni spin allo spin N+1
    for i in 1:N
        # Legame da i -> spin N+1
        push!(I_new, i)
        push!(J_new, idx_field)
        push!(V_new, abs(h_val)) # Wolff_hetero funziona per J > 0, posso gestire il segno alla fine
        
        #inverso
        push!(I_new, idx_field)
        push!(J_new, i)
        push!(V_new, abs(h_val)) 
    end
    
    #Ricompongo la matrice
    return sparse(I_new, J_new, V_new, N+1, N+1)
end


#in alternativa posso anche creare una J iniziale con una riga e una colonna in più con i valori del campo (es. 0.1)
function h_lattice_blocchi(L, J_val, h_val)
    J_base = clattice(L, L) * J_val
    N = L * L
    
    h_col = sparse(fill(abs(h_val), N))
    
    J_new = [J_base h_col; h_col'  0.0  ]
    return J_new
end


function h_ising(J_h, β, h_val, nsteps, n_therm)
    mags = Float64[]
    N_tot = size(J_h, 1) 
    N_spin = N_tot - 1 
    my_stats = (t, σ, Q) -> begin
        if t > n_therm
            # Non uso il valore assoluto
            M_field = sum(σ[1:N_spin]) * σ[end] * sign(h_val) #quando moltiplica per sigma end moltiplica per il campo e eventualmente ribalta. Procedimento validato dalla proprietà (a)
            push!(mags, M_field / N_spin)
        end
    end
    
    wolff_hetero(J_h, β, nsteps, stats=my_stats)
    
    return mean(mags)
end


function track_mag_zero_field(J_base, β, nsteps)
    mags_t = Float64[]
    N = size(J_base, 1)
    
    #non uso abs()
    my_stats = (t, σ, Q) -> push!(mags_t, sum(σ) / N) 
    
    wolff_hetero(J_base, β, nsteps, stats=my_stats)
    return mags_t
end


function track_mag_with_field(J_h, β, h_val, nsteps)
    mags_t = Float64[]
    N_spin = size(J_h, 1) - 1
    
    my_stats = (t, σ, Q) -> begin
        M_field = sum(σ[1:N_spin]) * σ[end] * sign(h_val)
        push!(mags_t, M_field / N_spin)
    end
    
    wolff_hetero(J_h, β, nsteps, stats=my_stats)
    return mags_t
end
