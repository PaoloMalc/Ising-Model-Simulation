function wolff_spin_glass(J, β, nsteps; stats = (t, σ, Q)->nothing) 
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
            for idx in J.colptr[i]:(J.colptr[i+1]-1)
                j = J.rowval[idx]
                J_ij = J.nzval[idx]
                
                # La probabilità dipende dal valore assoluto del legame
                p_ij = 1 - exp(-2 * β * abs(J_ij))
                
                if !R[j] && (σ[i] * σ[j] == sign(J_ij)) && rand() < p_ij #modifico le condizioni di wolff per rimanere coerente con la DBC
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

#magnetizzazione media = 0 ? 
function energy(J, σ)
    E = 0.0
    N = length(σ)
    for i in 1:N
        for idx in J.colptr[i]:(J.colptr[i+1]-1)
            j = J.rowval[idx]
            J_ij = J.nzval[idx]
            E += - J_ij * σ[i] * σ[j]
        end
    end
    # Dividiamo per 2 perché il ciclo passa due volte su ogni legame (i->j e j->i)
    return E / 2.0 
end


function randomh_lattice(L, J_val)
    J_base = clattice(L, L) * J_val
    N = L * L
    
    I_1, J_1, v_1 = findnz(J_base)
    
    I_new = copy(I_1)
    J_new = copy(J_1)
    V_new = copy(v_1)
    
    idx_field = N + 1
    
    #collego ogni spin allo spin N+1
    for i in 1:N
        h_i = rand() * 2.0 - 1.0
        # Legame da i -> spin N+1
        push!(I_new, i)
        push!(J_new, idx_field)
        push!(V_new, h_i) #lascio h_i normalmente
        
        #inverso
        push!(I_new, idx_field)
        push!(J_new, i)
        push!(V_new, h_i) 
    end
    
    #Ricompongo la matrice
    return sparse(I_new, J_new, V_new, N+1, N+1)
end