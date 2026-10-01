using LinearAlgebra, SparseArrays

lattice(n) = spdiagm(-1 => fill(1,n-1), 1 => fill(1,n-1))
clattice(n) = spdiagm(-1 => fill(1,n-1), 1 => fill(1,n-1), -(n-1) => [1], n-1 => [1])
A⊕B = kron(A,I(size(B,1))) + kron(I(size(A,1)),B)
lattice(dims...) = lattice(dims[1])⊕lattice(dims[2:end]...)
clattice(dims...) = clattice(dims[1])⊕clattice(dims[2:end]...)



neighbors(J,i) = @view rowvals(J)[nzrange(J, i)]

function wolff(J, β, nsteps; stats = (t, σ, Q)->nothing)
    n = size(J, 1)
    σ = rand((-1,1), n)
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