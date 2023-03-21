
function direct_product(symtext, irrep_1, irrep_2)
    c1 = symtext.ctab[irrep_1]
    c2 = symtext.ctab[irrep_2]
    co = symtext.ctab.class_orders
    v = c1.*c2.*co
    out = zeros(Int64, length(symtext.ctab.irreps))
    for (idx,irrep) in enumerate(symtext.ctab.irreps)
        out[idx] = round(Int64,sum(v.*symtext.ctab[irrep])/symtext.order)
    end
    return out
end

function Dtriple(intcol, α, β, γ, i, j, k, l, m, n)
    suum = 0.0
    irrm = intcol.irrm
    for S = 1:intcol.g
        #c = irrm[α][S][i,j] * irrm[β][S][k,l] * irrm[γ][S][m,n]
        #println(c)
        #suum += c
        suum += irrm[α][S][i,j] * irrm[β][S][k,l] * irrm[γ][S][m,n]
    end
    return suum * (intcol.irr_dims[γ]/intcol.g)
end

function build_ClebschGordan(intcol)
    # C[irrep_A, pfxn_A, g_A, irrep_i, pfxn_i, irrep_j, pfxn_j]
    #println(direct_product(intcol.symtext, "T2", "T2"))
    # Find a non-zero CG coefficient for each unique irrep triple such that i=j, k=l, and m=n
    max_irr_dim = maximum([i for i in values(intcol.irr_dims)])
    number_irreps = length(intcol.irr_dims)
    CG = zeros(Float64, (number_irreps, number_irreps, number_irreps, max_irr_dim, max_irr_dim, max_irr_dim, 2))
    #CG = zeros(Array{Vector{Float64},3}, intcol.g, intcol.g, intcol.g)
    for (Midx,Mirr) in enumerate(intcol.symtext.ctab.irreps)
        for (Nidx,Nirr) in enumerate(intcol.symtext.ctab.irreps)
            dir_prod = direct_product(intcol.symtext, Mirr, Nirr)
            for (Aidx,Airr) in enumerate(intcol.symtext.ctab.irreps)
                #CG[Midx,Nidx,Aidx] = zeros(Vector{Float64}, intcol.irr_dims[Mirr], intcol.irr_dims[Nirr], intcol.irr_dims[Airr])
                # Ignore where one of M or N are the totally symmetric irrep because the coefficient is one
                # Similarly ignore A if A ∉ M ⊗ N, and only look at unique pairs of M and N since M ⊗ N = N ⊗ M
                if dir_prod[Aidx] == 0 || Midx == 1 || Nidx == 1 || Midx > Nidx
                    continue
                end
                for i = 1:intcol.irr_dims[Mirr], k = 1:intcol.irr_dims[Nirr], m = 1:intcol.irr_dims[Airr]
                    c2 = Dtriple(intcol, Mirr, Nirr, Airr, i, i, k, k, m, m)
                    if abs(c2) < 1e-10
                        continue
                    end
                    if dir_prod[Aidx] == 1
                        C_ref = sqrt(c2)
                        CG[Midx,Nidx,Aidx,i,k,m,1] = sqrt(c2)
                        b = ""
                    elseif dir_prod[Aidx] == 2
                        C_ref = sqrt(c2)/2
                        CG[Midx,Nidx,Aidx,i,k,m,1] = c
                        CG[Midx,Nidx,Aidx,i,k,m,2] = c
                        b = "Look here: "
                    else
                        throw(ErrorException("Invalid result for direct product of two irreps: $(dir_prod[Aidx])"))
                    end
                    println(b*"$Mirr ⊗ $Nirr = $Airr: $C_ref")
                    for j = 1:intcol.irr_dims[Mirr], l = 1:intcol.irr_dims[Nirr], n = 1:intcol.irr_dims[Airr]
                        c1c2 = Dtriple(intcol, Mirr, Nirr, Airr, i, j, k, l, m, n)
                        if dir_prod[Aidx] == 1
                            CG[Midx,Nidx,Aidx,j,l,n,1] = c1c2 / C_ref
                        else
                            CG[Midx,Nidx,Aidx,j,l,n,1] = c1c2/(2*C_ref)
                            CG[Midx,Nidx,Aidx,j,l,n,2] = c1c2/(2*C_ref)
                        end
                    end
                    break
                end
                CG[Nidx,Midx,Aidx,:,:,:,:] = CG[Midx,Nidx,Aidx,:,:,:,:]
            end
        end
    end
    return CG
end

function transform_integral(CG, integral, M, N, P, Q, m, n, p, q, mp, np, pp, qp)
    return 0
end
