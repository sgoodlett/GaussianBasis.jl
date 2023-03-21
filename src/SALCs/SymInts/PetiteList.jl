function poopoo(salc1, salc2)
    chk1 = salc1.irrep == salc2.irrep
    chk2 = salc1.sh == salc2.sh
    chk3 = salc1.atom == salc2.atom
    chk4 = salc1.bfxn == salc2.bfxn
    chk5 = salc1.r == salc2.r
    if chk1 && chk2 && chk3 && chk4 && chk5
        return true
    end
    return false
end

function find_unique_salcs(salcs)
    # Group partner functions together
    unique_salcs = [(1,salcs[1])]
    out = [1]
    for (sidx, salc) in enumerate(salcs[2:end])
        chork = true
        for (usidx,usalc) in enumerate(unique_salcs)
            if poopoo(salc, usalc[2])
                push!(out, usalc[1])
                chork = false
                break
            end
        end
        if chork
            push!(unique_salcs, (sidx+1,salc))
            push!(out, sidx+1)
        end
    end
    return unique_salcs, out
end

function get_petite_idxs(salcs)
    nsalcs = length(salcs)
    unique_salcs = find_unique_salcs(salcs)[2]
    petite_list = []
    for i = 1:nsalcs
        for j = i:nsalcs
            for k = 1:nsalcs
                for l = k:nsalcs
                    if GaussianBasis.index2(i,j) < GaussianBasis.index2(k,l)
                        break
                    end
                    if i == unique_salcs[i] && j == unique_salcs[j] && k == unique_salcs[k] && l == unique_salcs[l]
                        push!(petite_list, [i,j,k,l])
                    end
                end
            end
        end
    end
    return petite_list
end
