
function get_atom_subgroup(atomidx, symtext)
    subgroup = Int64[]
    for (i,v) in enumerate(symtext.atom_map[atomidx,:])
        if v == atomidx
            push!(subgroup, i)
        end
    end
    return subgroup
end

function dummy_test(fn, a1, a2)
    to = TimerOutput()
    mol, symtext = Molecules.Symmetry.CharacterTables.symtext_from_file(fn)
    U = get_atom_subgroup(a1, symtext)
    V = get_atom_subgroup(a2, symtext)
    dcs = get_double_cosets_legacy(to, U, V, symtext)
    println(dcs)
    @timeit to "New" dcs = get_double_cosets2(U, V, symtext)
    println(dcs)
    show(to)
end

function get_double_cosets(U, V, symtext)
    found = Int64[]
    double_cosets = Vector{Int64}[]
    for g = 1:length(symtext.symels)
        if symtext.mult_table[U[1],symtext.mult_table[g, V[1]]] in found
            continue
        end
        double_coset = Int64[]
        for u in U, v in V
            ugv = symtext.mult_table[u,symtext.mult_table[g,v]]
            push!(found, ugv)
            push!(double_coset, ugv)
        end
        push!(double_cosets, double_coset)
    end
    return double_cosets
end

function get_double_cosets_legacy(to, U, V, symtext)
    @timeit to "All double cosets" begin
    double_cosets = []
    for g = 1:length(symtext.symels)
        double_coset = []
        for u in U, v in V
            push!(double_coset, symtext.mult_table[u,symtext.mult_table[g,v]])
        end
        push!(double_cosets, double_coset)
    end
    end #time
    @timeit to "Unique double cosets" begin 
    unique_double_cosets = [sort(double_cosets[1])]
    for dcoset = 2:length(double_cosets)
        chk = true
        for udcoset in unique_double_cosets
            srtd_dcoset = sort(double_cosets[dcoset])
            if srtd_dcoset == udcoset
                chk = false
                break
            end
        end
        if chk
            push!(unique_double_cosets, sort(double_cosets[dcoset]))
        end
    end
    end # time
    return unique_double_cosets
end

function get_Rs(to, U, V, symtext)
    Rs = []
    λs = []
    #U = get_atom_subgroup(atom1idx, symtext)
    #V = get_atom_subgroup(atom2idx, symtext)
    @timeit to "Double Coset" double_cosets = get_double_cosets(U, V, symtext)
    for i in double_cosets
        push!(Rs, i[1])
        λ = count(==(i[1]), i)
        push!(λs, λ)
    end
    return Rs, λs
end

function special_Rs(U, symtext)
    h = length(symtext.symels)
    sats = []
    for g = 1:h
        ginv = 0
        u_g_u = []
        u_ginv_u = []
        for g2 = 1:h
            if symtext.mult_table[g,g2] == 1
                ginv = g2
                break
            end
        end
        for u in U
            push!(u_g_u, symtext.mult_table[u,symtext.mult_table[g,u]])
            push!(u_ginv_u, symtext.mult_table[u,symtext.mult_table[ginv,u]])
        end
        push!(sats, union(u_g_u,u_ginv_u))
    end
    println(sats)
    unique_sats = [sort(sats[1])]
    for i = 2:length(sats)
        chk = true
        srt_sat = sort(sats[i])
        for usat in unique_sats
            if srt_sat == usat
                chk = false
                break
            end
        end
        if chk
            push!(unique_sats, srt_sat)
        end
    end
    println(unique_sats)

    return 0
end

function convert()

end

struct IntGarbage
    symtext::Molecules.Symmetry.CharacterTables.SymText
    irrm::Dict{String, Vector{AbstractArray{Float64}}}
    fxnmap::Vector{Vector{Matrix{Float64}}}
    stat_subgrps::Vector{Vector{Int64}}
    g::Int64
    irr_dims::Dict{String, Int64}
end

function Λfxn(garbo::IntGarbage, salc, ab, p, R::Int64)
    out = 0.0
    for U in garbo.stat_subgrps[salc.atom]
        RU = garbo.symtext.mult_table[R, U]
        out += garbo.fxnmap[RU][salc.sh][ab, salc.ml]*garbo.irrm[salc.irrep][RU][p,salc.r]
    end
    return out * garbo.irr_dims[salc.irrep] / garbo.g
end

# Same function but R=1
function Λfxn(garbo::IntGarbage, salc, ab, p, printstuff=false)
    out = 0.0
    for U in garbo.stat_subgrps[salc.atom]
        if printstuff
            println("C[$U, $(salc.sh), $ab, $(salc.ml)] = $(garbo.fxnmap[U][salc.sh][ab, salc.ml])")
        end
        out += garbo.fxnmap[U][salc.sh][ab, salc.ml]*garbo.irrm[salc.irrep][U][p,salc.r]
    end
    return out * garbo.irr_dims[salc.irrep] / garbo.g
end

function Γfxn(garbo::IntGarbage, salc1, salc2, ab, bb, R, λr, printstuff::Bool)
    out = 0.0
    for p = 1:garbo.irr_dims[salc1.irrep]
        porque1 = Λfxn(garbo, salc1, ab, p, printstuff)
        porque2 = Λfxn(garbo, salc2, bb, p, R)
        if printstuff
            println("p: $p, Λ: $porque1, Λ(R): $porque2")
        end
        out += Λfxn(garbo, salc1, ab, p) * Λfxn(garbo, salc2, bb, p, R)
    end
    return out * garbo.g/(garbo.irr_dims[salc1.irrep] * λr)
end

function get_garbage(mol, symtext, bset)
    irrm = eval(Meta.parse("Molecules.Symmetry.CharacterTables.irrm_"*symtext.pg))
    subgrps = Vector{Int64}[]
    for A in 1:length(mol)
        push!(subgrps, get_atom_subgroup(A, symtext))
    end
    aotoso, ignore_this, stalcs, doodad = GaussianBasis.SALCs.ProjectionOp(mol, bset, symtext)[2:5]
    return IntGarbage(symtext, irrm, doodad, subgrps, symtext.order, symtext.ctab.irrep_dims), stalcs, aotoso
end

function get_Stalcs(mol, bset)
    stalcs, doodad = GaussianBasis.SALCs.ProjectionOp(mol, bset)[4:5]
    return stalcs, doodad
end

function build_thang(bset::BasisSet)
    thang = []
    sh_ctr = 0
    bfxn_ctr = 0
    for atom = 1:length(bset.atoms)
        thang2 = []
        for sh_idx = 1:bset.shells_per_atom[atom]
            l = bset[sh_idx+sh_ctr].l
            for l_idx = 1:2l+1
                push!(thang2, (sh_idx+sh_ctr, collect(1+bfxn_ctr:2l+1+bfxn_ctr)))
            end
            bfxn_ctr += 2l+1
        end
        sh_ctr += bset.shells_per_atom[atom]
        push!(thang, thang2)
    end
    return thang
end

function check_salcs(salcs, aotoso)
    results = []
    for (i, salc) in enumerate(salcs)
        m = findmax(broadcast(abs, salc.coeffs-aotoso[:,i]))
        push!(results,m[1])
    end
    return results
end

function direct_product(irr_list, symtext)
    r = ones(Float64, length(symtext.ctab.classes))
    for (i, irr) = enumerate(irr_list)
        r .*= symtext.ctab[irr]
    end
    return r
end

function tran_as_A1(irr_list, symtext)
    if sum(direct_product(irr_list,symtext).*symtext.ctab.class_orders) > 1E-6
        return true
    else
        return false
    end
end

function sintegrals(molfn::String, basis::String, int_type::Function)
    mol = Molecules.parse_file(molfn)
    mol, symtext = Molecules.Symmetry.CharacterTables.symtext_from_mol(mol)
    return sintegrals(mol, basis, int_type, symtext)
end

#function sintegrals(mol::Vector{Molecules.Atom}, basis::String, int_type::Function)
#    bset = BasisSet(basis, mol)
#    return sintegrals(mol, bset, int_type)
#end
#
#function sintegrals(mol::Vector{Molecules.Atom}, basis::Vector{SphericalShell{Molecules.Atom{T,T},1,T}}, int_type::Function) where {T}
#    bset = BasisSet("NewBasis", mol, basis)
#    return sintegrals(mol, bset, int_type)
#end

function sintegrals(mol::Vector{Molecules.Atom}, basis, int_type, symtext)
    to = TimerOutput()
    bset = BasisSet(basis, mol)
    display(mol)
    @timeit to "Get stuff" begin
    abars = build_thang(bset)
    int_garb, stalcs, aotoso = get_garbage(mol, symtext, bset) end
    @timeit to "Ints hard" begin
    Sboy = int_type(bset)
    sym_S = transpose(aotoso) * Sboy * aotoso end
    #for i = 1:length(stalcs)
    #    println("Irrep.: $(stalcs[i].irrep), Atom: $(stalcs[i].atom), Bsfxn: $(stalcs[i].bfxn), sh, ml: $(stalcs[i].sh), $(stalcs[i].ml), i,r: $(stalcs[i].i), $(stalcs[i].r), Factor:$(stalcs[i].γ)")
    #end
    println("Max dif. in salcs: ", findmax(broadcast(abs, check_salcs(stalcs, aotoso)))[1])
    slength = length(stalcs)
    S = zeros(Float64, slength, slength)
    @timeit to "My Ints" begin
    for s1 = 1:slength
        for s2 = s1:slength
            salc1 = stalcs[s1]
            salc2 = stalcs[s2]
            printstuff = false
            if salc1.irrep == salc2.irrep && salc1.i == salc2.i
                U = int_garb.stat_subgrps[salc1.atom]
                if salc1.atom == salc2.atom && salc1.bfxn == salc2.bfxn && false # I don't want to deal with this yet, so we can calculate redundant integrals and I won't feel bad
                    # Do the other thing
                    Rs = special_Rs(U, symtext)
                else
                    V = int_garb.stat_subgrps[salc2.atom]
                    Rs, λs = get_Rs(U, V, int_garb.symtext)
                    for (Ridx,R) in enumerate(Rs)
                        Rb = int_garb.symtext.atom_map[salc2.atom,R]
                        for (abidx,ab) in enumerate(abars[salc1.atom][salc1.bfxn][2])
                            for (bbidx,bb) in enumerate(abars[salc2.atom][salc2.bfxn][2])
                                @timeit to "Γ" Γ = Γfxn(int_garb, salc1, salc2, abidx, bbidx, R, λs[Ridx], printstuff)
                                @timeit to "AO ints" out = int_type(bset, abars[salc1.atom][salc1.bfxn][1], abars[Rb][salc2.bfxn][1])[abidx,bbidx] # shell_a, shell_b, ml_ab, ml_bb
                                S[s1, s2] += Γ*out*salc1.γ*salc2.γ
                                if s1 != s2
                                    S[s2, s1] += Γ*out*salc1.γ*salc2.γ
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    end #time
    #display(sym_S)
    println(findmax(broadcast(abs, S-sym_S)))
    show(to)
    return 0
end

function calcH(H, to, int_garb, U, V, W, X, R, S, T, λTs, Tidx, salc1, salc2, salc3, salc4, abidx, bbidx, cbidx, dbidx)
    #Hs = [0.0 for i = 1:Threads.nthreads()]
    #@sync begin
    for v in V
        #Ht = Hs[Threads.threadid()]
        rv = int_garb.symtext.mult_table[R,v]
        for w in W
            tw = int_garb.symtext.mult_table[T,w]
            for x in X
                tsx = int_garb.symtext.mult_table[T,int_garb.symtext.mult_table[S,x]]
                for u in U
                    @timeit to "Sum over G" begin
                    #Threads.@spawn begin
                    #Ht = Hs[Threads.threadid()]
                    suuum = 0.0
                    for g = 1:int_garb.g
                        gu = int_garb.symtext.mult_table[g,u]
                        grv = int_garb.symtext.mult_table[g,rv]
                        gtw = int_garb.symtext.mult_table[g,tw]
                        gtsx = int_garb.symtext.mult_table[g,tsx]
                        #suuum += int_garb.irrm[salc1.irrep][int_garb.symtext.mult_table[g,u]][salc1.i,salc1.r]*int_garb.irrm[salc2.irrep][int_garb.symtext.mult_table[g,rv]][salc2.i,salc2.r]*int_garb.irrm[salc3.irrep][int_garb.symtext.mult_table[g,tw]][salc3.i,salc3.r]*int_garb.irrm[salc4.irrep][int_garb.symtext.mult_table[g,tsx]][salc4.i,salc4.r]
                        suuum += int_garb.irrm[salc1.irrep][gu][salc1.i,salc1.r]*int_garb.irrm[salc2.irrep][grv][salc2.i,salc2.r]*int_garb.irrm[salc3.irrep][gtw][salc3.i,salc3.r]*int_garb.irrm[salc4.irrep][gtsx][salc4.i,salc4.r]
                    end
                    #Hs[Threads.threadid()] += sum(Dgu.*Dgrv.*Dgtw.*Dgtsx) * int_garb.fxnmap[u][salc1.sh][abidx, salc1.ml]*int_garb.fxnmap[rv][salc2.sh][bbidx, salc2.ml]*int_garb.fxnmap[tw][salc3.sh][cbidx, salc3.ml]*int_garb.fxnmap[tsx][salc4.sh][dbidx, salc4.ml]
                    H += suuum * int_garb.fxnmap[u][salc1.sh][abidx, salc1.ml]*int_garb.fxnmap[rv][salc2.sh][bbidx, salc2.ml]*int_garb.fxnmap[tw][salc3.sh][cbidx, salc3.ml]*int_garb.fxnmap[tsx][salc4.sh][dbidx, salc4.ml]
                    end
                end
            end
        end
    end
    #end
    #H = sum(Hs)
    H *= (int_garb.irr_dims[salc1.irrep]*int_garb.irr_dims[salc2.irrep]*int_garb.irr_dims[salc3.irrep]*int_garb.irr_dims[salc4.irrep])/(λTs[Tidx]*(int_garb.g^4))
    return H
end

function calcHp(H, to, Da, Db, Dc, Dd, sDa, sDb, sDc, sDd, Ca, Cb, Cc, Cd)
        for u = 1:sDa
            for v = 1:sDb
                for w = 1:sDc
                    @simd for x = 1:sDd
                        @inbounds H += sum(Da[:,u] .* Db[:,v] .* Dc[:,w] .* Dd[:,x]) * Ca[u] * Cb[v] * Cc[w] * Cd[x]
                    end
                end
            end
        end
    return H
end

function calcHpp(H, to, Da, Db, Dc, Dd, Ca, Cb, Cc, Cd)
    @ein H[] := Da[g,u]*Db[g,v]*Dc[g,w]*Dd[g,x]*Ca[u]*Cb[v]*Cc[w]*Cd[x]
    return H[1]
end



function genCmat(fxnmap, U, abidx, sh, ml)
    return Float64[fxnmap[u][sh][abidx, ml] for u in U]
end

function genD2mat(mtab, irrm, irrep, i, r, G, U)
    return Float64[irrm[irrep][mtab[g,u]][i,r] for g = 1:G, u in U]
end

function genD3mat(mtab, irrm, irrep, i, r, G, U, V)
    return Float64[irrm[irrep][mtab[g,mtab[u,v]]][i,r] for g = 1:G, u in U, v in V]
end

function genC2matp(fxnmap, U, nf, sh, ml)
    return Float64[fxnmap[u][sh][abidx, ml] for u in U, abidx =1:nf]
end

function genC3matp(mtab, fxnmap, U, V, nf, sh, ml)
    return Float64[fxnmap[mtab[u,v]][sh][abidx, ml] for u in U, v in V, abidx =1:nf]
end

function stwintegrals(mol::String, bset)
    mol, symtext = Molecules.Symmetry.CharacterTables.symtext_from_file(mol)
    return stwintegrals(mol, bset, symtext)
end

function stwintegrals(mol::Vector{Molecules.Atom}, basis, symtext)
    to = TimerOutput()
    bset = BasisSet(basis, mol)
    @timeit to "Abars and IntGarb" begin
    abars = build_thang(bset)
    int_garb, salcs, aotoso = get_garbage(mol, symtext, bset)
    println("Number of basis functions $(length(salcs))")
    end
    @timeit to "Ints the hard way" begin
    ERI_AO = ERI_2e4c(bset)
    SERI = zeros(bset.nbas, bset.nbas, bset.nbas, bset.nbas)
    @tensoropt SERI[i,j,k,l] =  aotoso[μ, i]*aotoso[ν, j]*ERI_AO[μ, ν, ρ, σ]*aotoso[ρ, k]*aotoso[σ, l]
    end
    #for i = 1:length(salcs)
    #    println("Irrep.: $(salcs[i].irrep), Atom: $(salcs[i].atom), Bsfxn: $(salcs[i].bfxn), sh, ml: $(salcs[i].sh), $(salcs[i].ml), i,r: $(salcs[i].i), $(salcs[i].r)")
    #end
    println(check_salcs(salcs, aotoso))
    porque = zeros(Float64, size(aotoso))
    for i = 1:length(salcs)
        porque[:,i] = aotoso[:,i]-salcs[i].coeffs
    end
    println(findmax(broadcast(abs, porque)))
    slength = length(salcs)
    eri = zeros(Float64, slength, slength, slength, slength)
    #eris = [zeros(Float64, slength, slength, slength, slength) for i in 1:Threads.nthreads()]
    @timeit to "My Ints" begin
    #@sync begin
    for s1 = 1:slength
        salc1 = salcs[s1]
        U = int_garb.stat_subgrps[salc1.atom]
        Da = genDmat(int_garb.symtext.mult_table, int_garb.irrm, salc1.irrep, salc1.i, salc1.r, int_garb.g, U)
        #sDa = size(Da,2)
        s2runs = 0
        s3runs = 0
        s4runs = 0
        for s2 = s1:slength
            s2runs += 1
            salc2 = salcs[s2]
            V = int_garb.stat_subgrps[salc2.atom]
            @timeit to "get_Rs" Rs, λRs = get_Rs(to, U, V, int_garb.symtext)
            for R in Rs
                RV = Int64[int_garb.symtext.mult_table[R,v] for v in V]
                Db = genDmat(int_garb.symtext.mult_table, int_garb.irrm, salc2.irrep, salc2.i, salc2.r, int_garb.g, RV)
                #sDb = size(Db,2)
                Rb = int_garb.symtext.atom_map[salc2.atom,R]
                M = intersect(int_garb.stat_subgrps[salc1.atom], int_garb.stat_subgrps[int_garb.symtext.atom_map[salc2.atom, R]])
                for s3 = 1:slength
                    s3runs += 1
                    salc3 = salcs[s3]
                    W = int_garb.stat_subgrps[salc3.atom]
                    for s4 = s3:slength
                        salc4 = salcs[s4]
                        i12 = GaussianBasis.index2(s1, s2)
                        i34 = GaussianBasis.index2(s3, s4)
                        if i12 < i34
                            continue
                        end
                        if !tran_as_A1([salc1.irrep, salc2.irrep, salc3.irrep, salc4.irrep], symtext)
                            continue
                        end
                        s4runs += 1
                        X = int_garb.stat_subgrps[salc4.atom]
                        @timeit to "get_Ss" Ss, λSs = get_Rs(to, W, X, int_garb.symtext)
                        #Threads.@spawn begin
                            to2 = TimerOutput()
                            for S in Ss
                                N = intersect(int_garb.stat_subgrps[salc3.atom], int_garb.stat_subgrps[int_garb.symtext.atom_map[salc4.atom, S]])
                                @timeit to2 "get_Ts" Ts, λTs = get_Rs(to2, M, N, int_garb.symtext)
                                for (Tidx,T) in enumerate(Ts)
                                    TW = Int64[int_garb.symtext.mult_table[T,w] for w in W]                                
                                    TSX = Int64[int_garb.symtext.mult_table[T,int_garb.symtext.mult_table[S,x]] for x in X]                                
                                    @timeit to2 "Build Dc and Dd" begin
                                        Dc = genDmat(int_garb.symtext.mult_table, int_garb.irrm, salc3.irrep, salc3.i, salc3.r, int_garb.g, TW)
                                        Dd = genDmat(int_garb.symtext.mult_table, int_garb.irrm, salc4.irrep, salc4.i, salc4.r, int_garb.g, TSX)
                                        #sDc = size(Dc,2)
                                        #sDd = size(Dd,2)
                                    end
                                    #for (abidx,ab) in enumerate(abars[salc1.atom][salc1.bfxn][2])
                                    #Rb = int_garb.symtext.atom_map[salc2.atom,R]
                                    Tc = int_garb.symtext.atom_map[salc3.atom,T]
                                    TSd = int_garb.symtext.atom_map[salc4.atom,int_garb.symtext.mult_table[T,S]]
                                    @timeit to2 "AO int calc" out = ERI_2e4c(bset, abars[salc1.atom][salc1.bfxn][1], abars[Rb][salc2.bfxn][1], abars[Tc][salc3.bfxn][1], abars[TSd][salc4.bfxn][1])
                                    @timeit to2 "Gen. C mats" begin
                                        Ca = genCmatp(int_garb.fxnmap, U, 2*(salc1.sh-1)+1, salc1.sh, salc1.ml)
                                        Cb = genCmatp(int_garb.fxnmap, RV, 2*(salc2.sh-1)+1, salc2.sh, salc2.ml)
                                        Cc = genCmatp(int_garb.fxnmap, TW, 2*(salc3.sh-1)+1, salc3.sh, salc3.ml)
                                        Cd = genCmatp(int_garb.fxnmap, TSX, 2*(salc4.sh-1)+1, salc4.sh, salc4.ml)
                                    end
                                    @timeit to2 "Calc H and Contract" begin
                                        @timeit to2 "Contract" @tensoropt trout[g1,g2,g3,g4] := Da[g1,u]*Db[g2,v]*Dc[g3,w]*Dd[g4,x]*Ca[u,ab]*Cb[v,bb]*Cc[w,cb]*Cd[x,db]*out[ab,bb,cb,db]
                                        summa = 0.0
                                        for g = 1:int_garb.g
                                            summa += trout[g,g,g,g]
                                        end
                                        #@ein summa[] := Da[g,u]*Db[g,v]*Dc[g,w]*Dd[g,x]*Ca[u,ab]*Cb[v,bb]*Cc[w,cb]*Cd[x,db]*out[ab,bb,cb,db]
                                        #summa = sum(summa)
                                    end
                                    #summa = summa[1]
                                    #@timeit to2 "Idx loops" begin
                                    #summa = 0.0
                                    #for abidx = 1:2*(salc1.sh-1)+1
                                    #    Ca = genCmat(int_garb.fxnmap, U, abidx, salc1.sh, salc1.ml)
                                    #    for bbidx = 1:2*(salc2.sh-1)+1
                                    #        Cb = genCmat(int_garb.fxnmap, RV, bbidx, salc2.sh, salc2.ml)
                                    #        for cbidx = 1:2*(salc3.sh-1)+1
                                    #            Cc = genCmat(int_garb.fxnmap, TW, cbidx, salc3.sh, salc3.ml)
                                    #            for dbidx = 1:2*(salc4.sh-1)+1
                                    #                @timeit to2 "Build Cd" Cd = genCmat(int_garb.fxnmap, TSX, dbidx, salc4.sh, salc4.ml)
                                    #                @timeit to2 "Calc H" begin
                                    #                    H = 0.0
                                    #                    H = calcHpp(H, to, Da, Db, Dc, Dd, Ca, Cb, Cc, Cd)
                                    #                    summa += H*out[abidx,bbidx,cbidx,dbidx]#*salc1.γ*salc2.γ*salc3.γ*salc4.γ
                                    #                end
                                    #            end
                                    #        end
                                    #    end
                                    #end
                                    #end
                                    @timeit to2 "Int assignment" begin
                                        summa *= (int_garb.irr_dims[salc1.irrep]*int_garb.irr_dims[salc2.irrep]*int_garb.irr_dims[salc3.irrep]*int_garb.irr_dims[salc4.irrep])/(λTs[Tidx]*(int_garb.g^4))
                                        troot = summa*salc1.γ*salc2.γ*salc3.γ*salc4.γ
                                        #eri = eris[Threads.threadid()]
                                        eri[s1,s2,s3,s4] += troot
                                        if s1 != s2
                                            eri[s2,s1,s3,s4] += troot
                                            if s3 != s4
                                                eri[s1,s2,s4,s3] += troot
                                                eri[s2,s1,s4,s3] += troot
                                                if i12 != i34
                                                    eri[s3,s4,s1,s2] += troot
                                                    eri[s3,s4,s2,s1] += troot
                                                    eri[s4,s3,s1,s2] += troot
                                                    eri[s4,s3,s2,s1] += troot
                                                end
                                            elseif i12 != i34
                                                eri[s3,s4,s1,s2] += troot
                                                eri[s3,s4,s2,s1] += troot
                                            end
                                        elseif s3 != s4
                                            eri[s1,s2,s4,s3] += troot
                                            if i12 != i34
                                                eri[s3,s4,s1,s2] += troot    
                                                eri[s4,s3,s1,s2] += troot    
                                            end
                                        elseif i12 != i34
                                            eri[s3,s4,s1,s2] += troot
                                        end
                                    end

                                end
                            end
                            merge!(to, to2, tree_point = ["My Ints"])
                        #end
                        #merge!(to, to2)
                    end
                end
            end
        end
        println("*Salc 1: $s1, S2s: $s2runs, S3s: $s3runs, S4s: $s4runs")

    end
    end
    #end
    #display()
    #eri = sum(eris)
    println(findmax(broadcast(abs, eri-SERI)))
    show(to)
    return 0
    #return eri
end

function stwintegrals_fromD(mol::String, basis::String)
    mol, symtext = Molecules.Symmetry.CharacterTables.symtext_from_file(mol)
    return stwintegrals_fromD(mol, basis, symtext)
end

function bigD_setup(to, mol, basis, symtext)
    bset = BasisSet(basis, mol)
    @timeit to "Abars and IntGarb" begin
    abars = build_thang(bset)
    int_garb, salcs, aotoso = get_garbage(mol, symtext, bset)
    println("Number of basis functions $(length(salcs))")
    end
    println(check_salcs(salcs, aotoso))
    porque = zeros(Float64, size(aotoso))
    for i = 1:length(salcs)
        porque[:,i] = aotoso[:,i]-salcs[i].coeffs
    end
    println(findmax(broadcast(abs, porque)))

    return bset, int_garb, salcs, abars, aotoso
end

function stwintegrals_fromD(mol, basis, symtext)
    to = TimerOutput()
    bset, int_garb, salcs, abars, aotoso = bigD_setup(to, mol, basis, symtext)
    @timeit to "Ints the hard way" begin
        ERI_AO = ERI_2e4c(bset)
        SERI = zeros(bset.nbas, bset.nbas, bset.nbas, bset.nbas)
        @tensoropt SERI[i,j,k,l] =  aotoso[μ, i]*aotoso[ν, j]*ERI_AO[μ, ν, ρ, σ]*aotoso[ρ, k]*aotoso[σ, l]
    end
    slength = length(salcs)
    eri = zeros(Float64, slength, slength, slength, slength)
    load_bar = "Calculating Integrals\n"*"-"^slength
    println(load_bar)
    @timeit to "My Ints" begin
    for s1 = 1:slength
        print("*")
        for s2 = s1:slength
            for s3 = 1:slength
                for s4 = s3:slength
                    @timeit to "Filter" begin
                    i12 = GaussianBasis.index2(s1, s2)
                    i34 = GaussianBasis.index2(s3, s4)
                    if i12 < i34
                        continue
                    end
                    if !tran_as_A1([salcs[s1].irrep, salcs[s2].irrep, salcs[s3].irrep, salcs[s4].irrep], symtext)
                        continue
                    end
                    end
                    @timeit to "Direct Symint" troot = stwintegrals_bigD(to, bset, salcs, int_garb, abars, s1, s2, s3, s4)
                    @timeit to "Fill ERI Array" begin
                    eri[s1,s2,s3,s4] += troot
                    if s1 != s2
                        eri[s2,s1,s3,s4] += troot
                        if s3 != s4
                            eri[s1,s2,s4,s3] += troot
                            eri[s2,s1,s4,s3] += troot
                            if i12 != i34
                                eri[s3,s4,s1,s2] += troot
                                eri[s3,s4,s2,s1] += troot
                                eri[s4,s3,s1,s2] += troot
                                eri[s4,s3,s2,s1] += troot
                            end
                        elseif i12 != i34
                            eri[s3,s4,s1,s2] += troot
                            eri[s3,s4,s2,s1] += troot
                        end
                    elseif s3 != s4
                        eri[s1,s2,s4,s3] += troot
                        if i12 != i34
                            eri[s3,s4,s1,s2] += troot    
                            eri[s4,s3,s1,s2] += troot    
                        end
                    elseif i12 != i34
                        eri[s3,s4,s1,s2] += troot
                    end
                    end
                end
            end
        end
    end
    end
    print("\n")
    println(findmax(broadcast(abs, eri-SERI)))
    show(to)
end

function stwintegrals_bigD(to, bset, salcs, int_garb, abars, s1, s2, s3, s4)
    salc1 = salcs[s1]
    salc2 = salcs[s2]
    salc3 = salcs[s3]
    salc4 = salcs[s4]
    U = int_garb.stat_subgrps[salc1.atom]
    Da = genD2mat(int_garb.symtext.mult_table, int_garb.irrm, salc1.irrep, salc1.i, salc1.r, int_garb.g, U)
    V = int_garb.stat_subgrps[salc2.atom]
    W = int_garb.stat_subgrps[salc3.atom]
    X = int_garb.stat_subgrps[salc4.atom]
    @timeit to "get_Rs" Rs, λRs = get_Rs(to, U, V, int_garb.symtext)
    @timeit to "get_Ss" Ss, λSs = get_Rs(to, W, X, int_garb.symtext)
    troot = 0.0
    for R in Rs
        RV = Int64[int_garb.symtext.mult_table[R,v] for v in V]
        Db = genD2mat(int_garb.symtext.mult_table, int_garb.irrm, salc2.irrep, salc2.i, salc2.r, int_garb.g, RV)
        Rb = int_garb.symtext.atom_map[salc2.atom,R]
        M = intersect(int_garb.stat_subgrps[salc1.atom], int_garb.stat_subgrps[int_garb.symtext.atom_map[salc2.atom, R]])
        for S in Ss
            N = intersect(int_garb.stat_subgrps[salc3.atom], int_garb.stat_subgrps[int_garb.symtext.atom_map[salc4.atom, S]])
            @timeit to "get_Ts" Ts, λTs = get_Rs(to, M, N, int_garb.symtext)
            for (Tidx,T) in enumerate(Ts)
                TW = Int64[int_garb.symtext.mult_table[T,w] for w in W]                                
                TSX = Int64[int_garb.symtext.mult_table[T,int_garb.symtext.mult_table[S,x]] for x in X]                                
                @timeit to "Build Dc and Dd" begin
                    Dc = genD2mat(int_garb.symtext.mult_table, int_garb.irrm, salc3.irrep, salc3.i, salc3.r, int_garb.g, TW)
                    Dd = genD2mat(int_garb.symtext.mult_table, int_garb.irrm, salc4.irrep, salc4.i, salc4.r, int_garb.g, TSX)
                end
                Tc = int_garb.symtext.atom_map[salc3.atom,T]
                TSd = int_garb.symtext.atom_map[salc4.atom,int_garb.symtext.mult_table[T,S]]
                @timeit to "AO int calc" out = ERI_2e4c(bset, abars[salc1.atom][salc1.bfxn][1], abars[Rb][salc2.bfxn][1], abars[Tc][salc3.bfxn][1], abars[TSd][salc4.bfxn][1])
                @timeit to "Gen. C mats" begin
                    Ca = genC2matp(int_garb.fxnmap, U, 2*(salc1.sh-1)+1, salc1.sh, salc1.ml)
                    Cb = genC2matp(int_garb.fxnmap, RV, 2*(salc2.sh-1)+1, salc2.sh, salc2.ml)
                    Cc = genC2matp(int_garb.fxnmap, TW, 2*(salc3.sh-1)+1, salc3.sh, salc3.ml)
                    Cd = genC2matp(int_garb.fxnmap, TSX, 2*(salc4.sh-1)+1, salc4.sh, salc4.ml)
                end
                @timeit to "Calc H and Contract" begin
                    @timeit to "Contract" @tensoropt trout[g1,g2,g3,g4] := Da[g1,u]*Db[g2,v]*Dc[g3,w]*Dd[g4,x]*Ca[u,ab]*Cb[v,bb]*Cc[w,cb]*Cd[x,db]*out[ab,bb,cb,db]
                    summa = 0.0
                    for g = 1:int_garb.g
                        summa += trout[g,g,g,g]
                    end
                end
                @timeit to "Int assignment" begin
                    summa *= (int_garb.irr_dims[salc1.irrep]*int_garb.irr_dims[salc2.irrep]*int_garb.irr_dims[salc3.irrep]*int_garb.irr_dims[salc4.irrep])/(λTs[Tidx]*(int_garb.g^4))
                    troot += summa*salc1.γ*salc2.γ*salc3.γ*salc4.γ
                end
            end
        end
    end
    return troot
end