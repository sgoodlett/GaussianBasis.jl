
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

struct IntCollect
    symtext::Molecules.Symmetry.CharacterTables.SymText
    irrm::Dict{String, Vector{AbstractArray{Float64}}}
    fxnmap::Vector{Vector{Matrix{Float64}}}
    stat_subgrps::Vector{Vector{Int64}}
    g::Int64
    irr_dims::Dict{String, Int64}
end

function Λfxn(intcol::IntCollect, salc, ab, p, R::Int64)
    out = 0.0
    irrm = intcol.irrm[salc.irrep]
    for U in intcol.stat_subgrps[salc.atom]
        RU = intcol.symtext.mult_table[R, U]
        out += intcol.fxnmap[RU][salc.sh][ab, salc.ml]*irrm[RU][p,salc.r]
    end
    return out * intcol.irr_dims[salc.irrep] / intcol.g
end

# Same function but R=1
function Λfxn(intcol::IntCollect, salc, ab, p, printstuff=false)
    out = 0.0
    irrm = intcol.irrm[salc.irrep]
    for U in intcol.stat_subgrps[salc.atom]
        out += intcol.fxnmap[U][salc.sh][ab, salc.ml]*irrm[U][p,salc.r]
    end
    return out * intcol.irr_dims[salc.irrep] / intcol.g
end

function Γfxn(intcol::IntCollect, salc1, salc2, ab, bb, R, λr, printstuff::Bool)
    out = 0.0
    for p = 1:intcol.irr_dims[salc1.irrep]
        l1 = Λfxn(intcol, salc1, ab, p, printstuff)
        l2 = Λfxn(intcol, salc2, bb, p, R)
        if printstuff
            println("p: $p, Λ: $l1, Λ(R): $l2")
        end
        out += Λfxn(intcol, salc1, ab, p) * Λfxn(intcol, salc2, bb, p, R)
    end
    return out * intcol.g/(intcol.irr_dims[salc1.irrep] * λr)
end

function get_IntCollect(mol, symtext, bset)
    irrm = eval(Meta.parse("Molecules.Symmetry.CharacterTables.irrm_"*symtext.pg))
    subgrps = Vector{Int64}[]
    for A in 1:length(mol)
        push!(subgrps, get_atom_subgroup(A, symtext))
    end
    aotoso, ignore_this, salcs, bfxn_map = GaussianBasis.SALCs.ProjectionOp(mol, bset, symtext)[2:5]
    return IntCollect(symtext, irrm, bfxn_map, subgrps, symtext.order, symtext.ctab.irrep_dims), salcs, aotoso
end

function get_Stalcs(mol, bset)
    salcs, bfxn_map = GaussianBasis.SALCs.ProjectionOp(mol, bset)[4:5]
    return salcs, bfxn_map
end

function build_abars(bset::BasisSet)
    abars = []
    sh_ctr = 0
    bfxn_ctr = 0
    for atom = 1:length(bset.atoms)
        abars2 = []
        for sh_idx = 1:bset.shells_per_atom[atom]
            l = bset[sh_idx+sh_ctr].l
            for l_idx = 1:2l+1
                push!(abars2, (sh_idx+sh_ctr, collect(1+bfxn_ctr:2l+1+bfxn_ctr)))
            end
            bfxn_ctr += 2l+1
        end
        sh_ctr += bset.shells_per_atom[atom]
        push!(abars, abars2)
    end
    return abars
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
    abars = build_abars(bset)
    int_garb, salcs, aotoso = get_IntCollect(mol, symtext, bset) end
    @timeit to "Ints hard" begin
    S_hard = int_type(bset)
    sym_S = transpose(aotoso) * S_hard * aotoso end
    #for i = 1:length(salcs)
    #    println("Irrep.: $(salcs[i].irrep), Atom: $(salcs[i].atom), Bsfxn: $(salcs[i].bfxn), sh, ml: $(salcs[i].sh), $(salcs[i].ml), i,r: $(salcs[i].i), $(salcs[i].r), Factor:$(salcs[i].γ)")
    #end
    println("Max dif. in salcs: ", findmax(broadcast(abs, check_salcs(salcs, aotoso)))[1])
    slength = length(salcs)
    S = zeros(Float64, slength, slength)
    @timeit to "My Ints" begin
    for s1 = 1:slength
        for s2 = s1:slength
            salc1 = salcs[s1]
            salc2 = salcs[s2]
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
    abars = build_abars(bset)
    int_garb, salcs, aotoso = get_IntCollect(mol, symtext, bset)
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

function stwintegrals_fromDirect(mol::String, basis::String, Halg::Int=1)
    mol, symtext = Molecules.Symmetry.CharacterTables.symtext_from_file(mol)
    return stwintegrals_fromDirect(mol, basis, symtext, Halg)
end

function bigDirect_setup(to, mol, basis, symtext)
    bset = BasisSet(basis, mol)
    @timeit to "Abars and IntGarb" begin
    abars = build_abars(bset)
    int_garb, salcs, aotoso = get_IntCollect(mol, symtext, bset)
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

function isequiv1(salc1, salc2)
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

function find_unique_salcs1(salcs)
    # Group partner functions together
    out = [[(1,salcs[1])]]
    groop = [1]
    for (sidx,salc) in enumerate(salcs[2:end])
        chk = false
        for done_salcs in out
            if isequiv1(salc, done_salcs[1][2])
                push!(done_salcs, (sidx+1,salc))
                push!(groop, done_salcs[1][1])
                chk = true
            end
        end
        if !chk
            push!(out, [(sidx+1,salc)])
            push!(groop, sidx+1)
        end
    end
    return out, groop
end

function stwintegrals_fromDirect(mol, basis, symtext, Halg::Int)
    to = TimerOutput()
    bset, int_garb, salcs, abars, aotoso = bigDirect_setup(to, mol, basis, symtext)
    #@timeit to "Ints the hard way" begin
    #    ERI_AO = ERI_2e4c(bset)
    #    SERI = zeros(bset.nbas, bset.nbas, bset.nbas, bset.nbas)
    #    @tensoropt SERI[i,j,k,l] =  aotoso[μ, i]*aotoso[ν, j]*ERI_AO[μ, ν, ρ, σ]*aotoso[ρ, k]*aotoso[σ, l]
    #end
    slength = length(salcs)
    eri = zeros(Float64, slength, slength, slength, slength)
    #println(salcs)
    unique_salcs, groops = find_unique_salcs1(salcs)
    #display(salcs[29:end])
    #display(unique_salcs[end])
    #println(groops)
    #println(get_petite_idxs(salcs))
    #println(length(get_petite_idxs(salcs)))
    #cg = build_ClebschGordan(int_garb)
    #return cg
    load_bar = "Calculating Integrals\n"
    if Halg < 10
        load_bar *= "-"^slength
    else
        load_bar *= "-"^length(unique_salcs)
    end
    @timeit to "My Ints" begin
    println(load_bar)
    if Halg < 10
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
                        if (Halg%10) == 1
                            @timeit to "Direct Symint" troot = stwintegrals_bigDirect(to, bset, salcs, int_garb, abars, s1, s2, s3, s4)
                        elseif (Halg%10) == 2
                            @timeit to "Direct Symint" troot = stwintegrals_bigDirect_2(to, bset, salcs, int_garb, abars, s1, s2, s3, s4)
                        else
                            X = compute_X(to, bset, salcs, int_garb, abars, s1, s2, s3, s4)
                            troot = 0.0
                            salc1 = salcs[s1]
                            salc2 = salcs[s2]
                            salc3 = salcs[s3]
                            salc4 = salcs[s4]
                            if Halg == 3
                                D1 = build_Darr(int_garb, salcs[s1].irrep, salc1.i)
                                D2 = build_Darr(int_garb, salcs[s2].irrep, salc2.i)
                                D3 = build_Darr(int_garb, salcs[s3].irrep, salc3.i)
                                D4 = build_Darr(int_garb, salcs[s4].irrep, salc4.i)
                                suum = zeros(Float64, (int_garb.irr_dims[salc1.irrep],int_garb.irr_dims[salc2.irrep],int_garb.irr_dims[salc3.irrep],int_garb.irr_dims[salc4.irrep]))
                                for g = 1:int_garb.g
                                    suum += trout[g,g,g,g]
                                end
                                @timeit to "Contract" @tensoropt trout[g1,g2,g3,g4] := D1[g1,ib]*D2[g2,jb]*D3[g3,kb]*D4[g4,lb]*X[ib,jb,kb,lb]
                                #suum = 0.0
                                #for g = 1:int_garb.g
                                #    suum += trout[g,g,g,g]
                                #end
                                troot += suum
                            else
                                @timeit to "Contract" begin
                                for ib = 1:int_garb.irr_dims[salc1.irrep]
                                    for jb = 1:int_garb.irr_dims[salc2.irrep]
                                        for kb = 1:int_garb.irr_dims[salc3.irrep]
                                            for lb = 1:int_garb.irr_dims[salc4.irrep]
                                                trout = 0.0
                                                for g = 1:int_garb.g
                                                    trout += int_garb.irrm[salc1.irrep][g][salc1.i,ib]*int_garb.irrm[salc2.irrep][g][salc2.i,jb]*int_garb.irrm[salc3.irrep][g][salc3.i,kb]*int_garb.irrm[salc4.irrep][g][salc4.i,lb]
                                                end
                                                troot += trout*X[ib,jb,kb,lb]
                                            end
                                        end
                                    end
                                end
                                end
                            end
                            troot *= salc1.γ * salc2.γ * salc3.γ * salc4.γ
                        end
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
    elseif Halg > 10
        for s1_tup in unique_salcs
            s1 = s1_tup[1][1]
            print("*")
            for s2_tup in unique_salcs
                s2 = s2_tup[1][1]
                s1 > s2 ? continue : nothing
                for s3_tup in unique_salcs
                    s3 = s3_tup[1][1]
                    for s4_tup in unique_salcs
                        s4 = s4_tup[1][1]
                        s3 > s4 ? continue : nothing
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
                        if Halg >= 15
                            # Compute D(G,i,ib) tensors
                            # Compute X(ib,jb,kb,lb) tensor
                            # I(i,j,k,l) = sum over G: contract A(G1,G2,G3,G4) = Da(G1,i,ib)Db(G2,j,jb)Dc(G3,k,kb)Dd(G4,l,lb)X(ib,jb,kb,lb)
                            salc1 = salcs[s1]
                            salc2 = salcs[s2]
                            salc3 = salcs[s3]
                            salc4 = salcs[s4]
                            s1irrdim = int_garb.irr_dims[salc1.irrep]
                            s2irrdim = int_garb.irr_dims[salc2.irrep]
                            s3irrdim = int_garb.irr_dims[salc3.irrep]
                            s4irrdim = int_garb.irr_dims[salc4.irrep]
                            D1 = int_garb.irrm[salc1.irrep]
                            D2 = int_garb.irrm[salc2.irrep]
                            D3 = int_garb.irrm[salc3.irrep]
                            D4 = int_garb.irrm[salc4.irrep]
                            troot = zeros(Float64, (s1irrdim, s2irrdim, s3irrdim, s4irrdim))
                            #troot = zeros(Float64, (int_garb.irr_dims[salcs[s1].irrep], int_garb.irr_dims[salcs[s2].irrep], int_garb.irr_dims[salcs[s3].irrep], int_garb.irr_dims[salcs[s4].irrep]))
                            @timeit to "Compute X" X = compute_X(to, bset, salcs, int_garb, abars, s1, s2, s3, s4)
                            if Halg == 15
                                for i = 1:s1irrdim
                                    for ib = 1:s1irrdim
                                        for j = 1:s2irrdim
                                            for jb = 1:s2irrdim
                                                for k = 1:s3irrdim
                                                    for kb = 1:s3irrdim
                                                        for l = 1:s4irrdim
                                                            for lb = 1:s4irrdim
                                                                @timeit to "Contract over G" begin
                                                                trout = 0.0
                                                                for g = 1:int_garb.g
                                                                    #trout += int_garb.irrm[salc1.irrep][g][i,ib]*int_garb.irrm[salc2.irrep][g][j,jb]*int_garb.irrm[salc3.irrep][g][k,kb]*int_garb.irrm[salc4.irrep][g][l,lb]
                                                                    trout += D1[g][i,ib]*D2[g][j,jb]*D3[g][k,kb]*D4[g][l,lb]
                                                                end
                                                                end
                                                                troot[i,j,k,l] += trout*X[ib,jb,kb,lb]
                                                            end
                                                        end
                                                    end
                                                #troot[i,j,k,l] *= salcs[s1].γ * salcs[s2].γ * salcs[s3].γ * salcs[s4].γ
                                                end
                                            end
                                        end
                                    end
                                end
                                #if s1 == 4 && s2 == 4 && s3 == 4 && s4 ==4
                                #    println(troot[1,1,1,1])
                                #end
                            elseif Halg == 16
                                @timeit to "Contract over G" begin
                                if s1irrdim > 1
                                    if s2irrdim > 1
                                        if s3irrdim > 1
                                            if s4irrdim > 1
                                                for g = 1:int_garb.g
                                                    Di = D1[g]
                                                    Dj = D2[g]
                                                    Dk = D3[g]
                                                    Dl = D4[g]
                                                    @tensoropt troot[i,j,k,l] += Di[i,ib]*Dj[j,jb]*Dk[k,kb]*Dl[l,lb]*X[ib,jb,kb,lb]
                                                end
                                            else
                                                for g = 1:int_garb.g
                                                    Di = D1[g]
                                                    Dj = D2[g]
                                                    Dk = D3[g]
                                                    Dl = D4[g][1,1]
                                                    @tensoropt troot[i,j,k,l] += Di[i,ib]*Dj[j,jb]*Dk[k,kb]*Dl*X[ib,jb,kb,l]
                                                end
                                            end
                                        else
                                            if s4irrdim > 1
                                                for g = 1:int_garb.g
                                                    Di = D1[g]
                                                    Dj = D2[g]
                                                    Dk = D3[g][1,1]
                                                    Dl = D4[g]
                                                    @tensoropt troot[i,j,k,l] += Di[i,ib]*Dj[j,jb]*Dk*Dl[l,lb]*X[ib,jb,k,lb]
                                                end
                                            else
                                                for g = 1:int_garb.g
                                                    Di = D1[g]
                                                    Dj = D2[g]
                                                    Dk = D3[g][1,1]
                                                    Dl = D4[g][1,1]
                                                    @tensoropt troot[i,j,k,l] += Di[i,ib]*Dj[j,jb]*Dk*Dl*X[ib,jb,k,l]
                                                end
                                            end
                                        end
                                    else
                                        if s3irrdim > 1
                                            if s4irrdim > 1
                                                for g = 1:int_garb.g
                                                    Di = D1[g]
                                                    Dj = D2[g][1,1]
                                                    Dk = D3[g]
                                                    Dl = D4[g]
                                                    @tensoropt troot[i,j,k,l] += Di[i,ib]*Dj*Dk[k,kb]*Dl[l,lb]*X[ib,j,kb,lb]
                                                end
                                            else
                                                for g = 1:int_garb.g
                                                    Di = D1[g]
                                                    Dj = D2[g][1,1]
                                                    Dk = D3[g]
                                                    Dl = D4[g][1,1]
                                                    @tensoropt troot[i,j,k,l] += Di[i,ib]*Dj*Dk[k,kb]*Dl*X[ib,j,kb,l]
                                                end
                                            end
                                        else
                                            if s4irrdim > 1
                                                for g = 1:int_garb.g
                                                    Di = D1[g]
                                                    Dj = D2[g][1,1]
                                                    Dk = D3[g][1,1]
                                                    Dl = D4[g]
                                                    @tensoropt troot[i,j,k,l] += Di[i,ib]*Dj*Dk*Dl[l,lb]*X[ib,j,k,lb]
                                                end
                                            else
                                                for g = 1:int_garb.g
                                                    Di = D1[g]
                                                    Dj = D2[g][1,1]
                                                    Dk = D3[g][1,1]
                                                    Dl = D4[g][1,1]
                                                    @tensoropt troot[i,j,k,l] += Di[i,ib]*Dj*Dk*Dl*X[ib,j,k,l]
                                                end
                                            end
                                        end
                                    end
                                else
                                    if s2irrdim > 1
                                        if s3irrdim > 1
                                            if s4irrdim > 1
                                                for g = 1:int_garb.g
                                                    Di = D1[g][1,1]
                                                    Dj = D2[g]
                                                    Dk = D3[g]
                                                    Dl = D4[g]
                                                    @tensoropt troot[i,j,k,l] += Di*Dj[j,jb]*Dk[k,kb]*Dl[l,lb]*X[i,jb,kb,lb]
                                                end
                                            else
                                                for g = 1:int_garb.g
                                                    Di = D1[g][1,1]
                                                    Dj = D2[g]
                                                    Dk = D3[g]
                                                    Dl = D4[g][1,1]
                                                    @tensoropt troot[i,j,k,l] += Di*Dj[j,jb]*Dk[k,kb]*Dl*X[i,jb,kb,l]
                                                end
                                            end
                                        else
                                            if s4irrdim > 1
                                                for g = 1:int_garb.g
                                                    Di = D1[g][1,1]
                                                    Dj = D2[g]
                                                    Dk = D3[g][1,1]
                                                    Dl = D4[g]
                                                    @tensoropt troot[i,j,k,l] += Di*Dj[j,jb]*Dk*Dl[l,lb]*X[i,jb,k,lb]
                                                end
                                            else
                                                for g = 1:int_garb.g
                                                    Di = D1[g][1,1]
                                                    Dj = D2[g]
                                                    Dk = D3[g][1,1]
                                                    Dl = D4[g][1,1]
                                                    @tensoropt troot[i,j,k,l] += Di*Dj[j,jb]*Dk*Dl*X[i,jb,k,l]
                                                end
                                            end
                                        end
                                    else
                                        if s3irrdim > 1
                                            if s4irrdim > 1
                                                for g = 1:int_garb.g
                                                    Di = D1[g][1,1]
                                                    Dj = D2[g][1,1]
                                                    Dk = D3[g]
                                                    Dl = D4[g]
                                                    @tensoropt troot[i,j,k,l] += Di*Dj*Dk[k,kb]*Dl[l,lb]*X[i,j,kb,lb]
                                                end
                                            else
                                                for g = 1:int_garb.g
                                                    Di = D1[g][1,1]
                                                    Dj = D2[g][1,1]
                                                    Dk = D3[g]
                                                    Dl = D4[g][1,1]
                                                    @tensoropt troot[i,j,k,l] += Di*Dj*Dk[k,kb]*Dl*X[i,j,kb,l]
                                                end
                                            end
                                        else
                                            if s4irrdim > 1
                                                for g = 1:int_garb.g
                                                    Di = D1[g][1,1]
                                                    Dj = D2[g][1,1]
                                                    Dk = D3[g][1,1]
                                                    Dl = D4[g]
                                                    @tensoropt troot[i,j,k,l] += Di*Dj*Dk*Dl[l,lb]*X[i,j,k,lb]
                                                end
                                            else
                                                for g = 1:int_garb.g
                                                    Di = D1[g][1,1]
                                                    Dj = D2[g][1,1]
                                                    Dk = D3[g][1,1]
                                                    Dl = D4[g][1,1]
                                                    troot[1,1,1,1] += Di*Dj*Dk*Dl*X[1,1,1,1]
                                                end
                                            end
                                        end
                                    end
                                end
                                end

                            elseif Halg == 17
                                for ib = 1:s1irrdim
                                    for jb = 1:s2irrdim
                                        for kb = 1:s3irrdim
                                            for lb = 1:s4irrdim
                                                @timeit to "Contract over G" begin
                                                trout = 0.0
                                                for g = 1:int_garb.g
                                                    trout += D1[g][salc1.i,ib]*D2[g][salc2.i,jb]*D3[g][salc3.i,kb]*D4[g][salc4.i,lb]
                                                end
                                                end
                                                troot[salc1.i,salc2.i,salc3.i,salc4.i] += trout*X[ib,jb,kb,lb]
                                            end
                                        end
                                    end
                                end
                                @timeit to "Distribute to subspecies" begin
                                ref = troot[salc1.i,salc2.i,salc3.i,salc4.i]
                                for i = 1:s1irrdim
                                    for j = 1:s2irrdim
                                        for k = 1:s3irrdim
                                            for l = 1:s4irrdim
                                                boot = 0.0
                                                for g = 1:int_garb.g
                                                    boot += D1[g][salc1.i,i]*D2[g][salc2.i,j]*D3[g][salc3.i,k]*D4[g][salc4.i,l]
                                                end
                                                troot[i,j,k,l] = boot*ref
                                            end
                                        end
                                    end
                                end
                                #troot *= sqrt(s1irrdim*s2irrdim*s3irrdim*s4irrdim)/int_garb.g
                                #troot *= s1irrdim/int_garb.g
                                if true #s1 == s2 && s3 == s4 && s1 == s3
                                    troot *= s1irrdim / int_garb.g
                                else
                                    troot *= sqrt(s1irrdim*s2irrdim*s3irrdim*s4irrdim)/int_garb.g
                                end
                            end
                            end
                            #troot .*= salcs[s1].γ * salcs[s2].γ * salcs[s3].γ * salcs[s4].γ
                            @timeit to "Int Assignment" begin
                            for (s1_d_idx,s1_d) in enumerate(s1_tup)
                                for (s2_d_idx,s2_d) in enumerate(s2_tup)
                                    for (s3_d_idx,s3_d) in enumerate(s3_tup)
                                        for (s4_d_idx,s4_d) in enumerate(s4_tup)
                                            s1_d[1] > s2_d[1] ? continue : nothing
                                            s3_d[1] > s4_d[1] ? continue : nothing
                                            i12 = GaussianBasis.index2(groops[s1_d[1]],groops[s2_d[1]])
                                            i34 = GaussianBasis.index2(groops[s3_d[1]],groops[s4_d[1]])
                                            i12 < i34 ? continue : nothing
                                            troots = troot[s1_d_idx,s2_d_idx,s3_d_idx,s4_d_idx]
                                            troots *= s1_d[2].γ * s2_d[2].γ * s3_d[2].γ * s4_d[2].γ
                                            eri[s1_d[1],s2_d[1],s3_d[1],s4_d[1]] += troots
                                            if s1_d[1] != s2_d[1]
                                                eri[s2_d[1],s1_d[1],s3_d[1],s4_d[1]] += troots
                                                if s3_d[1] != s4_d[1]
                                                    eri[s1_d[1],s2_d[1],s4_d[1],s3_d[1]] += troots
                                                    eri[s2_d[1],s1_d[1],s4_d[1],s3_d[1]] += troots
                                                    if i12 != i34
                                                        eri[s3_d[1],s4_d[1],s1_d[1],s2_d[1]] += troots
                                                        eri[s3_d[1],s4_d[1],s2_d[1],s1_d[1]] += troots
                                                        eri[s4_d[1],s3_d[1],s1_d[1],s2_d[1]] += troots
                                                        eri[s4_d[1],s3_d[1],s2_d[1],s1_d[1]] += troots
                                                    end
                                                elseif i12 != i34
                                                    eri[s3_d[1],s4_d[1],s1_d[1],s2_d[1]] += troots
                                                    eri[s3_d[1],s4_d[1],s2_d[1],s1_d[1]] += troots
                                                end
                                            elseif s3_d[1] != s4_d[1]
                                                eri[s1_d[1],s2_d[1],s4_d[1],s3_d[1]] += troots
                                                if i12 != i34
                                                    eri[s3_d[1],s4_d[1],s1_d[1],s2_d[1]] += troots
                                                    eri[s4_d[1],s3_d[1],s1_d[1],s2_d[1]] += troots
                                                end
                                            elseif i12 != i34
                                                eri[s3_d[1],s4_d[1],s1_d[1],s2_d[1]] += troots
                                            end
                                        end
                                    end
                                end
                            end
                            end
                        end
                    end
                end
            end
        end
        end
    end
    print("\n")
    #println(findmax(broadcast(abs, eri-SERI)))
    #println(eri[9,6,7,7])
    #println(broadcast(abs, eri[4:6,4:6,7:9,7:9]-SERI[4:6,4:6,7:9,7:9]))
    show(to)
end

function RestrictedIntegrals(mol::String, basis::String)
    mol, symtext = Molecules.Symmetry.CharacterTables.symtext_from_file(mol)
    return RestrictedIntegrals(mol, basis, symtext)
end

function RestrictedIntegrals(mol, basis, symtext)
    to = TimerOutput()
    bset, int_garb, salcs, abars, aotoso = bigDirect_setup(to, mol, basis, symtext)
    @timeit to "Ints the hard way" begin
        ERI_AO = ERI_2e4c(bset)
        SERI = zeros(bset.nbas, bset.nbas, bset.nbas, bset.nbas)
        @tensoropt SERI[i,j,k,l] =  aotoso[μ, i]*aotoso[ν, j]*ERI_AO[μ, ν, ρ, σ]*aotoso[ρ, k]*aotoso[σ, l]
    end
    slength = length(salcs)
    eri = zeros(Float64, slength, slength, slength, slength)
    unique_salcs, groops = find_unique_salcs1(salcs)
    display(unique_salcs[end])
    println(groops)
    load_bar = "Calculating Integrals\n"
    load_bar *= "-"^length(unique_salcs)
    @timeit to "My Ints" begin
    println(load_bar)
    for s1_tup in unique_salcs
        s1 = s1_tup[1][1]
        print("*")
        for s2_tup in unique_salcs
            s2 = s2_tup[1][1]
            s1 > s2 ? continue : nothing
            for s3_tup in unique_salcs
                s3 = s3_tup[1][1]
                for s4_tup in unique_salcs
                    s4 = s4_tup[1][1]
                    s3 > s4 ? continue : nothing
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
                    # Compute D(G,i,ib) tensors
                    # Compute X(ib,jb,kb,lb) tensor
                    # I(i,j,k,l) = sum over G: contract A(G1,G2,G3,G4) = Da(G1,i,ib)Db(G2,j,jb)Dc(G3,k,kb)Dd(G4,l,lb)X(ib,jb,kb,lb)
                    salc1 = salcs[s1]
                    salc2 = salcs[s2]
                    salc3 = salcs[s3]
                    salc4 = salcs[s4]
                    s1irrdim = int_garb.irr_dims[salc1.irrep]
                    s2irrdim = int_garb.irr_dims[salc2.irrep]
                    s3irrdim = int_garb.irr_dims[salc3.irrep]
                    s4irrdim = int_garb.irr_dims[salc4.irrep]
                    D1 = int_garb.irrm[salc1.irrep]
                    D2 = int_garb.irrm[salc2.irrep]
                    D3 = int_garb.irrm[salc3.irrep]
                    D4 = int_garb.irrm[salc4.irrep]
                    troot = zeros(Float64, (s1irrdim, s2irrdim, s3irrdim, s4irrdim))
                    #troot = zeros(Float64, (int_garb.irr_dims[salcs[s1].irrep], int_garb.irr_dims[salcs[s2].irrep], int_garb.irr_dims[salcs[s3].irrep], int_garb.irr_dims[salcs[s4].irrep]))
                    @timeit to "Compute X" X = compute_X(to, bset, salcs, int_garb, abars, s1, s2, s3, s4)
                    @timeit to "Compute J" begin
                    J = 0.0
                    if salc1.irrep == salc2.irrep && salc1.i == salc2.i && salc3.irrep == salc4.irrep && salc3.i == salc4.i
                        for ib = 1:s1irrdim
                            for kb = 1:s3irrdim
                                J += X[ib,ib,kb,kb]
                            end
                        end
                        J *= 1.0/(s2irrdim)
                        end
                    end
                    @timeit to "Compute K" begin
                    if salc1.irrep == salc4.irrep && salc1.i == salc4.i && salc2.irrep == salc3.irrep && salc2.i == salc3.i
                        J = 0.0
                        for ib = 1:s1irrdim
                            for jb = 1:s2irrdim
                                J += X[ib,jb,jb,ib]
                            end
                        end
                        J *= 1.0/(s4irrdim)
                        end
                    end
                    @timeit to "Int Assignment" begin
                    for (s1_d_idx,s1_d) in enumerate(s1_tup)
                        for (s2_d_idx,s2_d) in enumerate(s2_tup)
                            for (s3_d_idx,s3_d) in enumerate(s3_tup)
                                for (s4_d_idx,s4_d) in enumerate(s4_tup)
                                    s1_d[1] > s2_d[1] ? continue : nothing
                                    s3_d[1] > s4_d[1] ? continue : nothing
                                    i12 = GaussianBasis.index2(groops[s1_d[1]],groops[s2_d[1]])
                                    i34 = GaussianBasis.index2(groops[s3_d[1]],groops[s4_d[1]])
                                    i12 < i34 ? continue : nothing
                                    troots = troot[s1_d_idx,s2_d_idx,s3_d_idx,s4_d_idx]
                                    troots *= s1_d[2].γ * s2_d[2].γ * s3_d[2].γ * s4_d[2].γ
                                    eri[s1_d[1],s2_d[1],s3_d[1],s4_d[1]] += troots
                                    if s1_d[1] != s2_d[1]
                                        eri[s2_d[1],s1_d[1],s3_d[1],s4_d[1]] += troots
                                        if s3_d[1] != s4_d[1]
                                            eri[s1_d[1],s2_d[1],s4_d[1],s3_d[1]] += troots
                                            eri[s2_d[1],s1_d[1],s4_d[1],s3_d[1]] += troots
                                            if i12 != i34
                                                eri[s3_d[1],s4_d[1],s1_d[1],s2_d[1]] += troots
                                                eri[s3_d[1],s4_d[1],s2_d[1],s1_d[1]] += troots
                                                eri[s4_d[1],s3_d[1],s1_d[1],s2_d[1]] += troots
                                                eri[s4_d[1],s3_d[1],s2_d[1],s1_d[1]] += troots
                                            end
                                        elseif i12 != i34
                                            eri[s3_d[1],s4_d[1],s1_d[1],s2_d[1]] += troots
                                            eri[s3_d[1],s4_d[1],s2_d[1],s1_d[1]] += troots
                                        end
                                    elseif s3_d[1] != s4_d[1]
                                        eri[s1_d[1],s2_d[1],s4_d[1],s3_d[1]] += troots
                                        if i12 != i34
                                            eri[s3_d[1],s4_d[1],s1_d[1],s2_d[1]] += troots
                                            eri[s4_d[1],s3_d[1],s1_d[1],s2_d[1]] += troots
                                        end
                                    elseif i12 != i34
                                        eri[s3_d[1],s4_d[1],s1_d[1],s2_d[1]] += troots
                                    end
                                end
                            end
                        end
                    end
                    end
                end
            end
        end
    end
    end
    print("\n")
    println(findmax(broadcast(abs, eri-SERI)))
    #println(eri[9,6,7,7])
    #println(broadcast(abs, eri[4:6,4:6,7:9,7:9]-SERI[4:6,4:6,7:9,7:9]))
    show(to)
end

function stwintegrals_bigDirect(to, bset, salcs, int_garb, abars, s1, s2, s3, s4)
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
function build_Λarr(int_garb, salc, abidx_l, ib_l,R)
    return Float64[Λfxn(int_garb, salc, ab, ib, R) for ab=1:abidx_l, ib=1:ib_l]
end

function build_Λarr(int_garb,salc,abidx_l,ib_l)
    return Float64[Λfxn(int_garb, salc, ab, ib) for ab=1:abidx_l, ib=1:ib_l]
end

function build_Darr(int_garb, irrep, i)
    return transpose(reduce(hcat,[int_garb.irrm[irrep][g][i,:] for g=1:int_garb.g]))
end

function stwintegrals_bigDirect_2(to, bset, salcs, int_garb, abars, s1, s2, s3, s4)
    salc1 = salcs[s1]
    salc2 = salcs[s2]
    salc3 = salcs[s3]
    salc4 = salcs[s4]
    U = int_garb.stat_subgrps[salc1.atom]
    V = int_garb.stat_subgrps[salc2.atom]
    W = int_garb.stat_subgrps[salc3.atom]
    X = int_garb.stat_subgrps[salc4.atom]
    s1irrdim = int_garb.irr_dims[salc1.irrep]
    s2irrdim = int_garb.irr_dims[salc2.irrep]
    s3irrdim = int_garb.irr_dims[salc3.irrep]
    s4irrdim = int_garb.irr_dims[salc4.irrep]
    D1 = build_Darr(int_garb, salc1.irrep, salc1.i)
    D2 = build_Darr(int_garb, salc2.irrep, salc2.i)
    D3 = build_Darr(int_garb, salc3.irrep, salc3.i)
    D4 = build_Darr(int_garb, salc4.irrep, salc4.i)
    Λ1 = build_Λarr(int_garb, salc1, 2*salc1.sh-1, s1irrdim)
    @timeit to "get_Rs" Rs, λRs = get_Rs(to, U, V, int_garb.symtext)
    @timeit to "get_Ss" Ss, λSs = get_Rs(to, W, X, int_garb.symtext)
    troot = 0.0
    for R in Rs
        Λ2 = build_Λarr(int_garb, salc2, 2*salc2.sh-1, s2irrdim, R)
        M = intersect(int_garb.stat_subgrps[salc1.atom], int_garb.stat_subgrps[int_garb.symtext.atom_map[salc2.atom, R]])
        Rb = int_garb.symtext.atom_map[salc2.atom,R]
        for S in Ss
            N = intersect(int_garb.stat_subgrps[salc3.atom], int_garb.stat_subgrps[int_garb.symtext.atom_map[salc4.atom, S]])
            @timeit to "get_Ts" Ts, λTs = get_Rs(to, M, N, int_garb.symtext)
            for (Tidx,T) in enumerate(Ts)
                Tc = int_garb.symtext.atom_map[salc3.atom,T]
                TSd = int_garb.symtext.atom_map[salc4.atom,int_garb.symtext.mult_table[T,S]]
                @timeit to "AO int calc" out = ERI_2e4c(bset, abars[salc1.atom][salc1.bfxn][1], abars[Rb][salc2.bfxn][1], abars[Tc][salc3.bfxn][1], abars[TSd][salc4.bfxn][1])
                @timeit to "Lambda build" begin
                Λ3 = build_Λarr(int_garb, salc3, 2*salc3.sh-1, s3irrdim, T)
                Λ4 = build_Λarr(int_garb, salc4, 2*salc4.sh-1, s4irrdim, int_garb.symtext.mult_table[T,S])
                end
                @timeit to "Calc H and Contract" begin
                    @timeit to "Contract" @tensoropt trout[g1,g2,g3,g4] := D1[g1,ib]*D2[g2,jb]*D3[g3,kb]*D4[g4,lb]*Λ1[ab,ib]*Λ2[bb,jb]*Λ3[cb,kb]*Λ4[db,lb]*out[ab,bb,cb,db]
                    summa = 0.0
                    for g = 1:int_garb.g
                        summa += trout[g,g,g,g]
                    end
                end
                @timeit to "Int assignment" begin
                    troot += summa*salc1.γ*salc2.γ*salc3.γ*salc4.γ/λTs[Tidx]
                end
            end
        end
    end
    return troot
end

function compute_X(to, bset, salcs, int_garb, abars, s1, s2, s3, s4)
    salc1 = salcs[s1]
    salc2 = salcs[s2]
    salc3 = salcs[s3]
    salc4 = salcs[s4]
    U = int_garb.stat_subgrps[salc1.atom]
    V = int_garb.stat_subgrps[salc2.atom]
    W = int_garb.stat_subgrps[salc3.atom]
    X = int_garb.stat_subgrps[salc4.atom]
    s1irrdim = int_garb.irr_dims[salc1.irrep]
    s2irrdim = int_garb.irr_dims[salc2.irrep]
    s3irrdim = int_garb.irr_dims[salc3.irrep]
    s4irrdim = int_garb.irr_dims[salc4.irrep]
    Xp = zeros(Float64, (s1irrdim, s2irrdim, s3irrdim, s4irrdim))
    Λ1 = build_Λarr(int_garb, salc1, 2*salc1.sh-1, s1irrdim)
    @timeit to "get_Rs" Rs, λRs = get_Rs(to, U, V, int_garb.symtext)
    @timeit to "get_Ss" Ss, λSs = get_Rs(to, W, X, int_garb.symtext)
    for R in Rs
        Λ2 = build_Λarr(int_garb, salc2, 2*salc2.sh-1, s2irrdim, R)
        M = intersect(int_garb.stat_subgrps[salc1.atom], int_garb.stat_subgrps[int_garb.symtext.atom_map[salc2.atom, R]])
        Rb = int_garb.symtext.atom_map[salc2.atom,R]
        for S in Ss
            N = intersect(int_garb.stat_subgrps[salc3.atom], int_garb.stat_subgrps[int_garb.symtext.atom_map[salc4.atom, S]])
            @timeit to "get_Ts" Ts, λTs = get_Rs(to, M, N, int_garb.symtext)
            for (Tidx,T) in enumerate(Ts)
                Tc = int_garb.symtext.atom_map[salc3.atom,T]
                TSd = int_garb.symtext.atom_map[salc4.atom,int_garb.symtext.mult_table[T,S]]
                @timeit to "AO int calc" out = ERI_2e4c(bset, abars[salc1.atom][salc1.bfxn][1], abars[Rb][salc2.bfxn][1], abars[Tc][salc3.bfxn][1], abars[TSd][salc4.bfxn][1])
                @timeit to "Lambda build" begin
                Λ3 = build_Λarr(int_garb, salc3, 2*salc3.sh-1, s3irrdim, T)
                Λ4 = build_Λarr(int_garb, salc4, 2*salc4.sh-1, s4irrdim, int_garb.symtext.mult_table[T,S])
                end
                λTi = 1/λTs[Tidx]
                @timeit to "Contract for X" @tensoropt Xp[ib,jb,kb,lb] += λTi * Λ1[ab,ib]*Λ2[bb,jb]*Λ3[cb,kb]*Λ4[db,lb]*out[ab,bb,cb,db]
            end
        end
    end
    return Xp
end