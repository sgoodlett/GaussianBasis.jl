# Shell is probs a bad name, actually bfxn idx
struct Stalc
    coeffs::Vector{Float64}
    irrep::String
    atom::Int64
    sh::Int64
    i::Int64
    j::Int64
    γ::Float64
end

function get_atom_subgroup(atomidx, symtext)
    subgroup = []
    for (i,v) in enumerate(symtext.atom_map[atomidx,:])
        if v == atomidx
            push!(subgroup, i)
        end
    end
    return subgroup
end

function get_double_cosets(U, V, symtext)
    double_cosets = []
    for g = 1:length(symtext.symels)
        double_coset = []
        for u in U, v in V
            push!(double_coset, symtext.mult_table[u,symtext.mult_table[g,v]])
        end
        push!(double_cosets, double_coset)
    end
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
    return unique_double_cosets
end

function get_Rs(U, V, symtext)
    Rs = []
    λs = []
    #U = get_atom_subgroup(atom1idx, symtext)
    #V = get_atom_subgroup(atom2idx, symtext)
    double_cosets = get_double_cosets(U, V, symtext)
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

struct IntGarbage
    symtext
    irrm
    fxnmap
    stat_subgrps
    g
    irr_dims
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
    subgrps = []
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
    end
    @timeit to "Ints the hard way" begin
    ERI_AO = ERI_2e4c(bset)
    SERI = zeros(bset.nbas, bset.nbas, bset.nbas, bset.nbas)
    @tensoropt SERI[i,j,k,l] =  aotoso[μ, i]*aotoso[ν, j]*ERI_AO[μ, ν, ρ, σ]*aotoso[ρ, k]*aotoso[σ, l]
    end
    for i = 1:length(salcs)
        println("Irrep.: $(salcs[i].irrep), Atom: $(salcs[i].atom), Bsfxn: $(salcs[i].bfxn), sh, ml: $(salcs[i].sh), $(salcs[i].ml), i,r: $(salcs[i].i), $(salcs[i].r)")
    end
    println(check_salcs(salcs, aotoso))
    porque = zeros(Float64, size(aotoso))
    for i = 1:length(salcs)
        porque[:,i] = aotoso[:,i]-salcs[i].coeffs
    end
    println(findmax(broadcast(abs, porque)))
    slength = length(salcs)
    eri = zeros(Float64, slength, slength, slength, slength)
    @timeit to "My Ints" begin
    for (s1,salc1) = enumerate(salcs)
        U = int_garb.stat_subgrps[salc1.atom]
        println("*Salc 1: $s1")
        for (s2,salc2) = enumerate(salcs)
            V = int_garb.stat_subgrps[salc2.atom]
            @timeit to "get_Rs" Rs, λRs = get_Rs(U, V, int_garb.symtext)
            for R in Rs
                M = intersect(int_garb.stat_subgrps[salc1.atom], int_garb.stat_subgrps[int_garb.symtext.atom_map[salc2.atom, R]])
                for (s3,salc3) = enumerate(salcs)
                    W = int_garb.stat_subgrps[salc3.atom]
                    for (s4,salc4) = enumerate(salcs)
                        if !tran_as_A1([salc1.irrep, salc2.irrep, salc3.irrep, salc4.irrep], symtext)
                            continue
                        end
                        X = int_garb.stat_subgrps[salc4.atom]
                        @timeit to "get_Ss" Ss, λSs = get_Rs(W, X, int_garb.symtext)
                        for S in Ss
                            N = intersect(int_garb.stat_subgrps[salc3.atom], int_garb.stat_subgrps[int_garb.symtext.atom_map[salc4.atom, S]])
                            @timeit to "get_Ts" Ts, λTs = get_Rs(M, N, int_garb.symtext)
                            for (Tidx,T) in enumerate(Ts)
                                for (abidx,ab) in enumerate(abars[salc1.atom][salc1.bfxn][2])
                                    for (bbidx,bb) in enumerate(abars[salc2.atom][salc2.bfxn][2])
                                        for (cbidx,cb) in enumerate(abars[salc3.atom][salc3.bfxn][2])
                                            for (dbidx,db) in enumerate(abars[salc4.atom][salc4.bfxn][2])
    # The fun begins
    # This looks really weird but I don't like how far to the right we've gone...
    @timeit to "Calc H" begin
    suum = 0.0
    for v in V
        rv = int_garb.symtext.mult_table[R,v]
        for w in W
            tw = int_garb.symtext.mult_table[T,w]
            for x in X
                tsx = int_garb.symtext.mult_table[T,int_garb.symtext.mult_table[S,x]]
                for u in U
                    @timeit to "Sum over G" begin
                    suuum = 0.0
                    for g = 1:int_garb.g
                        gu = int_garb.symtext.mult_table[g,u]
                        grv = int_garb.symtext.mult_table[g,rv]
                        gtw = int_garb.symtext.mult_table[g,tw]
                        gtsx = int_garb.symtext.mult_table[g,tsx]
                        suuum += int_garb.irrm[salc1.irrep][gu][salc1.i,salc1.r]*int_garb.irrm[salc2.irrep][grv][salc2.i,salc2.r]*int_garb.irrm[salc3.irrep][gtw][salc3.i,salc3.r]*int_garb.irrm[salc4.irrep][gtsx][salc4.i,salc4.r]
                    end
                    suum += suuum * int_garb.fxnmap[u][salc1.sh][abidx, salc1.ml]*int_garb.fxnmap[rv][salc2.sh][bbidx, salc2.ml]*int_garb.fxnmap[tw][salc3.sh][cbidx, salc3.ml]*int_garb.fxnmap[tsx][salc4.sh][dbidx, salc4.ml]
                    end
                end
            end
        end
    end
    Rb = int_garb.symtext.atom_map[salc2.atom,R]
    Tc = int_garb.symtext.atom_map[salc3.atom,T]
    TSd = int_garb.symtext.atom_map[salc4.atom,int_garb.symtext.mult_table[T,S]]

    suum *= (int_garb.irr_dims[salc1.irrep]*int_garb.irr_dims[salc2.irrep]*int_garb.irr_dims[salc3.irrep]*int_garb.irr_dims[salc4.irrep])/(λTs[Tidx]*(int_garb.g^4))
    end
    @timeit to "AO int calc" out = ERI_2e4c(bset, abars[salc1.atom][salc1.bfxn][1], abars[Rb][salc2.bfxn][1], abars[Tc][salc3.bfxn][1], abars[TSd][salc4.bfxn][1])[abidx,bbidx,cbidx,dbidx]
    eri[s1,s2,s3,s4] += suum*out*salc1.γ*salc2.γ*salc3.γ*salc4.γ
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
    #display()
    println(findmax(broadcast(abs, eri-SERI)))
    show(to)
    return 0
    #return eri
end
