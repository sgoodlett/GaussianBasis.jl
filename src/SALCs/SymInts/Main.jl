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

function Soverlap(stalcs, symtext, irrm)
    l = length(stalcs)
    S = zeros(Float64, l, l)
    for s1 = 1:l
        for s2 = 1:l
            if stalcs[s1].irrep == stalcs[s2].irrep && stalcs[s1].i == stalcs[s2].i
                S[s1, s2] = 1.0
            else
                S[s1, s2] = 0.0
            end
        end
    end
    return S
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

function Λfxn(garbo::IntGarbage, salc, ab, p, R)
    out = 0.0
    for U in garbo.stat_subgrps[salc.atom]
        RU = garbo.symtext.mult_table[R, U]
        out += garbo.fxnmap[RU][salc.sh][ab, salc.ml]*garbo.irrm[salc.irrep][RU][p,salc.r]
    end
    return out * garbo.irr_dims[salc.irrep] / garbo.g
end

# Same function but R=1
function Λfxn(garbo::IntGarbage, salc, ab, p)
    out = 0.0
    for U in garbo.stat_subgrps[salc.atom]
        out += garbo.fxnmap[U][salc.sh][ab, salc.ml]*garbo.irrm[salc.irrep][U][p,salc.r]
    end
    return out * garbo.irr_dims[salc.irrep] / garbo.g
end

function Γfxn(garbo::IntGarbage, salc1, salc2, ab, bb, R, λr)
    out = 0.0
    for p = 1:garbo.irr_dims[salc1.irrep]
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
    aotoso, ignore_this, stalcs, doodad = GaussianBasis.SALCs.ProjectionOp(mol, bset)[2:5]
    return IntGarbage(symtext, irrm, doodad, subgrps, symtext.order, symtext.ctab.irrep_dims), stalcs, aotoso
end

function get_Stalcs(mol, bset)
    #rt2 = 1/sqrt(2)
    #rt3 = 1/sqrt(3)
    #rt6 = 1/sqrt(6)
    #n1s  = [1,0,0,0,0,0,0,0]
    #n2s  = [0,1,0,0,0,0,0,0]
    #npz  = [0,0,0,0,1,0,0,0]
    #nxy1 = [0,0,1,0,0,0,0,0]
    #nxy2 = [0,0,0,1,0,0,0,0]
    #ha1  = [0,0,0,0,0,rt3,rt3,rt3]
    #he1  = [0,0,0,0,0,2rt6,-rt6,-rt6]
    #he2  = [0,0,0,0,0,0,rt2,-rt2]
    #stalcs = [Stalc(n1s, "A1", 1, 1, 1, 1, 1.0),
    #          Stalc(n2s, "A1", 1, 2, 1, 1, 1.0),
    #          Stalc(npz, "A1", 1, 5, 1, 1, 1.0),
    #          Stalc(ha1, "A1", 2, 1, 1, 1, sqrt(3)), # (2/6 * rt3)^-1
    #          Stalc(nxy1, "E", 1, 3, 1, 1, 2.0),
    #          Stalc(nxy2, "E", 1, 3, 2, 1, 2.0),
    #          Stalc(he1,  "E", 2, 1, 1, 1, sqrt(6)),
    #          Stalc(he2,  "E", 2, 1, 2, 1, sqrt(6))]#(3/2)*sqrt(2)*(2/sqrt(3)))]
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

function sintegrals(molfn::String, basis, int_type::Function)
    mol = Molecules.parse_file(molfn)
    return sintegrals(mol, basis, int_type)
end

function sintegrals(mol::Vector{Molecules.Atom}, basis, int_type::Function)
    #mol = Molecules.parse_file(molfn)
    #bset = BasisSet(basis, mol)
    bset = BasisSet(basis, mol)
    symtext = Molecules.Symmetry.CharacterTables.symtext_from_mol(mol)
    display(symtext)
    abars = build_thang(bset)
    int_garb, stalcs, aotoso = get_garbage(mol, symtext, bset)
    display(aotoso)
    Sboy = int_type(bset)
    sym_S = transpose(aotoso) * Sboy * aotoso
    display(stalcs)
    porque = zeros(Float64, size(aotoso))
    for i = 1:length(stalcs)
        porque[:,i] = aotoso[:,i]-stalcs[i].coeffs
    end
    println(findmax(broadcast(abs, porque)))
    slength = length(stalcs)
    S = zeros(Float64, slength, slength)
    for s1 = 1:slength
        for s2 = s1:slength
            salc1 = stalcs[s1]
            salc2 = stalcs[s2]
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
                        for (abidx,ab) in enumerate(abars[salc1.atom][salc1.bfxn][2]) # TODO
                            for (bbidx,bb) in enumerate(abars[salc2.atom][salc2.bfxn][2])
                                Γ = Γfxn(int_garb, salc1, salc2, abidx, bbidx, R, λs[Ridx])
                                out = int_type(bset, abars[salc1.atom][salc1.bfxn][1], abars[Rb][salc2.bfxn][1])[abidx,bbidx] # shell_a, shell_b, ml_ab, ml_bb
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
    display(sym_S)
    println(findmax(broadcast(abs, S-sym_S)))
    return S
end