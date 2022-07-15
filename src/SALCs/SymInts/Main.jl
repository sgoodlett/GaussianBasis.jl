# Shell is probs a bad name, actually bfxn idx
struct Stalc
    coeffs
    irrep
    atom
    sh
    i
    j
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

function c3v_test()
    c3v_irr = Molecules.Symmetry.CharacterTables.irrm_C3v
    bs = "sto-3g"
    mol = Molecules.parse_file("/home/smg13363/Molecules.jl/test/xyz/ammonia.xyz")
    bset = BasisSet(bs, mol)
    #println(bset)
    symtext = Molecules.Symmetry.CharacterTables.symtext_from_mol(mol)
    r32 = sqrt(3)/2
    fxn_sets = [[1],[2],[3,4],[5]] # By shell
    set_lens = [1,1,2,2,1] # By shell
    abars = [[[1],[2],[3,4],[3,4],[5]],[[6]],[[7]],[[8]]] # Atom, fxn
    bsfxn_maps = [
        # Atom, fxn, GroupOperation, fxn vector    
        # N    
        [[[1,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0]],
         [[0,1,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,1,0,0,0,0,0,0]],
         [[0,0,1,0,0,0,0,0],[0,0,-0.5,r32,0,0,0,0],[0,0,-0.5,-r32,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,-0.5,-r32,0,0,0,0],[0,0,-0.5,r32,0,0,0,0]],
         [[0,0,1,0,0,0,0,0],[0,0,-0.5,r32,0,0,0,0],[0,0,-0.5,-r32,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,-0.5,-r32,0,0,0,0],[0,0,-0.5,r32,0,0,0,0]],
         [[0,0,0,0,1,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,1,0,0,0]]],
        # H1
        [[[0,0,0,0,0,1,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,1,0,0]]]
    ]
    rt6 = 1/sqrt(6)
    rt3 = 1/sqrt(3)
    rt2 = 1/sqrt(2)
    n1s  = [1,0,0,0,0,0,0,0]
    n2s  = [0,1,0,0,0,0,0,0]
    npz  = [0,0,0,0,1,0,0,0]
    nxy1 = [0,0,1,0,0,0,0,0]
    nxy2 = [0,0,0,1,0,0,0,0]
    ha1  = [0,0,0,0,0,rt3,rt3,rt3]
    he1  = [0,0,0,0,0,2rt6,-rt6,-rt6]
    he2  = [0,0,0,0,0,0,rt2,-rt2]
    stalcs = [Stalc(n1s, "A1", 1, 1, 1, 1),
              Stalc(n2s, "A1", 1, 2, 1, 1),
              Stalc(npz, "A1", 1, 5, 1, 1),
              Stalc(ha1, "A1", 2, 1, 1, 1),
              Stalc(nxy1, "E", 1, 3, 1, 1),
              Stalc(nxy2, "E", 1, 3, 2, 1),
              Stalc(he1,  "E", 2, 1, 1, 1),
              Stalc(he2,  "E", 2, 1, 2, 1)]
    l = length(stalcs)
    S = zeros(Float64, l, l)
    g = 6
    ns = Dict("A1"=>1, "A2"=>1, "E"=>2)
    if true
    for s1 = 1:l
        for s2 = s1:l
            α = stalcs[s1].irrep
            β = stalcs[s2].irrep
            A = stalcs[s1].atom
            B = stalcs[s2].atom
            a = stalcs[s1].sh
            b = stalcs[s2].sh
            i = stalcs[s1].i
            j = stalcs[s2].i
            r = stalcs[s1].j
            s = stalcs[s2].j
            if α == β && i == j
                U = get_atom_subgroup(A, symtext)
                V = get_atom_subgroup(B, symtext)
                Rs, λs = get_Rs(U, V, symtext)
                println("SALCs $(s1) and $(s2). Rs:$(Rs) λs:$(λs)")
                for (Ridx,R) in enumerate(Rs)
                    for ab in abars[A][a]
                        for bb in abars[B][b]
                            Γ = 0.0
                            for p = 1:ns[α]
                                uΛ = 0.0
                                for u in U
                                    println("\t A:$(A), a:$(a), α:$(α), u:$(u), ab:$(ab), p:$(p), r:$(r)")
                                    println("\t\t$(bsfxn_maps[A][a][u][ab])")
                                    println("\t\t$(c3v_irr[α][u][p][r])")
                                    uΛ += bsfxn_maps[A][a][u][ab]*c3v_irr[α][u][p][r]
                                end
                                uΛ *= (ns[α]/g)
                                vΛ = 0.0
                                for v in V
                                    rv = symtext.mult_table[R,v]
                                    println("\t B:$(B), b:$(b), rv:$(rv), bb:$(bb)")
                                    println("\t\t$(bsfxn_maps[B][b][rv][bb])")
                                    println("\t\t$(c3v_irr[α][rv][p][s])")
                                    vΛ += bsfxn_maps[B][b][rv][bb]*c3v_irr[α][rv][p][s]
                                end
                                vΛ *= (ns[α]/g)
                                println("\t Final values u:$(uΛ), v:$(vΛ)")
                                Γ += uΛ * vΛ
                                println("\t Γ:$(Γ)")
                            end
                            println("\t Γ:$(Γ), g:$(g), nα:$(ns[α]), λR:$(λs[Ridx])")
                            Γ *= g / (ns[α]*λs[Ridx])
                            println("Double check: Γ=$(Γ) and S[$(s1),$(s2)]=$(S[s1,s2])")
                            s = 1.0
                            S[s1, s2] += Γ*s
                            if s1 != s2
                                S[s2, s1] += Γ
                            end
                        end
                    end
                end
            end
        end
    end
    return S
    end
    U = get_atom_subgroup(1, symtext)
    V = get_atom_subgroup(1, symtext)
    println(U)
    println(V)
    Rs, λs = get_Rs(U, V, symtext)
end