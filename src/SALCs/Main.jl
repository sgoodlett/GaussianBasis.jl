
#struct SALC{F<:Real,I<:Integer}
#    coeffs::Vector{F}
#    irrep::String
#    atom::I
#    bfxn::I
#    sh::I # This is actually l+1, kinda weird but it's so we can index into things, I'll change this later, but I don't want to reset julia rn
#    ml::I # 1 <= ml <= 2l+1
#    i::I
#    r::I
#    γ::F
#end

mutable struct SALCblock
    irrep::String
    lcao::Array{Float64,2}
end

#for a basis function that is mapped onto a SEA (atom1 -> atom2)
#we need a way of finding the index of that basis function on atom 2

#the purpose of this function is to count the indicies between the range
#of indices spanned by these two atoms
function obstruct(atom1, atom2, bset)
    if abs(atom1 - atom2) == 1
        obstruction = 0
    else
        obstruction = 0
        slice = bset.basis_per_atom
        A = atom1 + 1
        B = atom2 - 1
        for x in slice[A:B]
            obstruction += x
        end 
    end

    return obstruction
end


#returns the max am value within the basis set object
#for the SH recursion function

function maxamcheck(bset)
    maxamval = 0
    for i in 1:bset.nshells
        am = bset.basis[i].l
        if am > maxamval
            maxamval = am
        end
    end
    return maxamval
end

#build the vector of SH rotation arrays 
function collectRotations(maxam,symels)
    rotated = Vector{Matrix{Float64}}[]
    #rotated = Vector{Matrix{Float64}}
    for i in 1:length(symels)
        push!(rotated, generateRotations(maxam, symels[i].rrep))
    end
    return rotated
end

#creates a vector of vectors for l values of the SH on atom centers
function basis_am(bset)
    molecule_am = []
    basis = 0
    for i in bset.shells_per_atom
        amperatom = []
        for j in 1:i
            push!(amperatom, bset.basis[j + basis].l)
        end
        basis += i
        push!(molecule_am, amperatom) 
    end
    #println("molcule_am $molecule_am")
    return molecule_am
end

#initialize the salc structure per irrep of pg

function salc_irreps(ct, bset)
    salcs = []
    for i in ct.irreps
        #println("Irrep $i")
        push!(salcs, SALCblock(i, Array{Float64}(undef,bset.nbas,0)))
    end
    return salcs
end

#function that adds to the salc struct if lcao is unique
function addlcao!(salcs, salc, ir, irrep, salc_vector, atomidx, bfxnidx, l, ml, gammas)
    #New = []
    rank_atol = 1e-10
    for r = 1:size(salc)[1]
    for (i,s) in enumerate(salc[:,r]) # i is for salc_vector, r = 1 always! WRONG STUPID
        check = true
        if gammas[i,r] == 0.0
            check = false
        end
        for y in eachcol(salcs[ir].lcao)
            if isapprox(s, y, atol = 1e-6) #|| isapprox(s, -y, atol = 1e-6)
                check = false
                break
            elseif length(salcs[ir].lcao) > 1 && rank(hcat(salcs[ir].lcao, s), atol=rank_atol) <= rank(salcs[ir].lcao, atol=rank_atol)
                check = false
                break
            end
        end
        if check
            #push!(salcs[ir].lcao, s)
            salcs[ir].lcao = hcat(salcs[ir].lcao, s)
            push!(salc_vector, SALC{Float64,Int64}(s, irrep, atomidx, bfxnidx, l+1, ml, i, r, gammas[i,r]))
        end
    end
    end 
    return salcs
end


#non-abelian projection operator for real-spherical harmonics
function ProjectionOp(mol, bset, symtext)
    salc_vector = Vector{SALC{Float64,Int64}}()
    D = Molecules.Symmetry.buildD(mol)
    SEAs = Molecules.Symmetry.findSEA(D, 5)
    maxam = maxamcheck(bset)
    rotated = collectRotations(maxam, symtext.symels)
    salcs = salc_irreps(symtext.ctab, bset)
    nbas_vec = bset.basis_per_atom
    outers = basis_am(bset)
    
    #loop over SEA sets
    basis = 1
    sea_chk = []
    for (atom_idx, atom) in enumerate(mol)
        for (sea_idx,sea) in enumerate(SEAs)
            if atom_idx in sea.set
                if sea_idx in sea_chk
                    basis += nbas_vec[atom_idx]
                else
                    push!(sea_chk, sea_idx)
                    equivatom = SEAs[sea_idx].set[1]
                    bsfxn_counter = 0
                    for (k, l) in enumerate(outers[equivatom])
                        for ml in 1:2*l + 1
                            bsfxn_counter += 1
                            #now we loop over irreps of the character Table
                            #grab the rotation slice, 1 x 2*l + 1 shape
                            #broadcast char * slice (previously called result)
                            for (ir, irrep) in enumerate(symtext.ctab.irreps)
                                irrmat = eval(Meta.parse("Molecules.Symmetry.CharacterTables.irrm_$(symtext.pg)"))[irrep]
                                dims = size(irrmat[1])
                                salc = [zeros(bset.nbas) for x in 1:dims[1], y in 1:dims[1]]
                                salc = projection(salc, bset, symtext, irrep, irrmat, rotated, l, ml, equivatom, basis, nbas_vec, sea)
                                chk = false
                                for brodx = 1:dims[1]
                                    if sum(broadcast(abs, salc[1, brodx])) > 1e-8
                                        chk = true
                                        break
                                    end
                                end
                                if chk #sum(broadcast(abs, salc[1])) > 1e-14
                                    gammas = zeros(Float64, dims)
                                    for (i, salcy) in enumerate(salc)
                                        #find largest element by magnitude
                                        #convention for these salcs: largest element in salc needs to be positive, apply to pfs
                                        index = findmax(broadcast(abs, salcy))
                                        if index[1] < 1e-7
                                            continue
                                        end
                                        factor = sum(x -> x^2, salcy)
                                        if abs(factor) < 1e-8
                                            println("Warning: Zero SALC detected")
                                        end
                                        if salcy[index[2]] < 0
                                            γ = -1/sqrt(factor) # This is for salc_vector
            			                    salc[i] = (-1.0/sqrt(factor))*salcy
            			                else
                                            γ = 1/sqrt(factor) # This is for salc_vector
                                            salc[i] = (1.0/sqrt(factor))*salcy
                                        end
                                        gammas[i] = γ
                                    end
                                    salcs = addlcao!(salcs, salc, ir, irrep, salc_vector, equivatom, bsfxn_counter, l, ml, gammas)
                                end 
                            end
                        end 
                        basis += l*2 +1
                    end
                end
            end
        end
    end
    #now reshape into ao by so array
    bigg = []
    so_irrep = []
    for (i, x) in enumerate(salcs)
	    sal = reshape([(salcs[i].lcao...)...], bset.nbas, :)
	    bigg = vcat(bigg, sal...)
	    so_irrep = vcat(so_irrep, fill(i, size(salcs[i].lcao)[1]))
	    salcs[i].lcao = sal
    end
    excess_salc_num = floor(Int64,size(bigg)[1]//bset.nbas)-bset.nbas
    bigg = reshape(bigg, bset.nbas, bset.nbas+excess_salc_num) 
    bigg = convert(Array{Float64}, bigg)
    #sort!(salc_vector, by = x->(Molecules.Symmetry.CharacterTables.irrep_sort_idx(x.irrep), x.atom, x.bfxn, x.i, x.r))
    sort!(salc_vector, by = x->(Molecules.Symmetry.CharacterTables.irrep_sort_idx(x.irrep)))
    if excess_salc_num == 0
        salc_chk = LinearAlgebra.det(bigg)
        if abs(salc_chk) < 1e-8
            display(bigg)
            throw(ErrorException("Error, the SALCs are linearly dependent!!!"))
        end
    else
        println("Allowing too many SALCs through: ($excess_salc_num), cannot check SALCs for linear dependence")
    end
    return salcs, bigg, so_irrep, salc_vector, rotated#, symtext
end

function garbo()
    molfn = "/home/smg13363/Molecules.jl/test/xyz/water.xyz"
    basisnm = "cc-pvdz"
    mol = Molecules.parse_file(molfn)
    bset = GaussianBasis.BasisSet(basisnm, mol)
    mol,symtext = Molecules.Symmetry.CharacterTables.symtext_from_mol(mol)    
    #mol = Molecules.translate(mol, Molecules.center_of_mass(mol))
    #D = Molecules.Symmetry.buildD(mol)
    #SEAs = Molecules.Symmetry.findSEA(D, 5)
    irrmat = eval(Meta.parse("Molecules.Symmetry.CharacterTables.irrm_$(symtext.pg)"))["A1"]
    dims = size(irrmat[1])
    salc = [zeros(bset.nbas) for x in 1:dims[1], y in 1:dims[1]]
    maxam = maxamcheck(bset)
    rotated = collectRotations(maxam, symtext.symels)
    display(rotated[2][3])
    salcs = salc_irreps(symtext.ctab)
    nbas_vec = bset.basis_per_atom
    #outers = basis_am(bset)
    #span = Any[]
    #irrep = 1
    #salc = []
    for i = 1:5
        projection(salc, bset, symtext, "A1", irrmat, rotated, 2, i, 1, 10, nbas_vec, 0)
        display(normalize(salc[1][10:14]))
    end
end

function projection(salc, bset, symtext, irrep, irrmat, rotated, l, ml, equivatom, basis, nbas_vec, sea)
    for (op, class) in enumerate(symtext.class_map)
        char = symtext.ctab[irrep][class]
        rotate_result = rotated[op]
        #if the am = 0, the operation has no effect on the spherical harmonic, hence result = 1
        #result on equivalent atom (itself or within SEA set)  
        if l + 1 == 1
            result = [1]
        else
            #result = rotate_result[l + 1][ml,:]
            result = rotate_result[l + 1][:,ml]
        end 
         
        #multiply the character by the rotation vector, and check what atom its centered on
        
        #NEED TO INCLUDE CLASS ORDER WHEN CLASSES ARE FIXED IN CHARTABLE
        atom2 = symtext.atom_map[:,op][equivatom]
        if atom2 == equivatom
            for i in 1:length(irrmat[op])
                salc[i][basis:basis + 2*l] += irrmat[op][i] .*result
            end

        else
            obstruction = obstruct(equivatom, atom2, bset)
            offset = obstruction + nbas_vec[equivatom]
            for i in 1:length(irrmat[op])
                salc[i][basis + offset:basis + offset + 2*l] += irrmat[op][i] .*result
            end
        end
    end
    return salc * (symtext.ctab.irrep_dims[irrep] / symtext.order)
end
