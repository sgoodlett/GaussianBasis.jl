using Molecules


#indices 

mutable struct SALC
    irrep
    lcao
end


include("SHRotations.jl")


c3 = cos(π/3)
s3 = sin(π/3)
c5 = cos(π/5)
s5 = sin(π/5)
c25 = cos(2π/5)
s25 = sin(2π/5)
irrm_C2v = Dict(
    "A1" => [[1],[1],[1],[1]],
    "A2" => [[1],[1],[-1],[-1]],
    "B1" => [[1],[-1],[1],[-1]],
    "B2" => [[1],[-1],[-1],[1]])

irrm_C3v = Dict(
    "A1" => [[1],[1],[1],[1],[1],[1]],
    "A2" => [[1],[1],[1],[-1],[-1],[-1]],
    "E"  => [[1 0;0 1], [-c3 -s3; s3 -c3], [-c3 s3; -s3 -c3], [1 0;0 -1], [-c3 -s3; -s3 c3], [-c3 s3; s3 c3]])
    #"E"  => [[1 0;0 1], [-c3 -s3; s3 -c3], [-c3 s3; -s3 -c3], [-c3 -s3; -s3 c3], [1 0;0 -1],  [-c3 s3; s3 c3]])

    irrm_Td  = Dict(
        "A1" => [ [1], # E
                  [1], [1], [1], [1], [1], [1], [1], [1], # C3
                  [1], [1], [1], # C2
                  [1], [1], [1], [1], [1], [1], # σd
                  [1], [1], [1], [1], [1], [1]], # S4
        "A2" => [ [1], 
                  [1], [1], [1], [1], [1], [1], [1], [1], 
                  [1], [1], [1], 
                  [-1],[-1],[-1],[-1],[-1],[-1],
                  [-1],[-1],[-1],[-1],[-1],[-1]],
        "E"  => [[1 0;0 1],
                 [-c3 -s3; s3 -c3],[-c3 s3; -s3 -c3],[-c3 -s3; s3 -c3],[-c3 s3; -s3 -c3],[-c3 -s3; s3 -c3],[-c3 s3; -s3 -c3],[-c3 -s3; s3 -c3],[-c3 s3; -s3 -c3],
                 [1 0;0 1],[1 0;0 1],[1 0;0 1],
                 [1 0; 0 -1],[1 0; 0 -1],[-c3 s3; s3 c3],[-c3 s3; s3 c3],[-c3 -s3; -s3 c3],[-c3 -s3; -s3 c3],
                 [-c3 -s3; -s3 c3],[-c3 -s3; -s3 c3],[-c3 s3; s3 c3],[-c3 s3; s3 c3],[1 0;0 -1],[1 0;0 -1]],
        "T1" => [[1 0 0;0 1 0;0 0 1], # E
                 [0 0 1;1 0 0;0 1 0],[0 1 0;0 0 1;1 0 0],[0 0 1;-1 0 0;0 -1 0],[0 -1 0;0 0 -1;1 0 0], # C3 (α,β)
                 [0 0 -1;1 0 0;0 -1 0],[0 1 0;0 0 -1;-1 0 0],[0 0 -1;-1 0 0;0 1 0],[0 -1 0;0 0 1;-1 0 0], # C3 (γ,δ)
                 [1 0 0;0 -1 0;0 0 -1],[-1 0 0;0 1 0;0 0 -1],[-1 0 0;0 -1 0;0 0 1], # C2 (x,y,z)
                 [0 1 0;1 0 0;0 0 -1],[0 -1 0;-1 0 0;0 0 -1],[0 0 1;0 -1 0;1 0 0],[0 0 -1;0 -1 0;-1 0 0],[-1 0 0;0 0 1;0 1 0],[-1 0 0;0 0 -1;0 -1 0], # σd (xy,xz,yz)
                 [1 0 0;0 0 1;0 -1 0],[1 0 0;0 0 -1;0 1 0],[0 0 -1;0 1 0;1 0 0],[0 0 1;0 1 0;-1 0 0],[0 1 0;-1 0 0;0 0 1],[0 -1 0;1 0 0;0 0 1]], # S4 (x,y,z)
        "T2" => [[1 0 0;0 1 0;0 0 1], # E
                 [0 0 1;1 0 0;0 1 0],[0 1 0;0 0 1;1 0 0],[0 0 1;-1 0 0;0 -1 0],[0 -1 0;0 0 -1;1 0 0], # C3 (α,β)
                 [0 0 -1;1 0 0;0 -1 0],[0 1 0;0 0 -1;-1 0 0],[0 0 -1;-1 0 0;0 1 0],[0 -1 0;0 0 1;-1 0 0], # C3 (γ,δ)
                 [1 0 0;0 -1 0;0 0 -1],[-1 0 0;0 1 0;0 0 -1],[-1 0 0;0 -1 0;0 0 1], # C2 (x,y,z)
                 [0 -1 0;-1 0 0;0 0 1],[0 1 0;1 0 0;0 0 1],[0 0 -1;0 1 0;-1 0 0],[0 0 1;0 1 0;1 0 0],[1 0 0;0 0 -1;0 -1 0],[1 0 0;0 0 1;0 1 0], # σd (xy,xz,yz)
                 [-1 0 0;0 0 -1;0 1 0],[-1 0 0;0 0 1;0 -1 0],[0 0 1;0 -1 0;-1 0 0],[0 0 -1;0 -1 0;1 0 0],[0 -1 0;1 0 0;0 0 -1],[0 1 0;-1 0 0;0 0 -1]]) # S4 (x,y,z)
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
    rotated = []
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
    println("molcule_am $molecule_am")
    return molecule_am
end

#initialize the salc structure per irrep of pg

function salc_irreps(ct)
    salcs = []
    for i in ct.irreps
        println("Irrep $i")
        push!(salcs, SALC(i, []))
    end
    return salcs
end

#function that adds to the salc struct if lcao is unique

function addlcao(salcs, salc, ir, irrep)
    New = []
    for s in salc[:,1]
        check = true
        for y in salcs[ir].lcao
            if isapprox(s, y, atol = 1e-8)
                check = false
                break
            end
        end
        if check
	    #println("s $s")
	    #t = [(s...)...]
	    #println("t $t")
            push!(salcs[ir].lcao, s)
            push!(New, s)
        end
    end
        
    return salcs
end


#non-abelian projection operator for real-spherical harmonics

function ProjectionOp(mol, bset)
   
    symtext = Molecules.Symmetry.CharacterTables.symtext_from_mol(mol)    
    #functions below will be in symtext## 
    println("symtext $symtext")    
    mol = Molecules.translate(mol, Molecules.center_of_mass(mol))
    #pg = symtext.pg
    ##symels = symtext.symels
    #ct = symtext.ctab
    #class_map = symtext.class_map
    #amap = symtext.atom_map
    D = Molecules.Symmetry.buildD(mol)
    SEAs = Molecules.Symmetry.findSEA(D, 5)

    #pg, paxis, saxis = Molecules.Symmetry.find_point_group(mol)
    #symels = Molecules.Symmetry.CharacterTables.pg_to_symels(pg)
    #symels = Molecules.Symmetry.CharacterTables.rotate_symels_to_mol(symels, paxis, saxis)
    #ct = Molecules.Symmetry.CharacterTables.pg_to_chartab(pg)
    #class_map = Molecules.Symmetry.CharacterTables.generate_symel_to_class_map(symels, ct)
    #D = Molecules.Symmetry.buildD(mol)
    #SEAs = Molecules.Symmetry.findSEA(D, 5)
    #amap = Molecules.Symmetry.CharacterTables.get_atom_mapping(mol, symels)
    #functions above will be in symtext## 
    
    maxam = maxamcheck(bset)
    rotated = collectRotations(maxam, symtext.symels)
    salcs = salc_irreps(symtext.ctab)
    nbas_vec = bset.basis_per_atom
    outers = basis_am(bset)
    span = Any[]
    
    #loop over SEA sets
    basis = 1
    for sea in 1:length(SEAs)
        equivatom = SEAs[sea].set[1]
        #loop over basis function on center on each equivatom
        for (k, l) in enumerate(outers[equivatom])
            for ml in 1:2*l + 1
                #now we loop over irreps of the character Table

                #grab the rotation slice, 1 x 2*l + 1 shape
                #broadcast char * slice (previously called result)
                for (ir, irrep) in enumerate(symtext.ctab.irreps)
	
                    #until the irrmats become part of the symtext, they will be specified manually
		   
		    irrmat = irrm_C2v[irrep]
                    #irrmat = irrm_C3v[irrep]
                    #irrmat = irrm_Td[irrep]
                    
                    #dims = size(irrm_C3v[irrep][1])
                    dims = size(irrm_C2v[irrep][1])
                    
                    salc = [zeros(bset.nbas) for x in 1:dims[1], y in 1:dims[1]]
                    #salc = zeros(bset.nbas) for x in 1:dims[1], y in 1:dims[1]
                    salc = projection(salc, bset, symtext.class_map, ir, symtext.ctab, irrmat, rotated, l, ml, symtext.atom_map, equivatom, basis, nbas_vec, sea)
                    if sum(broadcast(abs, salc[1])) > 1e-14
                        #println("irrep, salc*** $irrep $salc")
                        for (i, salcy) in enumerate(salc)
                            #find largest element by magnitude
                            #convention for these salcs: largest element in salc needs to be positive, apply to pfs
                   
                            index = findmax(broadcast(abs, salcy))
                            factor = sum(x -> x^2, salcy)
                            if salcy[index[2]] < 0
			                    salc[i] = (-1.0/sqrt(factor))*salcy
			                else
                                salc[i] = (1.0/sqrt(factor))*salcy
                            end

                        end
                        salcs = addlcao(salcs, salc, ir, irrep)
			#end
                    end 
                   
                end
            end 
            basis += l*2 +1
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
    bigg = reshape(bigg, bset.nbas, bset.nbas) 
    bigg = convert(Array{Float64}, bigg)
    return salcs, bigg, so_irrep, symtext
end

function projection(salc, bset, class_map, ir, ct, irrmat, rotated, l, ml, amap, equivatom, basis, nbas_vec, sea)
    for (op, class) in enumerate(class_map)
        char = ct.characters[ir,:][class]
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
        atom2 = amap[:,op][equivatom]
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
    return salc
end
