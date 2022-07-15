#using Molecules
#using RowEchelon
struct SALC
    irrep
    lcao
end


include("SHRotations.jl")

#for a basis function that is mapped onto a SEA (atom1 -> atom2)
#we need a way of finding the index of that basis function on atom 2

#the purpose of this function is to count the indicies between the range
#of indices spanned by these two atoms
function obstruct(atom1, atom2, bset)
    if abs(atom1 - atom2) == 1
        #println("no chance")
        obstruction = 0
    else
        obstruction = 0
        slice = init_ao_to_so(bset)
        #println("slice")
        A = atom1 + 1
        B = atom2 - 1
        for x in slice[A:B]
            obstruction += x
        end 
    end

    return obstruction
end


#returns the max am value within the basis set object
#for the SH recursion relation

function amcheck(bset)
    maxamval = 0 
    for i in 1:length(bset.basis)
        for j in 1:length(bset.basis[i])
            am = bset.basis[i][j].l
            if am > maxamval
                maxamval = am
            end
        end
    end
    return maxamval
end

#creates a vector of vectors for l values of the SH on atom centers
function basis_am(bset)
    outer = []
    for i in 1:length(bset.atoms)
        inner = []
        for j in 1:length(bset.basis[i])
            shellam = bset.basis[i][j].l
            push!(inner, shellam)
        end
        push!(outer, inner)
    end
    return outer
end

#finds sum over each basis functions 2*l + 1 per atom
function init_ao_to_so(bset)
    outer = basis_am(bset)
    #println("outer $outer")
    count = []
    counts = []
    for i in 1:length(outer)
        counter = 0
        for j in outer[i]
            push!(counts, 2*j + 1)
            counter += 2*j + 1
        end
        push!(count, counter)
    end
    #println("Counts $counts count $count")
    return count
end

#build the vector of SH rotation arrays 
function collectRotations(maxam,symels)
    rotated = []
    for i in 1:length(symels)
        push!(rotated, generateRotations(maxam, symels[i].rrep))
    end
    return rotated
end

function salc_irreps(ct)
    salcs = []
    for i in ct.irreps
        println("Irrep $i")
        push!(salcs, SALC(i, Any[]))
    end
    return salcs
end

#non-abelian projection operator for real-spherical harmonics

function ProjectionOp(mol, bset)
    mol = Molecules.translate(mol, Molecules.center_of_mass(mol))
    pg, paxis, saxis = Molecules.Symmetry.find_point_group(mol)
    symels = Molecules.Symmetry.CharacterTables.pg_to_symels(pg)
    symels = Molecules.Symmetry.CharacterTables.rotate_symels_to_mol(symels, paxis, saxis)
    ct = Molecules.Symmetry.CharacterTables.pg_to_chartab(pg)
    class_map = Molecules.Symmetry.CharacterTables.generate_symel_to_class_map(symels, ct)
    println("mol $mol symels $symels ct $ct")
    maxam = amcheck(bset)
    rotated = collectRotations(maxam, symels)
    salcs = salc_irreps(ct)
    println(salcs)
    println(salcs[1].lcao) 
    nbas_vec = init_ao_to_so(bset)
    outers = basis_am(bset)
    println("These are the am values $outers")
    D = Molecules.Symmetry.buildD(mol)
    SEAs = Molecules.Symmetry.findSEA(D, 5)
    println("SEA Set ", SEAs)
    amap = Molecules.Symmetry.CharacterTables.get_atom_mapping(mol, symels)
    #amap = Molecules.atom_map(mol, d)
    #println("atommap $amap equivatom $equivatom")
    println("atom map $amap")
    println("nbas_vec $nbas_vec")
    span = Any[]
    #loop over SEA sets
    counters = Any[]
    basis = 1
    for sea in 1:length(SEAs)
        whytho = 0
        equivatom = SEAs[sea].set[1]
        #println("equivalent atom ", equivatom)
        #loop over basis function on center on each equivatom
        counter = 0
        #prelude is the number of salc positions iterated over before a RSH maps to another atom
        prelude = 1
        for (k, l) in enumerate(outers[equivatom])
            println("k $k, l $l")
            whytho += 1
            #println("This is the value ", l)
            for ml in 1:2*l + 1
                #println("This is projection ", ml)
                #now we loop over irreps of the character Table
                counter += 1

                #grab the rotation slice, 1 x 2*l + 1 shape
                #broadcast char * slice (previously called result)
                #for (op_index, class) in enumerate(class_map)
                #    println("operation $(symel[op_index]) class $class")
                #end
                for (ir, irrep) in enumerate(ct.irreps)
                    if ct.characters[ir,:][1] > 1
                        println(ct.characters[ir,:][1])
                        println("Partner function needed")
                        salc = zeros(bset.nbas)
                        salc = nondegenerate(salc, class_map, ir, ct, rotated, l, ml, amap, equivatom, basis, nbas_vec, sea)
                        #for (fop, fclass) in enumerate(class_map)
                        #   salc = zeros(bset.nbas)
                        #   psalc = zeros(bset.nbas)
                        #   salc, psalc = degenerate(salc, class_map, ir, ct, rotated, l, ml, amap, equivatom, basis, nbas_vec, psalc, fop)
                        #end
                    else
                        salc = zeros(bset.nbas)
                        salc = nondegenerate(salc, class_map, ir, ct, rotated, l, ml, amap, equivatom, basis, nbas_vec, sea)
                    end

                    if sum(broadcast(abs, salc)) > 1e-14
                        println("irrep, salc*** $irrep $salc")
                        
                        #append each lcao salc to the class object
                        
                        for i in 1:length(salcs)
                            if salcs[i].irrep == irrep
                                #println(irrep, " ", salc)
                                push!(salcs[i].lcao, salc)
                                #println("Salcs ", salcs[i].lcao)
                                if length(salcs[i].lcao) == 0
                                    println("right hurr")
                                    #push!(salcs[i].lcao, transpose(salc))
                                    push!(salcs[i].lcao, salc)
                                else
                                    println("error")
                                    #hcat(salcs[i].lcao, transpose(salc))
                                    #hcat(salcs[i].lcao, salc)
                                    vcat(salcs[i].lcao, salc)
                                end
                            end
                        end
                        
                        #println("basis val $basis")
                        #println("prelude val $prelude")
                        #if ct.characters[ir,:][1] > 1
                        #    println("Degenerate irrep, Partner Function needed")
                        #    span = unique(span)
                        #    println("span $span")
                        #    #osalc = pf(salc, sea, SEAs, k, l, basis, rotated, ct, outers, class_map, amap, nbas_vec)
                        #end 
                    end 
                end
            
            end 
            basis += l*2 +1
            prelude += l*2 +1
        end
    end
    println("""
    
    """)
    println("SALCS NOW $salcs")

    #println(salc)


    return salcs
end

function nondegenerate(salc, class_map, ir, ct, rotated, l, ml, amap, equivatom, basis, nbas_vec, sea)
    for (op, class) in enumerate(class_map)
        char = ct.characters[ir,:][class]
        println("char op class $char $op $class")
        rotate_result = rotated[op]
        
        #if the am = 0, the operation has no effect on the spherical harmonic, hence result = 1
        #result on equivalent atom (itself or within SEA set)  
        if l + 1 == 1
            result = [1]
        else
        result = rotate_result[l + 1][ml,:]
        end 
        
        #multiply the character by the rotation vector, and check what atom its centered on
        
        #NEED TO INCLUDE CLASS ORDER WHEN CLASSES ARE FIXED IN CHARTABLE
        atom2 = amap[:,op][equivatom]
        if atom2 == equivatom
            salc[basis:basis + 2*l] += char  .*result
        else
            obstruction = obstruct(equivatom, atom2, bset)
            offset = obstruction + nbas_vec[equivatom]
            salc[basis + offset:basis + offset + 2*l] += char .*result
        end
    end
    #THIS CODE NEEDS COMMENTED OUT UNTIL IRREDUCIBLE REPS ARE CREATED
    #Current method uses the salcs generated with trace projection operator to rotate the salc, check if it is degenerate, then project out w/gram-shmidt
    #check if the salc requires partner functions
    #if sum(broadcast(abs, salc)) > 1e-14 && ct.characters[ir,:][1] > 1
    #    for (pop, pclass) in enumerate(class_map)
    #        psalc = zeros(bset.nbas)
    #        rotate_result = rotated[pop]
    #        #if the am = 0, the operation has no effect on the spherical harmonic, hence result = 1
    #        #result on equivalent atom (itself or within SEA set)  
    #        if l + 1 == 1
    #            result = [1]
    #        else
    #        result = rotate_result[l + 1][ml,:]
    #        end
    #        #generate the partner function(s) 
    #        psalc = pfchang(result, pop, salc, class_map, ir, ct, rotated, l, ml, amap, equivatom, basis, nbas_vec, sea)
    #        if sum(broadcast(abs,(salc - psalc))) < 0.00001 || sum(broadcast(abs,(salc - -1*psalc))) < 0.00001
    #            println("keep looking for pf")
    #        else
    #            println("Need to orthogonalize")
    #            println(typeof(salc))
    #            println(typeof(psalc))
    #            salcpair = vcat(salc', psalc')
    #            #push!(salc, psalc)
    #            println("Salc $salcpair")
    #            salcpair = rref(salcpair)
    #            println("Salc $salcpair")
    #            #println(typeof(salcpair[2,:])) 
    #            salcpair[1,:] = salc'
    #            println("salcpair $salcpair") 
    #            break
    #        end
    #    end
    #end
    return salc
end

function pfchang(result,pop, salc, class_map, ir, ct, rotated, l, ml, amap, equivatom, basis, nbas_vec, sea)
    psalc = zeros(bset.nbas)
    for (x, y) in enumerate(SEAs[sea].set)
        #println("$x $y")
        atom3 = amap[:,pop][y]
        if length(SEAs[sea].set) == 1
            println("It has to map onto itself, just check the a.m. rotation")
            #println("$([basis:basis + 2*l])")
            #println("salc $(salc[basis:basis + 2*l])")
            coef = salc[basis:basis + 2*l]
            #coef .* result
            #take the vector vector product between coef and am rotation, grab the vector corresponding to ml to see how ml rotated
            tryit = coef * result'
            #psalc[basis:basis + 2*l] = salc[basis: basis + 2*l] .* result
            #psalc[basis:basis + 2*l] = tryit[ml,:]
            psalc[basis:basis + 2*l] = (salc[basis:basis + 2*l] * result')[ml,:]
        else
            println("$y mapped to $atom3")
            if equivatom == y == atom3
                println("CONDITION 0")
                #println("$equivatom $y $atom3")
                println("$([basis:basis + 2*l]) -> $([basis:basis + 2*l])")
                #psalc[basis:basis + 2*l] = salc[basis:basis + 2*l] #*result'
                psalc[basis:basis + 2*l] = (salc[basis:basis + 2*l] *result')[ml,:]
            #elseif atom3 == y
            elseif equivatom == y
                println("CONDITION 1")
                println("$y -> $atom3")
                #if equivatom == atom3
                #    println("DA SAME $equivatom $y $atom3")
                #else    
                #println("find the obstruction between this and the first atom of the list in [1, 3, 4] etc")
                obstruction = obstruct(equivatom,atom3, bset)
                offset = obstruction + nbas_vec[equivatom]
                #println("offset $offset")
                #println("$equivatom, but $y -> $atom3")
                #println("$([basis + offset:basis + offset + 2*l]) -> $([basis + offset:basis + offset + 2*l])")
                println("$([basis:basis + 2*l]) -> $([basis + offset:basis + offset + 2*l])")
                #psalc[basis + offset:basis + offset + 2*l] = salc[basis:basis + 2*l]
                psalc[basis + offset:basis + offset + 2*l] = (salc[basis:basis + 2*l] * result')[ml,:]
                #end
            elseif atom3 < y
                println("CONDITION 2")
                println("$y -> $atom3")
                #first, generate the offset from equivatom to y
                obstruction1 = obstruct(equivatom, y, bset)
                offset1 = obstruction1 + nbas_vec[equivatom]
                #println("obstruction1 offset1 $obstruction1 $offset1") 
                obstruction2 = - obstruct(atom3, y, bset)
                #offset2 = obstruction2 - nbas_vec[equivatom]
                #doesn't matter what you index nbas_vec with if all SEAs have equivalent basis functions
                offset2 = obstruction2 - nbas_vec[equivatom]
                #println("obstruction2 offset2 $obstruction2 $offset2") 
                println("$([basis + offset1:basis + offset1 + 2*l]) -> $([basis + offset1 + offset2:basis + offset1 + offset2 + 2*l])")
                #psalc[basis + offset1 + offset2:basis + offset1 + offset2 + 2*l] = salc[basis + offset1:basis + offset1 + 2*l] 
                psalc[basis + offset1 + offset2:basis + offset1 + offset2 + 2*l] = (salc[basis + offset1:basis + offset1 + 2*l] * result')[ml,:] 
                
                #println("$([basis:basis + 2*l]) -> $([basis + offset:basis + offset + 2*l])")
            elseif y == atom3
                println("CONDITION 3")
                println("$y -> $atom3")
                obstruction = obstruct(equivatom, atom3, bset)
                offset = obstruction + nbas_vec[equivatom]
                #println("offset $offset")
                #println("$([basis:basis + 2*l]) -> $([basis + offset:basis + offset + 2*l])")
                println("$([basis + offset:basis + offset + 2*l]) -> $([basis + offset:basis + offset + 2*l])")
                #psalc[basis + offset:basis + offset + 2*l] = salc[basis:basis + 2*l]
                #psalc[basis + offset:basis + offset + 2*l] = salc[basis + offset:basis + offset + 2*l]
                psalc[basis + offset:basis + offset + 2*l] = (salc[basis + offset:basis + offset + 2*l] * result')[ml,:]
            else
                println("CONDITION 4")
                println("$y -> $atom3")
                obstruction1 = obstruct(equivatom, y, bset)
                offset1 = obstruction1 + nbas_vec[equivatom]
                obstruction2 = obstruct(y, atom3, bset)
                offset2 = obstruction2 + nbas_vec[equivatom]
                println("$([basis + offset1:basis + offset1 + 2*l]) -> $([basis + offset1 + offset2:basis + offset1 + offset2 + 2*l])")
                #psalc[basis + offset1 + offset2:basis + offset1 + offset2 + 2*l] = salc[basis + offset1:basis + offset1 + 2*l] 
                psalc[basis + offset1 + offset2:basis + offset1 + offset2 + 2*l] = (salc[basis + offset1:basis + offset1 + 2*l] * result')[ml,:] 
            end 
        end
    end
    return psalc
end 

#function nondegenerate(salc, class_map, ir, ct, rotated, l, ml, amap, equivatom, basis, nbas_vec, sea)
#    for (op, class) in enumerate(class_map)
#        char = ct.characters[ir,:][class]
#        println("char op class $char $op $class")
#        rotate_result = rotated[op]
#        
#        #if the am = 0, the operation has no effect on the spherical harmonic, hence result = 1
#        #result on equivalent atom (itself or within SEA set)  
#        if l + 1 == 1
#            result = [1]
#        else
#        result = rotate_result[l + 1][ml,:]
#        end 
#        
#        #multiply the character by the rotation vector, and check what atom its centered on
#        
#        #NEED TO INCLUDE CLASS ORDER WHEN CLASSES ARE FIXED IN CHARTABLE
#        atom2 = amap[:,op][equivatom]
#        if atom2 == equivatom
#            salc[basis:basis + 2*l] += char  .*result
#        else
#            obstruction = obstruct(equivatom, atom2, bset)
#            offset = obstruction + nbas_vec[equivatom]
#            salc[basis + offset:basis + offset + 2*l] += char .*result
#        end
#    end
#    if sum(broadcast(abs, salc)) > 1e-14 && ct.characters[ir,:][1] > 1
#        println("This salc is degenerate and non-zero, pf needed $ir $salc")
#        println("SEAs[sea].set")
#        println(SEAs[sea].set)
#        println("this is where basis ended $basis")
#        println("this is where the am is at $l $ml")
#        for (pop, pclass) in enumerate(class_map)
#            psalc = zeros(bset.nbas)
#            rotate_result = rotated[pop]
#            #println("Rotate_result $rotate_result") 
#            #if the am = 0, the operation has no effect on the spherical harmonic, hence result = 1
#            #result on equivalent atom (itself or within SEA set)  
#            if l + 1 == 1
#                result = [1]
#            else
#            result = rotate_result[l + 1][ml,:]
#            end 
#            #println("RESULT $result")
#            for (x, y) in enumerate(SEAs[sea].set)
#                #println("$x $y")
#                atom3 = amap[:,pop][y]
#                if length(SEAs[sea].set) == 1
#                    println("It has to map onto itself, just check the a.m. rotation")
#                    #println("$([basis:basis + 2*l])")
#                    #println("salc $(salc[basis:basis + 2*l])")
#                    coef = salc[basis:basis + 2*l]
#                    #coef .* result
#                    #take the vector vector product between coef and am rotation, grab the vector corresponding to ml to see how ml rotated
#                    tryit = coef * result'
#                    #psalc[basis:basis + 2*l] = salc[basis: basis + 2*l] .* result
#                    psalc[basis:basis + 2*l] = tryit[ml,:]
#                    #psalc[basis:basis + 2*l] = [ salc[basis:basis + 2*l] * result'][ml,:]
#                else
#                    println("$y mapped to $atom3")
#                    if equivatom == y == atom3
#                        println("CONDITION 0")
#                        #println("$equivatom $y $atom3")
#                        println("$([basis:basis + 2*l]) -> $([basis:basis + 2*l])")
#                        psalc[basis:basis + 2*l] = salc[basis:basis + 2*l]
#                    #elseif atom3 == y
#                    elseif equivatom == y
#                        println("CONDITION 1")
#                        println("$y -> $atom3")
#                        #if equivatom == atom3
#                        #    println("DA SAME $equivatom $y $atom3")
#                        #else    
#                        #println("find the obstruction between this and the first atom of the list in [1, 3, 4] etc")
#                        obstruction = obstruct(equivatom,atom3, bset)
#                        offset = obstruction + nbas_vec[equivatom]
#                        #println("offset $offset")
#                        #println("$equivatom, but $y -> $atom3")
#                        #println("$([basis + offset:basis + offset + 2*l]) -> $([basis + offset:basis + offset + 2*l])")
#                        println("$([basis:basis + 2*l]) -> $([basis + offset:basis + offset + 2*l])")
#                        psalc[basis + offset:basis + offset + 2*l] = salc[basis:basis + 2*l]
#                        #end
#                    elseif atom3 < y
#                        println("CONDITION 2")
#                        println("$y -> $atom3")
#                        #first, generate the offset from equivatom to y
#                        obstruction1 = obstruct(equivatom, y, bset)
#                        offset1 = obstruction1 + nbas_vec[equivatom]
#                        #println("obstruction1 offset1 $obstruction1 $offset1") 
#                        obstruction2 = - obstruct(atom3, y, bset)
#                        #offset2 = obstruction2 - nbas_vec[equivatom]
#                        #doesn't matter what you index nbas_vec with if all SEAs have equivalent basis functions
#                        offset2 = obstruction2 - nbas_vec[equivatom]
#                        #println("obstruction2 offset2 $obstruction2 $offset2") 
#                        println("$([basis + offset1:basis + offset1 + 2*l]) -> $([basis + offset1 + offset2:basis + offset1 + offset2 + 2*l])")
#                        psalc[basis + offset1 + offset2:basis + offset1 + offset2 + 2*l] = salc[basis + offset1:basis + offset1 + 2*l] 
#                        
#                        #println("$([basis:basis + 2*l]) -> $([basis + offset:basis + offset + 2*l])")
#                    elseif y == atom3
#                        println("CONDITION 3")
#                        println("$y -> $atom3")
#                        obstruction = obstruct(equivatom, atom3, bset)
#                        offset = obstruction + nbas_vec[equivatom]
#                        #println("offset $offset")
#                        #println("$([basis:basis + 2*l]) -> $([basis + offset:basis + offset + 2*l])")
#                        println("$([basis + offset:basis + offset + 2*l]) -> $([basis + offset:basis + offset + 2*l])")
#                        #psalc[basis + offset:basis + offset + 2*l] = salc[basis:basis + 2*l]
#                        psalc[basis + offset:basis + offset + 2*l] = salc[basis + offset:basis + offset + 2*l]
#                        psalc
#                    else
#                        println("CONDITION 4")
#                        println("$y -> $atom3")
#                        obstruction1 = obstruct(equivatom, y, bset)
#                        offset1 = obstruction1 + nbas_vec[equivatom]
#                        obstruction2 = obstruct(y, atom3, bset)
#                        offset2 = obstruction2 + nbas_vec[equivatom]
#                        println("$([basis + offset1:basis + offset1 + 2*l]) -> $([basis + offset1 + offset2:basis + offset1 + offset2 + 2*l])")
#                        psalc[basis + offset1 + offset2:basis + offset1 + offset2 + 2*l] = salc[basis + offset1:basis + offset1 + 2*l] 
#                    end 
#                end
#            end
#            println("psalc $psalc")
#            if sum(broadcast(abs,(salc - psalc))) < 0.00001 || sum(broadcast(abs,(salc - -1*psalc))) < 0.00001
#                println("keep looking for pf")
#            else
#                break
#            end
#        end
#        #for atomb in SEAs[sea].set
#        #    println("atomb $atomb")
#        #end    
#        #for (pop, pclass) in enumerate(class_map)
#            #atom3 = amap[:op,1]
#    end
#    return salc
#end



#if sum(broadcast(abs, salc)) > 1e-14 && ct.characters[ir,:][1] > 1
#        println("This salc is degenerate and non-zero, pf needed $ir $salc")
#        println("SEAs[sea].set")
#        println(SEAs[sea].set)
#        println("this is where basis ended $basis")
#        for (pop, pclass) in enumerate(class_map)
#            for (x, y) in enumerate(SEAs[sea].set)
#                #println("$x $y")
#                atom3 = amap[:,pop][y]
#                if length(SEAs[sea].set) == 1
#                    println("It has to map onto itself, just check the a.m. rotation")
#                else
#                    println("$y mapped to $atom3")
#                    if atom3 == y
#                        println("CONDITION 1")
#                        if equivatom == atom3
#                            println("DA SAME $equivatom $y $atom3")
#                        else    
#                        #println("find the obstruction between this and the first atom of the list in [1, 3, 4] etc")
#                        obstruction = obstruct(equivatom,atom3, bset)
#                        offset = obstruction + nbas_vec[equivatom]
#                        #println("offset $offset")
#                        #println("$equivatom, but $y -> $atom3")
#                        #println("$([basis:basis + 2*l]) -> $([basis + offset:basis + offset + 2*l])")
#                        end
#                    elseif atom3 < y
#                        println("CONDITION 2")
#                        println("$y -> $atom3")
#
#                        obstruction = - obstruct(atom3, y, bset)
#                        println("obstruction $obstruction")
#                        offset = obstruction - nbas_vec[equivatom]
#                        println("offset $offset") 
#                        
#                        
#                        
#                        #println("$([basis:basis + 2*l]) -> $([basis + offset:basis + offset + 2*l])")
#                    else
#                        println("CONDITION 3")
#                        obstruction = obstruct(y, atom3, bset)
#                        offset = obstruction + nbas_vec[equivatom]
#                        println("offset $offset")
#                    end 
#                end
#            end
#            println("end of mapping")
#        end
#        #for atomb in SEAs[sea].set
#        #    println("atomb $atomb")
#        #end    
#        #for (pop, pclass) in enumerate(class_map)
#            #atom3 = amap[:op,1]
#    end