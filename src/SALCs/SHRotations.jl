#This code generates the rotation matrices for spherical harmonics by a recursion relation
#found in "Rotation Matrices for Real Spherical Harmonics. Direct Determination by Recursion"
#J. Ivanic and K. Rudenberg: doi/10.1021/jp953350u


# return identity matrix
function identity_matrix()
    a = zeros(Int8, (3,3))
    for i = 1:3
        a[i,i] = 1
    end    
    return a
end

#generates u, w, and v coefficients found in Table 1 of reference
function UWVCoefficient(l, m1, m2)
    δ=  ==(0,m1)
    if abs(m2) < l
        denom = (l + m2)*(l - m2)
    else
        denom = (2.0*l)*(2.0*l-1.0)
    end
    u = (((l + m1)*(l - m1)) /denom)^(1/2)
    v = 0.5*(((1 + δ)*(l + abs(m1) - 1)*(l + abs(m1))/ denom)^(1/2))*(1 - 2*δ)
    w = - 0.5*(((l - abs(m1) - 1)*(l - abs(m1))/ denom)^(1/2))*(1 - δ)
    return u, v, w
end

#generates function U found in Table 2 of reference
function Ufun(l, m1, m2, rot, Rsh)
    #println("reference 1")
    U = Pfun(l, 0, m1, m2, rot, Rsh)
    return U
end

#generates function V found in Table 2 of reference, with a sign correction (***DERIVE**)
function Vfun(l, m1, m2, rot, Rsh)
    
    if m1 == 0
        #println("reference 2")
        V = Pfun(l, 1, 1, m2, rot, Rsh) + Pfun(l, - 1, - 1, m2, rot, Rsh)
    elseif m1 == 1
        #println("reference 3")
        V = sqrt(2)*Pfun(l,  1, 0, m2, rot, Rsh)
    elseif m1 == -1
        #println("reference 4")
        V = sqrt(2)*Pfun(l, -1, 0, m2, rot, Rsh) 
    elseif m1 > 0
        #println("reference 5")
        V = Pfun(l, 1, m1 - 1, m2, rot, Rsh) - Pfun(l, -1, -m1 + 1, m2, rot, Rsh)
    else
        #println("reference 6")
        V = Pfun(l, 1, m1 + 1, m2, rot, Rsh) + Pfun(l, -1, -m1 - 1, m2, rot, Rsh)
    end
    return V
end

#generates function W found in Table 2 of reference
function Wfun(l, m1, m2, rot, Rsh)
    if m1 > 0

        #println("reference 7")
        W = Pfun(l, 1, m1 + 1, m2, rot, Rsh) + Pfun(l, -1, -m1 -1, m2, rot, Rsh)
    else m1 < 0
        #println("reference 8")
        W = Pfun(l, 1, m1 - 1, m2, rot, Rsh) - Pfun(l, -1, -m1 +1, m2, rot, Rsh) 
    end
    return W
end

#function P in terms of matrix R_ij and R^(l-1)
function Pfun(l, i, m1, m2, rot, Rsh)
    
    rsh = Rsh[l]
    dl = size(Rsh[l])[1]
    ol = Int((dl -1)/2)
    if m2 == l
        
        P1 = rot[i + 2, 3] * rsh[m1 + ol + 1, l  + ol]
        P2 = rot[i + 2, 1] * rsh[m1 + ol + 1, 2 - l + ol]
        P = P1 - P2
    elseif m2 == -l
        P1 = rot[i + 2, 3] * rsh[m1 + ol + 1, 2 - l + ol]
        P2 = rot[i + 2, 1] * rsh[m1 + ol + 1, l + ol]
        P = P1 + P2
    else
        P = rot[i + 2, 2] * rsh[m1 + ol + 1, m2 + ol + 1]
    end
    return P
end

function adapt(rot)
    rrot = zeros(Float64,3,3)
    #println(rot)
    rrot[1,1] = rot[2,2]
    rrot[1,2] = rot[2,3]
    rrot[1,3] = rot[2,1]
    rrot[2,1] = rot[3,2]
    rrot[2,2] = rot[3,3]
    rrot[2,3] = rot[3,1]
    rrot[3,1] = rot[1,2]
    rrot[3,2] = rot[1,3]
    rrot[3,3] = rot[1,1]
    #println(rot)
    return rrot
end
function generateRotations(Lmax, rot)
    Rsh = Matrix{Float64}[]
    rrot = adapt(rot)
    #rrot = rot
    #println("Rotates rot ", rrot)
    push!(Rsh, identity_matrix())
    l = 1
    while l < Lmax + 1
        if l == 1
            push!(Rsh,rrot) # rrot
        end
        if l > 1
        R = zeros(Float64, (2*l + 1,2 * l + 1))
            for m1 in -l:l
                for m2 in -l:l
                    u, v, w = UWVCoefficient(l, m1, m2)
                    if u != 0
                        u *= Ufun( l, m1, m2, rrot, Rsh)
                    end 
                    if v != 0
                        v *= Vfun( l, m1, m2, rrot, Rsh) 
                    end
                    if w != 0
                        w *= Wfun( l, m1, m2, rrot, Rsh)
                    end
                    R[m1 + l + 1, m2 + l + 1] = u + v + w
                end
            end
            push!(Rsh, R)
        end 
        l += 1
    end
    if Lmax >= 1
        Rsh[2] = rot
    end
    return Rsh
end
