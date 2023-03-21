struct SALC{F<:Real,I<:Integer}
    coeffs::Vector{F}
    irrep::String
    atom::I
    bfxn::I
    sh::I # This is actually l+1, kinda weird but it's so we can index into things, I'll change this later, but I don't want to reset julia rn
    ml::I # 1 <= ml <= 2l+1
    i::I
    r::I
    γ::F
end

function string_repr(salcs::Vector{SALC{F,I}}) where {F,I}
    strang = ""
    for (sidx,salc) in enumerate(salcs)
        strang *= "SALC $sidx: irrep = $(salc.irrep), subspecies = $(salc.i), r = $(salc.r), γ = $(salc.γ), "
        strang *= "atom = $(salc.atom), basis function = $(salc.bfxn), shell = $(salc.sh), m_l = $(salc.ml)\n"
        #strang *= "\n$(salc.coeffs)\n"
    end
    return strang
end

function show(io::IO, ::MIME"text/plain", salcs::Vector{SALC{F,I}}) where {F,I}
    print(io, string_repr(salcs))
end
