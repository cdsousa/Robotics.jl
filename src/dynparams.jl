import Base: getindex, +, *
getindex(A::Symmetric, i::Integer) = A.S[i]
getindex(A::Symmetric, r::Range) = A.S[r]
+(A::Symmetric, B::Array) = A.S + B
+(A::Array, B::Symmetric) = B + A.S
*(A::Symmetric, B::Array) = A.S * B
*(A::Array, B::Symmetric) = B * A.S


import Base.convert
import DataStructures.OrderedDict
export OrderedDict

export symm, dynp_vect2odict, dynp_odict2vect, dynp_formatof, dynp_odict

function symm{T}(U::Vector{T})
    # construct Symmetric matrix from Vector of unique (row by row) values
    if length(U) == 6
        Symmetric([U[1] U[2] U[3]
                   U[2] U[4] U[5]
                   U[3] U[5] U[6]])
    else
        throw("Do not know how to make a symmetric matrix from a vector of size $(length(U))")
    end
end

convert{T}(::Type{Symmetric{T}}, V::AbstractVector{T}) = symm(V)


function dynp_vect2odict{T}(v::AbstractVector{T}, ndof::Integer, format::AbstractVector{(Symbol, DataType, Int)})
    # construct a OrderedDict of dynamic parameters from a vector and a given format
    ndynpdof = length(v)/ndof
    if !isinteger(ndynpdof)
        throw("Lenght of dynamic parameter vector is not multiple of DOF number")
    end

    d = OrderedDict{Symbol, Any}()
    for (s,typ) in format d[s] = typ[] end

    offset = 0
    for (s,typ,len) in format
        for i in 1:ndof
            idx = ndynpdof*(i-1)+offset+1
            if len>1 idx = idx:idx+len-1 end
            push!(d[s], convert(typ, v[idx]))
        end
        offset += len
    end

    d
end



flat(S::Symmetric) = vcat([S.S[i:end, i] for i in 1:size(S,1)]...)
flat(M::Matrix) = vec(M)
flat(V::Vector) = V
flat{T}(x::T) = T[x]

function dynp_odict2vect(d::OrderedDict{Symbol})
    # construct a flat vector from an OrderedDict of dynamic parameters
    ndof = length(d[:m])
    typ = eltype(d[:m])
    v = Array{typ, 1}[[] for _ in 1:ndof]
    for (s,x) in d
        for i in 1:ndof append!(v[i], flat(x[i])) end
    end
    vcat(v...)
end

# get the format of an OrderedDict of dynamic parameters
dynp_formatof(d::OrderedDict{Symbol}) =
    (Symbol, DataType, Int)[(s, eltype(x), length(flat(x[1]))) for (s,x) in d]


dynp_odict(;args...) = OrderedDict([args...])

