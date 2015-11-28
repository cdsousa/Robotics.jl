export vectosymm, symmtovec, dynp_init, dynp_vect2dict, dynp_dict2vect


vectosymm(U::AbstractVector) =
    [U[1] U[2] U[3]
     U[2] U[4] U[5]
     U[3] U[5] U[6]]

symmtovec(S::AbstractMatrix) =
    [S[1,1], S[1,2], S[1,3], S[2,2], S[2,3], S[3,3]]


immutable DynParmType{T} end

flatter(::Type{DynParmType{:scalar}}) = x->[x]
flatter(::Type{DynParmType{:vec3}}) = identity
flatter(::Type{DynParmType{:symm3}}) = symmtovec

expander(::Type{DynParmType{:scalar}}) = v->v[1]
expander(::Type{DynParmType{:vec3}}) = identity
expander(::Type{DynParmType{:symm3}}) = vectosymm

dynptype(::Type{DynParmType{:scalar}}, T) = T
dynptype(::Type{DynParmType{:vec3}}, T) = Vector{T}
dynptype(::Type{DynParmType{:symm3}}, T) = Matrix{T}

dynplen(::Type{DynParmType{:scalar}}) = 1
dynplen(::Type{DynParmType{:vec3}}) = 3
dynplen(::Type{DynParmType{:symm3}}) = 6



function dynp_init(;args...)
    ndof = length(args[1][2][1])
    vtype = eltype(eltype(args[1][2][1]))

    format = Tuple{Symbol, Symbol}[]
    v = Array{vtype, 1}[[] for _ in 1:ndof]

    for (s, (x, typ)) in args
        flat = flatter(DynParmType{typ})
        for i in 1:ndof append!(v[i], flat(x[i])) end
        push!(format, (s, typ))
    end

    vcat(v...), format
end



function dynp_vect2dict{T}(v::AbstractVector{T}, ndof::Int, format::AbstractVector{Tuple{Symbol, Symbol}})
    ndynpdof = div(length(v), ndof)
    if !isinteger(ndynpdof)
        throw("Lenght of dynamic parameter vector is not multiple of DOF number")
    end

    d = Dict{Symbol, Any}()

    for (s, typ) in format
        d[s] = dynptype(DynParmType{typ}, T)[]
    end

    offset = 0
    for (s, typ) in format
        len = dynplen(DynParmType{typ})
        expand = expander(DynParmType{typ})
        for i in 1:ndof
            idx = ndynpdof*(i-1)+offset
            range = idx+1:idx+len
            push!(d[s], expand(v[range]))
        end
        offset += len
    end

    d
end


function dynp_dict2vect(d, format::AbstractVector{Tuple{Symbol, Symbol}})
    ndof = length(d[format[1][1]])
    vtype = eltype(eltype(d[format[1][1]]))

    v = Array{vtype, 1}[[] for _ in 1:ndof]

    for (s, typ) in format
        x = d[s]
        flat = flatter(DynParmType{typ})
        for i in 1:ndof append!(v[i], flat(x[i])) end
    end

    vcat(v...)
end
