
import DataStructures.OrderedDict

export dynparams2vect, dynparams_from_vect, dynparams_canonical


symm2uniq(A) = vcat([A[i:end, i] for i in 1:size(A,1)]...)

uniq2symm_3x3(x) =
    [x[1] x[2] x[3]
     x[2] x[4] x[5]
     x[3] x[5] x[6]]


function dynparams2vect(;args...)
    p = OrderedDict{Symbol, Any}(args)
    m = pop!(p, :m)
    dof = length(m)

    if haskey(p, :l)
        l = pop!(p, :l)
    else
        r = pop!(p, :r)
        l = [m[i]*r[i] for i in 1:dof]
    end

    if haskey(p, :sL)
        sL = pop!(p, :sL)
    else
        if haskey(p, :L)
            L = pop!(p, :L)
        else
            if haskey(p, :sI)
                sI = pop!(p, :sI)
                I = [uniq2symm_3x3(sI[i]) for i in 1:dof]
            else
                I = pop!(p, :I)
            end
            L = [I[i] + m[i]*skew(r[i])'*skew(r[i]) for i in 1:dof]
        end
        sL = [symm2uniq(L[i]) for i in 1:dof]
    end

    order = pop!(p, :order, :Khalil)
    if order == :Khalil
        vcat([vcat(sL[i]..., l[i]..., m[i]) for i in 1:dof]...)
    else
        throw("dynparams: unknown order '$order'")
    end
end


function dynparams_from_vect(paramsvect, dof, format=[:L, :l, :m])
    params = Dict()
    np = length(paramsvect)/dof
    if isinteger(np) np=int(np) else throw("number of parameters not multiple of dof number") end

    o = 0
    po = OrderedDict()
    for p in format
        if p==:L
            params[p] = [uniq2symm_3x3(paramsvect[ i*np+o+1 : i*np+o+6 ]) for i in 0:dof-1]
            o += 6
        elseif p==:l
            params[p] = [paramsvect[ i*np+o+1 : i*np+o+3 ] for i in 0:dof-1]
            o += 3
        else
            params[p] = [paramsvect[ i*np+o+1 ] for i in 0:dof-1]
            o += 1
        end
    end

    params
end

function dynparams_canonical(;args...)
    p = OrderedDict{Symbol, Any}(args)
    dof = length(p[:m])
    v = dynparams2vect(;args...)
    dynparams_from_vect(v, dof)
end

