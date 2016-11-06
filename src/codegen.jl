
export cse, jlcode, genfunc, sympi

using SymPy
import Base: dot, sign #, zero, one

sign(x::Sym) = sympy_meth(:sign, x)

# const _symzero = oftype(Sym, 0)
# zero(::Type{Sym}) = _symzero
#
# const _symone = oftype(Sym, 1)
# one(::Type{Sym}) = _symone

const sympi = Sym(sympy[:pi])


# import SymPy.cse
#
# function cse(e::Sym; symbols=nothing, order="none")
#     cseout = sympy_meth(:cse, e.x, symbols=symbols, order=order)
#     (convert(Vector{Tuple{Sym,Sym}}, cseout[1]), cseout[2][1])
# end
#
# function cse(A::Array{Sym}; symbols=nothing, order="none")
#     A = convert(SymMatrix, A)
#     cseout = sympy_meth(:cse, A.x, symbols=symbols, order=order)
#     (convert(Vector{Tuple{Sym,Sym}}, cseout[1]), convert(Array{Sym}, cseout[2][1]))
# end
retype_cse(symcode) = ([symcode[1]...], [x for x in symcode[2]])
export retype_cse


function jlcode(expr::Sym)
    codestr = SymPy.sympy_meth(:ccode, expr)
    codestr = replace(codestr, "pow", "^")
    parse(codestr)
end

jlcode(exprvec::Vector{Sym}) = Expr(:vcat, [jlcode(e) for e in exprvec]...)

jlcode(exprmat::Matrix{Sym}) =
    Expr(:vcat, [Expr(:row, [jlcode(exprmat[i,j]) for j in 1:size(exprmat,2)]...) for i in 1:size(exprmat,1)]...)


assigntoarray(array::Symbol, exprvec::Vector{Sym}) =
    [Expr(:(=), Expr(:ref, array, i), jlcode(exprvec[i])) for i in 1:length(exprvec)]

assigntoarray(array::Symbol, exprmat::Matrix{Sym}) =
    [Expr(:(=), Expr(:ref, array, i, j), jlcode(exprmat[i, j])) for j in 1:size(exprmat,2) for i in 1:size(exprmat,1)]


function genfunc(name, args::Vector, ret, ivs::Vector{Tuple{Sym,Sym}}; inbounds=false, overret=nothing)
    body_ivs = [Expr(:(=), jlcode(iv), jlcode(e)) for (iv,e) in ivs]
    if string(name)[end] == '!'
        out = Symbol(args[1])
        body_ret = assigntoarray(out, ret)
        push!(body_ret, Expr(:return, out))
    else
        if overret === nothing
            body_ret = [Expr(:return, jlcode(ret))]
        else
            body_ret = [Expr(:return, Expr(overret.head, overret.args..., jlcode(ret)))]
        end
    end
    body = Expr(:block, [body_ivs; body_ret]...)
    if inbounds
        body = Expr(:macrocall, Symbol("@inbounds"), body)
    end
    if name == :(->)
        Expr(:(->), Expr(:tuple, args...), body)
    else
        Expr(:function, Expr(:call, name, args...), body)
    end
end

genfunc{T}(name, args::Vector, ret::Tuple{Vector{Tuple{Sym,Sym}}, T}; inbounds=false, overret=nothing) = genfunc(name, args::Vector, ret[2], ret[1]; inbounds=inbounds, overret=overret)
