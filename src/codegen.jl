
export cse, jlcode, genfunc, sympi

using SymPy
import Base: dot, sign, zero, one

sign(x::Sym) = sympy_meth(:sign, x)

const _symzero = oftype(Sym, 0)
zero(::Type{Sym}) = _symzero

const _symone = oftype(Sym, 1)
one(::Type{Sym}) = _symone

const sympi = Sym(sympy.pi)


import SymPy.cse

function cse(e::Sym; symbols=nothing, order="none")
    cseout = sympy_meth(:cse, e.x, symbols=symbols, order=order)
    (convert(Vector{(Sym,Sym)}, cseout[1]), cseout[2][1])
end

function cse(A::Array{Sym}; symbols=nothing, order="none")
    A = convert(SymMatrix, A)
    cseout = sympy_meth(:cse, A.x, symbols=symbols, order=order)
    (convert(Vector{(Sym,Sym)}, cseout[1]), convert(Array{Sym}, cseout[2][1]))
end



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
    [Expr(:(=), Expr(:ref, array, i, j), jlcode(exprmat[i, j])) for (i,j) in
        [(i,j) for i in 1:size(exprmat,1), j in 1:size(exprmat,2)]]


function genfunc(name, args::Vector, ret, ivs::Vector{(Sym,Sym)})
    body_ivs = [Expr(:(=), jlcode(iv), jlcode(e)) for (iv,e) in ivs]
    if string(name)[end] == '!'
        out = symbol(args[1])
        body_ret = assigntoarray(out, ret)
        push!(body_ret, Expr(:return, out))
    else
        body_ret = [Expr(:return, jlcode(ret))]
    end
    body = [body_ivs; body_ret]
    if name == :(->)
        Expr(:(->), Expr(:tuple, args...), Expr(:block, body...))
    else
        Expr(:function, Expr(:call, name, args...), Expr(:block, body...))
    end
end

genfunc{T}(name, args::Vector, ret::(Vector{(Sym,Sym)}, T)) = genfunc(name, args::Vector, ret[2], ret[1])


