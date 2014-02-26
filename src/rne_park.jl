
export rne_park_forward, rne_park_backward, rne_park


function Adj(G, g)
    R = G[1:3, 1:3]
    p = G[1:3, 4]
    [        R  zeros(eltype(G),3,3)
     skew(p)*R                     R] * g
end

function Adjdual(G, g)
    R = G[1:3, 1:3]
    p = G[1:3, 4]
    [        R zeros(eltype(G),3,3)
     skew(p)*R                   R]' * g
end

function adj(g, h)
    wg = g[1:3, 1]
    vg = g[4:6, 1]
    [skew(wg) zeros(eltype(g),3,3)
     skew(vg)             skew(wg)] * h
end


function adjdual(g, h)
    wg = g[1:3, 1]
    vg = g[4:6, 1]
    [skew(wg) zeros(eltype(g),3,3)
     skew(vg)             skew(wg)]' * h
end


function rne_park_forward(q, dq, ddq, g, Ti_inv, S, ifunc=identity)

    typ = eltype(q)
    dof = length(q)

    V = [Array(typ, 6) for _ in 1:dof+1]
    dV = [Array(typ, 6) for _ in 1:dof+1]

    V[1] = zeros(typ, 6)
    dV[1] = [zeros(typ, 3); -g]

    # Forward
    for i in 1:dof
        V[i+1] = ifunc(Adj(Ti_inv[i], V[i])) + ifunc(S[i]*dq[i])
        V[i+1] = ifunc(V[i+1])

        dV[i+1] = ifunc(S[i]*ddq[i]) + ifunc(Adj(Ti_inv[i], dV[i])) +
                ifunc(adj(ifunc(Adj(Ti_inv[i], V[i])), ifunc(S[i]*dq[i])))
        dV[i+1] = ifunc(dV[i+1])
    end

    V, dV
end

function rne_park_backward(Ti_inv, S, L, l, m, V, dV, ifunc=identity)

    typ = eltype(m)
    dof = length(m)

    # extend Ti_inv so that Ti_inv[dof+1] return identity
    Ti_inv = array(Ti_inv..., eye(typ, 4))

    F = [Array(typ, 6) for _ in 1:dof+1]

    F[dof+1] = zeros(typ, 6)

    tau = zeros(typ, dof)

    # Backward
    for i in reverse(1:dof)
        Llm = [       L[i]        skew(l[i])
               -skew(l[i]) m[i].*eye(typ, 3)]

        F[i] = Adjdual(Ti_inv[i+1], F[i+1]) + Llm*dV[i+1] - adjdual(V[i+1],  Llm*V[i+1])
        F[i] = ifunc(F[i])

        tau[i] = ifunc( dot(S[i], F[i]) )
    end

    tau
end


function rne_park(q, dq, ddq, Ti_inv, S, g, L, l, m, ifunc=identity)
    V, dV = rne_park_forward(q, dq, ddq, g, Ti_inv, S, ifunc)
    tau = rne_park_backward(Ti_inv, S, L, l, m, V, dV, ifunc)
    tau
end
