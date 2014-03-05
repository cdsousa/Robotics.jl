
export observmat, dependcol


function observmat(T, m, n, func, otherargs=())
    H = zeros(T, m, n)
    params = zeros(T, n)
    for i in 1:n
        params[i] = one(T)
        H[:, i] = func(params, otherargs...)
        params[i] = zero(T)
    end
    H
end


function dependcol(Hs, prec=10)
    np = size(Hs, 2)

    Q1, R1 = qr(Hs)
    R1_diag = round(diag(R1), prec)

    dbi = Int[]
    ddi = Int[]
    for (i, e) in enumerate(R1_diag)
        if e != 0
            push!(dbi, i)
        else
            push!(ddi, i)
        end
    end
    dbn = length(dbi)

    P = eye(Int, np)[:, [dbi; ddi]]
    Pb = P[:, 1:dbn]
    Pd = P[:, dbn+1:end]

    _ , Rbd1 = qr(Hs*P)
    Rb1 = Rbd1[1:dbn, 1:dbn]
    Rd1 = Rbd1[1:dbn, dbn+1:end]

    Kd = round(inv(Rb1) * Rd1, prec)

    Pb, Pd, Kd
end

