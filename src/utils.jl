
array{T}(x::T...) = T[x...]

skew(v) =
    [    0 -v[3]  v[2]
      v[3]     0 -v[1]
     -v[2]  v[1]     0]

