
export standard_dh_transfmat, modified_dh_transfmat, inv_transfmat,
    interlink_transfs_dhstd, interlink_transfs_dhmod, transfs_tobase,
    interlink_lie_S_rotjoint_dhstd

standard_dh_transfmat(α,a,d,θ) =
    [cos(θ) -sin(θ)*cos(α)  sin(θ)*sin(α) a*cos(θ)
     sin(θ)  cos(θ)*cos(α) -cos(θ)*sin(α) a*sin(θ)
     0              sin(α)         cos(α)        d
     0                   0              0        1]


modified_dh_transfmat(α,a,d,θ) =
    [       cos(θ)       -sin(θ)       0         a
     sin(θ)*cos(α) cos(θ)*cos(α) -sin(α) -sin(α)*d
     sin(θ)*sin(α) cos(θ)*sin(α)  cos(α)  cos(α)*d
                 0             0       0         1]

inv_transfmat(A) =
    [A[1:3,1:3]'    -A[1:3,1:3]'*A[1:3,4]
                            A[end:end, :]]

interlink_transfs_dhstd(params) = [standard_dh_transfmat(dhp...) for dhp in params]
interlink_transfs_dhmod(params) = [modified_dh_transfmat(dhp...) for dhp in params]

function transfs_tobase(Ti)
    T = [Ti[1]]
    for i in 2:length(Ti)
        push!(T, T[end] * Ti[i])
    end
    T
end


### Lie algebra robotics:

function unskew_SE3(g)
    w = [g[3, 2], g[1, 3], g[2, 1]]
    v = g[1:3, 4]
    [w; v]
end

_Lie_S_rot(α,a,d,θ) = [      0 -cos(α) sin(α)         0
                        cos(α)       0      0  a*cos(α)
                       -sin(α)       0      0 -a*sin(α)
                             0       0      0         0]


lie_S_rotjoint(α,a,d,θ) = unskew_SE3(_Lie_S_rot(α,a,d,θ))

interlink_lie_S_rotjoint_dhstd(params) = [lie_S_rotjoint(dhp...) for dhp in params]
