export regressor_matrix, parameter_vector

function regressor_matrix_full(x::Rn, model::ASVTowingModel)
    # Unpack state
    ψ = x[3]
    θ = x[4]
    vx = x[5]
    vy = x[6]
    ψ_dot = x[7]
    θ_dot = x[8]

    cψ = cos(ψ)
    sψ = sin(ψ)
    cθ = cos(θ)
    sθ = sin(θ)

    L = model.L

    # Regressor matrix
    # Damping coefficients
    Y_d = [
        vx*cψ^2 + vy*cψ*sψ   -vy*sψ^2 - vx*cψ*sψ              0  vy*cψ^2 - vx*cψ*sψ  vx*sψ^2 - vy*cψ*sψ              0  ψ_dot*cψ  -ψ_dot*sψ      0                vx - L*θ_dot*sθ;
        vy*sψ^2 + vx*cψ*sψ    vx*cψ^2 + vy*cψ*sψ              0  vy*cψ*sψ - vx*sψ^2  vy*cψ^2 - vx*cψ*sψ              0  ψ_dot*sψ   ψ_dot*cψ      0                vy + L*θ_dot*cθ;
                         0                     0  vx*cψ + vy*sψ                   0                   0  vy*cψ - vx*sψ         0          0  ψ_dot                              0;
                         0                     0              0                   0                   0              0         0          0      0  θ_dot*L^2 - L*vx*sθ + L*vy*cθ;
    ]
    # Damping times ocean current x (last column of the damping matrix has no effect)
    Y_x = [
         -cψ^2  cψ*sψ    0  cψ*sψ  -sψ^2   0    -1;
        -cψ*sψ  -cψ^2    0   sψ^2  cψ*sψ   0     0;
             0      0  -cψ      0      0  sψ     0;
             0      0    0      0      0   0  L*sθ;
    ]
    # Damping times ocean current y (last column of the damping matrix has no effect)
    Y_y = [
        -cψ*sψ    sψ^2    0   -cψ^2  cψ*sψ    0      0;
         -sψ^2  -cψ*sψ    0  -cψ*sψ  -cψ^2    0     -1;
             0       0  -sψ       0      0  -cψ      0;
             0       0    0       0      0    0  -L*cθ;
    ]

    return Y_d, Y_x, Y_y
end

function regressor_matrix(x::Rn, model::ASVTowingModel, damping::ASVDamping)
    Y_d, Y_x, Y_y = regressor_matrix_full(x, model)  
    
    # Columns of Y_d: D_11, D_21, D_31, D_12, D_22, D_32, D_13, D_23, D_33, b_T
    # Columns of Y_x: D_11, D_21, D_31, D_12, D_22, D_32, b_T
    if damping == FULL
        return hcat(Y_d, Y_x, Y_y)
    elseif damping == SWAY_YAW
        cols_Y_d = [1,5,6,8,9,10]  # D_11, D_22, D_32, D_23, D_33, b_T
        cols_Y_x = [1,5,6,7]       # D_11, D_22, D_32, b_T
        return hcat(Y_d[:, cols_Y_d], Y_x[:, cols_Y_x], Y_y[:, cols_Y_x])
    elseif damping == DIAGONAL
        cols_Y_d = [1,5,9,10]      # D_11, D_22, D_33, b_T
        cols_Y_x = [1,5,7]         # D_11, D_22, b_T
        return hcat(Y_d[:, cols_Y_d], Y_x[:, cols_Y_x], Y_y[:, cols_Y_x])
    else
        error("Unknown damping type")
    end
end

function parameter_vector(model::ASVTowingModel, damping::ASVDamping, V_c::Rn)
    Vx = V_c[1]
    Vy = V_c[2]

    ζ_d = [vec(model.D); model.b_T]  # Damping matrix + payload damping (10 elements)
    ζ_x = ζ_d[[1,2,3,4,5,6,10]] * Vx  # Damping matrix columns for ocean current x (7 elements)
    ζ_y = ζ_d[[1,2,3,4,5,6,10]] * Vy  # Damping matrix columns for ocean current y (7 elements)

    if damping == FULL
        return vcat(ζ_d, ζ_x, ζ_y)
    elseif damping == SWAY_YAW
        inds_d = [1,5,6,8,9,10]  # D_11, D_22, D_32, D_23, D_33, b_T
        inds_xy = [1,5,6,7]     # D_11, D_22, D_32, b_T
        return vcat(ζ_d[inds_d], ζ_x[inds_xy], ζ_y[inds_xy])
    elseif damping == DIAGONAL
        inds_d = [1,5,9,10]     # D_11, D_22, D_33, b_T
        inds_xy = [1,5,7]    # D_11, D_22, b_T
        return vcat(ζ_d[inds_d], ζ_x[inds_xy], ζ_y[inds_xy])
    else
        error("Unknown damping type")
    end
end
