function path_angle(s::Real, path_fcn::Function)
    dp = ForwardDiff.derivative(path_fcn, s)
    return atan(dp[2], dp[1])
end

function path_curvature(s::Real, path_fcn::Function)
    dp = ForwardDiff.derivative(path_fcn, s)
    d2p = ForwardDiff.derivative(t -> ForwardDiff.derivative(path_fcn, t), s)
    dx, dy = dp
    d2x, d2y = d2p
    return (dx * d2y - dy * d2x) / (dx^2 + dy^2)^(3/2)    
end
