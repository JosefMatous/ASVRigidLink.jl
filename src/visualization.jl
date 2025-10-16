using Plots

export SystemDimensions, plot_system!

"""

Dimensions of the plotted system

# Fields
$(FIELDS)
"""
struct SystemDimensions
    "ASV main body length"
    asv_body_length::Real
    "ASV bow length"
    asv_bow_length::Real
    "ASV width"
    asv_width::Real
    "ASV center of gravity offset"
    asv_xcg::Real
    "Cable length"
    cable_length::Real
end

function SystemDimensions(asv_body_length::Real, asv_bow_length::Real, asv_width::Real, model::ASVTowingModel)
    return SystemDimensions(asv_body_length, asv_bow_length, asv_width, model.x_G, model.L)
end

function SystemDimensions(asv_body_length::Real, asv_bow_length::Real, asv_width::Real, sim::SimulationParameters)
    return SystemDimensions(asv_body_length, asv_bow_length, asv_width, sim.model)
end

function plot_system!(plt, x::Rn, dimensions::SystemDimensions;
    asv_params=(), cable_params=(), payload_params=())

    # Unpack state
    x_asv = x[1]
    y_asv = x[2]
    ψ = x[3]
    θ = x[4]

    # ASV vertices in body frame
    L_b = dimensions.asv_body_length
    L_f = dimensions.asv_bow_length
    W = dimensions.asv_width
    x_cg = dimensions.asv_xcg
    x_center = (3*L_b + L_f) / 5 + x_cg
    # Vertices of the polygon are in this order:
    #  [bow, front right, aft right, aft left, front left]
    X_body = [L_b + L_f - x_center, L_b - x_center, -x_center, -x_center, L_b - x_center]
    Y_body = [0, W/2, W/2, -W/2, -W/2]

    # Transform to inertial frame
    X_asv = x_asv .+ cos(ψ) .* X_body .- sin(ψ) .* Y_body
    Y_asv = y_asv .+ sin(ψ) .* X_body .+ cos(ψ) .* Y_body

    # Payload position
    L = dimensions.cable_length
    x_payload = x_asv + L * cos(θ)
    y_payload = y_asv + L * sin(θ)

    # Plot
    plot!(plt, Y_asv, X_asv; seriestype=:shape, asv_params...)
    plot!(plt, [y_asv, y_payload], [x_asv, x_payload]; cable_params...)
    scatter!(plt, [y_payload], [x_payload]; payload_params...)

    return plt
end
