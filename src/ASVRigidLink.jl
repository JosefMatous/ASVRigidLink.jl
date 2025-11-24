module ASVRigidLink

using LinearAlgebra
import ForwardDiff, DiffResults
import DocStringExtensions: FIELDS

const Rn = AbstractVector{<:Real}
const Rmxn = AbstractMatrix{<:Real}
export Rn, Rmxn

include("saturation.jl")
include("sign.jl")
include("model.jl")
include("asv_control.jl")
include("straightline/guidance.jl")
include("straightline/simulation.jl")
include("curved/path.jl")
include("curved/control_lin.jl")
include("curved/guidance.jl")
include("curved/simulation.jl")
include("visualization.jl")

end # module ASVRigidLink
