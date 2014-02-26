println("--------- statring Robotics.jl tests ----------")

using Base.Test

using Robotics
@test isdefined(:Robotics)
@test typeof(Robotics) === Module

include("checkpuma.jl")

println("✓ all Robotics.jl tests passed")

