# This script includes a bunch of the other files and libraries I need.  I split it off
# from the main script so the imports don't run over and over
print("Loading Julia modules\n")
using BenchmarkTools
using Profile
using PProf

using ControlSystems
using DifferentialEquations
using DiffEqCallbacks
using LinearAlgebra
using PlotlyJS
using Zygote
using Statistics
using Memoization
using LRUCache
using Distributions
using MatrixEquations
using StaticArrays
using ProgressMeter 

print("Loading my modules\n")
# My modules
include("Control.jl")
include("PlottingFuns.jl")
include("InsertDeleteRowsCols.jl")
include("SystemDynamics.jl")
include("RotationFuns.jl")
include("NoiseFunctions.jl")
include("Forces.jl")
include("QuadConstants.jl")
include("ReferenceSignals.jl")


##### Convenience Functions ################
function unit(d, idx)
    # Creates a unit vector in R^d with entry in position idx
    unitVec = zeros(d)
    unitVec[idx] = 1
    return unitVec
end
#############################################

# Run this to cause julia to compile/run all the important functions with the
# proper arguments
sInit = [0.0, 0.0, 0.0]
vInit = [0.0, 0.0, 0.0]
qInit = euler2quat(0, 0, 0)
ωInit = [0.0,  0.0,  0.0]
uInit = [0.0, 0.0, 0.0, 0.0]
xInit = vec([sInit; vInit; qInit; ωInit; uInit])
tInit0 = 0.0
tInit1 = 0.05
tSample = 0.01
initParams = (quadConsts=CFConsts,
                controller=controllerNone(),
                u=[0.0, 0.0, 0.0, 0.0],
                forces=externalForcesOnlyG,
                reference=referenceSignalConstant) 
initParams[:u] .= initParams[:controller](xInit, initParams, tInit0)
dx0 = autodiffDFun(xInit, initParams, tInit0)

AInit, BInit = linearize(autodiffDFun, xInit, initParams, tInit0)
microCallback = PeriodicCallback(microcontroller!, tSample, initial_affect=true)
prob = ODEProblem(innerDFun!, xInit, (tInit0, tInit1), initParams, callback=microCallback, save_everystep=false)
sol = solve(prob, AutoTsit5(Rosenbrock23()), progress=true)

initializedQuad = true;