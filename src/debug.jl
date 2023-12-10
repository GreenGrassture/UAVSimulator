using ControlSystems
using DifferentialEquations
using LinearAlgebra
using PlotlyJS
using Zygote
using Statistics
using Memoize
using LRUCache
using Distributions
using LinearAlgebra

# My modules
include("Control.jl")
include("PlottingFuns.jl")
include("InsertDeleteRowsCols.jl")
include("SystemDynamics.jl")
include("RotationFuns.jl")
include("NoiseFunctions.jl")
include("Forces.jl")
include("QuadConstants.jl")

sInit = vec([0. 0. 0.])
vInit = vec([0. 0. 0.])
qInit = euler2quat(0, 0, 0)
ωInit = vec([0. 0. 0.])
xInit = vec([sInit; vInit; qInit; ωInit])
tInit0 = 0.
tInit1 = 5.
tSample = 0.01
initParams = Dict("quadConsts"=>CFConsts,
                    "controller"=>controllerNone(),
                    "forces"=>externalForcesOnlyG,
                    "reference"=>referenceSignalConstant) 
uInit = initParams["controller"](xInit, initParams, tInit0)
dx0 = innerDFun(xInit, uInit, initParams, tInit0)

A = [0.0  0.0  0.0  1.0  0.0  0.0   0.0    0.0     0.0   0.0  0.0  0.0  0.0;
0.0  0.0  0.0  0.0  1.0  0.0   0.0    0.0     0.0   0.0  0.0  0.0  0.0;
0.0  0.0  0.0  0.0  0.0  1.0   0.0    0.0     0.0   0.0  0.0  0.0  0.0;
0.0  0.0  0.0  0.0  0.0  0.0   0.0    0.0   -19.62  0.0  0.0  0.0  0.0;
0.0  0.0  0.0  0.0  0.0  0.0   0.0   19.62    0.0   0.0  0.0  0.0  0.0;
0.0  0.0  0.0  0.0  0.0  0.0  39.24   0.0     0.0   0.0  0.0  0.0  0.0;
0.0  0.0  0.0  0.0  0.0  0.0   0.0    0.0     0.0   0.0  0.0  0.0  0.0;
0.0  0.0  0.0  0.0  0.0  0.0   0.0    0.0     0.0   0.0  0.5  0.0  0.0;
0.0  0.0  0.0  0.0  0.0  0.0   0.0    0.0     0.0   0.0  0.0  0.5  0.0;
0.0  0.0  0.0  0.0  0.0  0.0   0.0    0.0     0.0   0.0  0.0  0.0  0.5;
0.0  0.0  0.0  0.0  0.0  0.0   0.0    0.0     0.0   0.0  0.0  0.0  0.0;
0.0  0.0  0.0  0.0  0.0  0.0   0.0    0.0     0.0   0.0  0.0  0.0  0.0;
0.0  0.0  0.0  0.0  0.0  0.0   0.0    0.0     0.0   0.0  0.0  0.0  0.0]

B =[0.0       0.0              0.0      0.0;
    0.0       0.0              0.0      0.0;
    0.0       0.0              0.0      0.0;
    0.0       0.0              0.0      0.0;
    0.0       0.0              0.0      0.0;
    31.25     0.0              0.0      0.0;
    0.0       0.0              0.0      0.0;
    0.0       0.0              0.0      0.0;
    0.0       0.0              0.0      0.0;
    0.0       0.0              0.0      0.0;
    0.0   62500.0              0.0      0.0;
    0.0       0.0          62500.0      0.0;
    0.0       0.0              0.0  34482.8]

microCallback = PeriodicCallback(microcontroller!, tSample, initial_affect=true)
sol = simSys(xInit, initParams, (tInit0, tInit1), microCallback, false, false);

sQ = 100*[1., 1., 1.]
vQ = 1.0*[1., 1., 1.]
qQ = 1.0*[1., 1., 1., 1.]
ωQ = 1.0*[1., 1., 1.]
Q = Matrix(Diagonal([sQ; vQ; qQ; ωQ]))
R = 1.0*Matrix(I(4))

controllerLQR = tryLQR(A, B, Q, R, lqr)