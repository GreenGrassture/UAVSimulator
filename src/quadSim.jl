# Check if the initialization script has run.  If not, run it
print("Starting\n")
if (! @isdefined initializedQuad)||(initializedQuad==false)
    print("Initializing workspace\n")
    include("initScript.jl")
    print("Initialization complete\n")
end
# Notes:
# I'm running into a problem with the quaternion representation of the orientation - at hover conditions, the row of the A matrix corresponding to q0 is all zeros.
# This means that the state is not controllable, and it leads to a singularity exception when trying to solve the CARE for the purpose of LQR control.  Somehow I need to
# account for this issue.  My first thought is to delete q0 in the linear model when I'm using it for the CARE.  Then I can take the resulting K matrix and patch in an 
# additional row and column of zeros.  My biggest concern is that it might not be unique to q0.  Maybe the other elements of q become uncontrollable at other orientations.
# It seems like the underlying issue is that the quaternion representation is over-determined.  q0^2 = 1 - q1^2 + q2^2 + q3^2

# The nonlinearities are only in the orientation.  The position and velocity control could really be handled completely by a linear controller, with a nonlinear controller
# handling the leftover states.  Though velocity is still coupled to angle via mgsin(θ) etc.

# Steps:
# 1) Define nonlinear system along with a dummy controller to get open loop system
# 2) Use automatic differentiation to get linearized system (without control inputs)
# 3) Write LQR controller for linearized system
# 4) Plug LQR controller into nonlinear system

# Also, I need to implement a kalman filter and and a random noise generator that plays nicely with the ODE solver

# First define an array to hold all of the callbacks that might be created while initializing
# other things

# Define the nominal operating point and get the linearized model at that point
print("Deriving linear model\n")
sLin = [0.0, 0.0, 0.0,] # Displacement
vLin = [0.0, 0.0, 0.0] # Velocity
qLin = euler2quat(0,0,0) # Orientation
ωLin = [0.0, 0.0, 0.0] # Angular velocity
uLin = [0.0, 0.0, 0.0, 0.0] # Controller output
xLin = [sLin; vLin; qLin; ωLin; uLin] # Complete state vector

linParams = (quadConsts=CFConsts,
             controller=controllerNone(),
             forces=externalForcesOnlyG,
             reference=referenceSignalConstant)
tLin = 0.0

# Compute the initial control signal.  We need to have the full state to do this, which is why
# we initialize it with zeros above, then replace it with the actual value here
uLin = linParams[:controller](xLin, linParams, tLin)
xLin[14:17] = uLin

# Define the cost function J = Σ x'Qx + u'Ru
sQ = 10*[1., 1., 1.]
vQ = 1.0*[1., 1., 1.]
qQ = 1.0*[1., 1., 1., 1.]
ωQ = 1.0*[1., 1., 1.]
Q = Matrix(Diagonal([sQ; vQ; qQ; ωQ]))
R = 1.0*Matrix(I(4))
# The frequency of the microcontroller
tSample = 0.05
# Create a controller for that operating point
controllerLQR = makeLQR(autodiffDFun, xLin, linParams, tLin, Q, R, tSample)
print("Created controller\n")

# Define initial conditions for the ODE
print("Specifying initial conditions\n")
s0 =  [0.0, 0.0, 0.0] # position
v0 =  [0.0, 0.0, 0.0] # velocity
q0 = euler2quat(0.0, 0.0, 0.0) # orientation - defined with euler angles for convenience, converted to quaternions for internal use
ω0 =  [0.0, 0.0, 0.0] # angular velocity
u0 = [0.0, 0.0, 0.0, 0.0]
x0 = [s0; v0; q0; ω0; u0] # combine into initial state vector
# define simulation time
t0 = 0.0
tf = 5.0
tSpan = (t0, tf)
# Set the problem parameters
noiseMeans = 0.05*ones(6) #zeros(6)
noiseCovs = Diagonal(0.01*ones(Float64, 6))
noiseBandwidth = 20
numNoiseInterp = 11
externalForcesGravAndNoise = createForceFunctionGravAndNoise(noiseMeans, noiseCovs, noiseBandwidth, numNoiseInterp)
#
ref = makeReference_Constant([3, 0, 0], 
                             [0, 0, 0 ], 
                             euler2quat([0, 0, 0]), 
                             [0, 0, 0])
#
#=
ref = makeReference_Sine([0.0 0.0 0.0],
                         [1.0 0.0 1.0],
                         0.2)
=#  

ODEParams = (quadConsts=CFConsts,
             controller=controllerLQR,
             forces=externalForcesGravAndNoise,
             reference=ref)                                           
# Constructing callbacks
# Ensures that quaternions maintain unit length
quatCallback = ManifoldProjection(quaternionError)
# Periodically calls whatever code the microcontroller would run
microCallback = PeriodicCallback(microcontroller!, tSample, initial_affect=true)
# Early stopping condition, in case system becomes unstable and causes integrator to chug
function terminateCondition(u, t, integrator)
    # Check if any of the angular velocities have become too large
    @views bigVel = maximum(abs.(u[11:13])) > 100
    return bigVel
end
terminateCallback = DiscreteCallback(terminateCondition, terminate!)

cb = CallbackSet(quatCallback, microCallback, terminateCallback)

print("Solving ODE\n")
prob = ODEProblem(innerDFun!, x0, tSpan, ODEParams, callback=cb, save_everystep=false)
odeAlg = AutoTsit5(Rosenbrock23())
sol = solve(prob, odeAlg, progress=true)
print("Solution complete\n")
p = plotStates(sol)