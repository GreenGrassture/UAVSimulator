# Functions for computing the behavior of the system

function autodiffDFun(x, p, t)
# Wrapper function for the in-place derivative function to perform automatic differentiation
# Zygote doesn't allow array mutation, so we have to wrap the array in a Buffer
    dx = Zygote.Buffer(x)
    innerDFun!(dx, x, p, t)
    return copy(dx)
end

function innerDFun!(dx, x, p, t)
    # Has x and u as arguments to allow for taking derivatives wrt to each of them.  
    # Wrapped in an outer function to make it compatable with the ODE solver
    #
    # x: the current state vector.  Includes the output of the controller
    # p: parameters of the integrator.  Used for holding things like mass, moment of inertia, the controller etc.
    # State: [x1,  x2,  x3,  v1,  v2,  v3,  q0,  q1,  q2,  q3,  p,  q,  r,  tt,  nx,  ny,  nz]
    #        [1]             [4]            [7]                [11]        [14]           [17]
    #        position        velocity       orientation        angular vel  controller outputs

    quadConsts = p[:quadConsts]
    mass = quadConsts.mass
    J = quadConsts.J

    #s = @views x[1:3] # position
    v = @views x[4:6] # velocity
    q = @views x[7:10] # quaternion representation of rotation
    ω = @views x[11:13] # angular velocity
    u = @views x[14:17]

    # TODO motor dynamics can in principle be functions of the entire state.  Update function arguments to fit
    uLim = motorModel(u, p) # Limits the motor commands to the physically possible values
    tt = uLim[1] # Total thrust
    n = @views uLim[2:4] # Moments about x, y, z axes (defined in the body frame)

    T_BE = quat2mat(q) # rotation matrix form of rotation q
    T_EB = T_BE' # T_EB is orthogonal so its inverse is equal to its transpose

    # Compute accelerations, transforming them all into the earth frame E before combining
    motorForce_B = [0, 0, tt] # Body frame.  Thrust will always be along the vertical axis of the quadcopter in the body frame.
    motorForce_E = T_EB*motorForce_B # Converted to earth frame
    externalForcesAndMoments_E = p[:forces](x, p, t)
    @views externalForces_E = externalForcesAndMoments_E[1:3]
    @views externalMoments_E = externalForcesAndMoments_E[3:end]
    
    #ds = sDot(v)
    dx[1:3] = sDot(v)

    #dv = vDot((motorForce_E, externalForces_E), mass)
    dx[4:6] = vDot((motorForce_E, externalForces_E), mass)

    #dq = qDot(q, ω)
    dx[7:10] = qDot(q, ω)

    #dω = ωDot(ω, n, J)
    dx[11:13] = ωDot(ω, n, J)

    # du - these states are only changed by the microcontroller callback
    dx[14:17] = zeros(Float64, 4)

    #return [ds; dv; dq; dω]
    return nothing
end

@inline function qDot(q, ω)
    M = @views 0.5*[-q[2:4]'; q[1]*I(3)+skewMat(q[2:4])]
    dq = M*ω
    return dq
end

@inline function ωDot(ω, n, J)
    return inv(J)*(-skewMat(ω)*J*ω + n)
end

@inline function vDot(forces_E, mass)
    return sum(forces_E)/mass
end

@inline function sDot(v)
    return v
end

function quaternionError(resid, x, params, t)
# Function that measures the error in the norm of the rotation quaternion.  Used according to:
# https://diffeq.sciml.ai/stable/features/callback_library/#Manifold-Conservation-and-Projection
    #@views quatNorm = dot(x[7:10], x[7:10])
    quatNorm =  x[7]^2 + x[8]^2 + x[9]^2 + x[10]^2
    resid[1] = 1 - quatNorm
    resid[2:17] = zeros(16)
end

function motorModel(u, p)
    fMax = p[:quadConsts].maxMotorForce
    motorForces = p[:quadConsts].unmix*u
    # Limit the force that the motors can produce.  Really there should be a maximum motor torque, and the way that translates to force will be
    # different depending on whether the propeller is moving in the right direction.

    # should replace this with an implementation using clamp()
    for i in 1:size(motorForces)[1]
        force = motorForces[i]
        if abs(force) > fMax
            motorForces[i] = sign(force)*fMax
        end
    end
    @views uLim = p[:quadConsts].mix*motorForces
    return uLim
end

function speed2force(ω)
    # From values found in 136 lab
    return 4.254E-05*ω - 2.149E-02
end

function speed2PWM(ωDesired)
    # Convert a desired motor speed in radians/second into the necessary pwm command to achieve that speed 
    # (according to our linear model)
    a = -39.12526  # the zeroth order term
    b = 0.096503   # the first order term
    return Int64(round(a + b*ωDesired));
end

function force2speed(desiredForce)
    # Convert a desired force (in a single propeller) into the motor speed in radians/second that will produce that force
    propConstant = 1.0e-08
    # We implement a safety check,
    # (no sqrtf for negative numbers)
    if (desiredForce <= 0) 
        return 0.0
    else
        return sqrt(desiredForce/propConstant)
    end
end