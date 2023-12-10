##### Reference Signals ####################
# Return the reference signal as a function of system state, simulation parameters, and time
# This conforms to the convention used by the integrator.  Most reference signals won't
# actually care about most of those parameters.

function referenceSignalConstant(x, p, t)
    # This just returns a constant value, which is the equilibrium position for the quadcopter
    sRef = vec([0, 0, 0])
    vRef = vec([0, 0, 0])
    qRef = euler2quat(0, 0, 0)
    omegaRef = vec([0, 0, 0])
    xRef = [sRef; vRef; qRef; omegaRef]
    return xRef
end

function makeReference_Constant(x0, v0, q0, ω0)
# Builds a constant reference signal with the specified conditions.
    function referenceSignal(x, p, t)
        return [x0; v0; q0; ω0]
    end
    return referenceSignal
end

function makeReference_Ramp(s0, sf, t0, tf)
# Makes a reference signal that ramps linearly from position s0 at time t0
# to position sf at time tf
    function referenceSignal(x, p, t)
        # Returns the reference signal as a function of time, and possibly the system state.
        # This one ramps the desired position from s0 at time t0 to sF at time tf and holds it there
        if t < t0
            sRef = s0
        elseif t < tf
            # linearly interpolate between initial and final positions
            sRef = (tf - t)/(tf - t0)*s0 + (t - t0)/(tf - t0)*sf
        else
            sRef = sf
        end
        vRef = vec([0, 0, 0])
        qRef = euler2quat(0, 0, 0)
        ωRef = vec([0, 0, 0])
        xRef = [sRef; vRef; qRef; ωRef]
    end
    return referenceSignal
end

function makeReference_Sine(s0, sf, freq)
# Creates a position reference signal that oscillates between position s0 and sf with 
# frequency freq
    function referenceSignal(x, p, t)
        sRef = vec(s0 + sf.*sin(2*pi*freq*t))
        vRef = [0.0, 0.0, 0.0]
        qRef = euler2quat(0, 0, 0)
        ωRef = [0.0, 0.0, 0.0]
        xRef = [sRef; vRef; qRef; ωRef]
    end
    return referenceSignal
end

rampRef = makeReference_Ramp([0, 0, 0], [1, 0, 0], 0, 1)