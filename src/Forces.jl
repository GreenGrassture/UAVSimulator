using LinearAlgebra

##### External Forces ######################
function externalForcesNone(x, t, quadConsts)
    # Function that gives externally imposed forces as a function of time (in the earth frame)
    # This one applies no forces.
    zeroForce = [0.0, 0.0, 0.0]
    zeroMoment = [0.0, 0.0, 0.0]
    return vcat(zeroForce, zeroMoment)
end
function externalForcesOnlyG(x, p, t)
    # Function that gives externally imposed forces as a function of time (in the earth frame)
    # The usual situation, with no forces other than gravity
    mass = p[:quadConsts].mass
    #J = p["quadConsts"].J
    gravForce_E = [0.0, 0.0, -mass*g, ] # Earth Frame.  Gravity will always be vertical and downward in the earth frame.
    zeroMoment_E = [0.0, 0.0, 0.0]
    return vcat(gravForce_E, zeroMoment_E)
end

function createForceFunctionGravAndNoise(noiseMeans, noiseCovarianceMatrix, bandwidth, numInterp)
    #noiseDiscrete = makeNormal(noiseMeans, noiseCovarianceMatrix)
    #numInterp = 11 # the number of samples that will be included in the interpolation
    #bandwidth = 500 # Hz
    #noiseContinuous = makeInterpolatedNoiseFunction(noiseDiscrete, bandwidth, numInterp)
    noiseContinuous = makeContinuousNormal(noiseMeans, noiseCovarianceMatrix, bandwidth, numInterp)


    function externalForcesGravAndNoise(x, p, t)
        mass = p[:quadConsts].mass # theoretically this could be pulled out of the inner function right now, but I want to allow for non-constant mass
        #J = p["quadConsts"].J
        gravForce_E = [0.0, 0.0, -mass*g]
        n = noiseContinuous(t)
        forces_E = @views gravForce_E + n[1:3]
        moments_E = @views n[4:end]
        return vcat(forces_E, moments_E)
    end
    return externalForcesGravAndNoise
end

tSample = 0.05
means = zeros(Float64, 6)
covs = Diagonal(0.001*ones(Float64, 6))
externalForcesGravAndNoise = createForceFunctionGravAndNoise(means, covs, 1/tSample, 5)

        
############################################
