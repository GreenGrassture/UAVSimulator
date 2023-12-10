using Statistics
using Memoization
using LRUCache
using Distributions
using LinearAlgebra
using StaticArrays
using Random

include("noiseFunctionTypes.jl")

function sampleFrom(lcg, t::Float64)
# Implement a linear congruential generator
# Uses one of the constants reported in Steele and Vigna 2021
# The function is completely deterministic, so a particular instance
# of the function will always return the same value for a given time
# There are two modifications made - they look similar but are for different reasons.
# The initial xor with the private seed is to suppress correlations in the outputs
# that survived the LCG randomization.  The second xor is to make sure that
# numbers that enter the same point in the LCG period nevertheless exit 
# in different places.
    y0 = reinterpret(UInt64, t)

    # Improved version of MurmerHash3 tested
    y0 = xor(y0, (y0 >> 33));
    y0 *= 0xff51afd7ed558ccd;
    y0 = xor(y0, (y0 >> 33));
    y0 *= 0xc4ceb9fe1a85ec53;
    y0 = xor(y0, (y0 >> 33));

    y0 = xor(y0, lcg.privateSeed)
    y = UInt128(y0)
    for i in 1:10
        #y = mod(lcg.A*y, M)
        y = (A64*y)&(M-1)
    end
    y2 = UInt64(y)
    #y2 = xor(bitrotate(y2, 31), y2)
    y2 = xor(y2, lcg.privateSeed)
    u = Float64(y2/M)
    return u
end

#const A64 = UInt64(0xfc0072fa0b15f4fd)
#const M = UInt128(2)^64
lcg = LCG()

function uniformToStandardNormal(u)
    # Transforms a uniform distribution on [0, 1] into a standard normal
    return quantile(Normal(0,1), u)
end

function makeNormal(meanVec, covMat)
    lcgs = [LCG() for i in 1:length(meanVec)]
    return makeNormal(meanVec, covMat, lcgs)
end

function makeNormal(meanVec, covMat, lcgs)
# Constructs a function to sample from a potentially multivariate normal distribution with the given
# mean and covariance matrix
    #lcgs = [makeLCG(10, A64, M) for i in 1:length(meanVec)]
    #lcgs = [LCG() for i in 1:length(meanVec)]
    A = cholesky(covMat).L # AA* = Σ because K_xx ≡ E[x*x.T] = E[Az(Az).T] = A*E[z*z.T]A.T = A*I*A.T ≡ Σ

    function sampleFromUnivariateNormal(t)
    #@memoize LRU(maxsize=21) function sampleFromUnivariateNormal(t)
        z = uniformToStandardNormal(sampleFrom(lcgs[1],t))
        return covMat*z + meanVec
    end

    function sampleFromMultivariateNormal(t)
    #@memoize LRU{Float64, Array{Float64}}(maxsize=2000) function sampleFromMultivariateNormal(t)
        us = [sampleFrom(lcg,t) for lcg in lcgs]
        zs = uniformToStandardNormal.(us)
        xs = A*zs + meanVec
        return xs
    end
    if length(meanVec) == 1
        return sampleFromUnivariateNormal
    else
        return sampleFromMultivariateNormal
    end
end

function makeDiscreteNormal(meanVec::Vector{<:Number}, covMat::Matrix{<:Number}, lcgs::Array{LCG})
# Non-standard 
# multivariate
# RNG provided
    # Intended for the case where discrete time noise is the final intended result (e.g. sensor noise)
    A = cholesky(covMat).L # AA* = Σ because K_xx ≡ E[x*x.T] = E[Az(Az).T] = A*E[z*z.T]A.T = A*I*A.T ≡ Σ
    function sampleFromDiscreteNormal(t)
        us = [sampleFrom(lcg,t) for lcg in lcgs]
        zs = uniformToStandardNormal.(us)
        xs = A*zs + meanVec
        return xs
    end
    return sampleFromDiscreteNormal 
end

function makeDiscreteNormal(meanVec::Vector{<:Number}, covMat::Matrix{<:Number})
# Non-standard
# multivariate 
# RNG not provided
    lcgs = [LCG() for i in 1:length(meanVec)]
    return makeDiscreteNormal(meanVec, covMat, lcgs)
end

function makeDiscreteNormal(meanVal::Number, varVal::Number, lcg::LCG)
# Non-standard 
# univariate 
# RNG provided
    function sampleFromDiscreteNormal(t)
        z = uniformToStandardNormal(sampleFrom(lcg,t))
        return varVal*z + meanVal
    end
end

function makeDiscreteNormal(meanVal::Number, varVal::Number)
# Non-standard 
# univariate 
# RNG not provided
    lcg = LCG()
    return makeDiscreteNormal(meanVal, varVal, lcg)
end

function makeDiscreteNormal(lcgs::Array{LCG})
# Standard 
# multivariate 
# RNG provided
    # Intended for the case where discrete time noise is an intermediate result, to be used inside of makeContinuousNormal
    function sampleFromDiscreteNormal(t)
        us = [sampleFrom(lcg,t) for lcg in lcgs]
        zs = uniformToStandardNormal.(us)
        return zs
    end
    return sampleFromDiscreteNormal
end
function makeDiscreteNormal(lcg::LCG)
# Standard
# univariate 
# RNG provided
    function sampleFromDiscreteNormal(t)
        return uniformToStandardNormal(sampleFrom(lcg, t))
    end
    return sampleFromDiscreteNormal
end

function makeDiscreteNormal(dims::Int)
# Standard
# Univariate OR multivariate
# RNG not provided
    # Creates a function to sample from the standard normal of any dimension
    # I expect that this will be type unstable, but it should only run a couple times
    # and the function it returns should itself be type stable
    if dims < 1
        throw(DomainError(nothing))
    elseif dims == 1
        return makeDiscreteNormal(LCG())
    else
        return makeDiscreteNormal([LCG() for i in 1:dims])
    end
end

function makeContinuousNormal(meanVec, covMat, lcgs, bw, nSamples)
# RNG provided
    A = cholesky(covMat).L # AA* = Σ because K_xx ≡ E[x*x.T] = E[Az(Az).T] = A*E[z*z.T]A.T = A*I*A.T ≡ Σ
    standardContNormal = makeInterpolatedNoiseFunction(makeDiscreteNormal(lcgs), bw, nSamples)
    function sampleFromContinuousNormal(t)
        zs = standardContNormal(t)
        return A*zs + meanVec
    end
    return sampleFromContinuousNormal
end

function makeContinuousNormal(meanVec, covMat, bw, nSamples)
# RNG not provided
    lcgs = [LCG() for i in 1:length(meanVec)]
    return makeContinuousNormal(meanVec, covMat, lcgs, bw, nSamples)
end

function makeInterpolatedNoiseFunction(discreteSampler, bw, nSamples::Int)
# Constructs a function to sample from a continuous-time function with
# statistical properties defined by discreteSampler, with maximum bandwidth bw 
    dtSample = 1/(2*bw)
    # The case with no smoothing at all.  Just return the closest sample.
    if nSamples == 1
        function noiseSingleSample(t)
            tInTermsOfK = t/dtSample
            kClosest = Int(floor(tInTermsOfK))
            return discreteSampler(dtSample*kClosest)
        end
        return noiseSingleSample
    else
        mid = Int(floor(nSamples/2))+1
        minIdx = mid-nSamples
        maxIdx = nSamples-mid
        idxs = [j for j in minIdx:maxIdx]
        function noiseMultipleSamples(t)
            # There is some unnecessary conversion back and forth between integers and floats here
            # as a result of the rng taking a float (t) as an input.  Really we could rewrite the 
            # rng to accept an int and save some trouble.
            tInTermsOfK = t/dtSample
            kClosest = Int(floor(tInTermsOfK))
            offset = tInTermsOfK - kClosest
            ks = idxs .+ kClosest
            # The input times for the samples need to be exactly the same every time
            # in order for the cache to work, so I'm using a comprehension here to make
            # sure nothing weird happens with range arithmetic
            tsSampler = [dtSample*k for k in ks] # The times at which the samples are drawn from the rng
            # Using a list comprehension here so that the arguments to discreteSampler
            # are scalars, and we can cache them individually.  There will be overlap
            # between calls, but not entirely, and if we call it with a vector there
            # won't be any overlap from the point of view of the cache
            samples = [discreteSampler(tS) for tS in tsSampler]
            tsKernel = Float64.(idxs) .- offset
            #kernelVals = sinc.(tsKernel).*sinc.(tsKernel./5)
            kernelVals = sinc.(tsKernel)
            return sum(samples.*kernelVals)
        end
        return noiseMultipleSamples
    end
end

meanVec = [1.0, 2.0, 3.0]
covMat =  [1.0 0.5 0.2;
           0.5 2.0 0.0;
           0.2 0.0 3.0]
lcgs3 = [LCG(), LCG(), LCG()]
myNormalDisc = makeNormal(meanVec, covMat, lcgs3)
bw = 50
nSamples = 5
myNormalCont1 = makeInterpolatedNoiseFunction(myNormalDisc, bw, nSamples)
myNormalCont2 = makeContinuousNormal(meanVec, covMat, lcgs3, bw, nSamples)
ts = [t for t in -10.0:0.001:10.0]
samplesOld = reduce(hcat, myNormalCont1.(ts));
samplesNew = reduce(hcat, myNormalCont2.(ts));
