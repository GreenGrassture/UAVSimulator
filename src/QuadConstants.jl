
const g = 9.81 # m/s^2

JxxCF = 16e-6 # kg*m^2
JyyCF = 16e-6 # kg*m^2
JzzCF = 29e-6 # kg*m^2
JMatCF = [JxxCF 0 0
          0 JyyCF 0
          0 0 JzzCF]
massCF =  32e-3 # kg
LCF = 33e-3 # m
kappaCF = 0.01 #
maxMotorForceCF = 0.58 # N

struct QuadParams
    mass::Float64
    J::Matrix{Float64} # moment of inertia matrix
    L::Float64
    kappa::Float64
    maxMotorForce::Float64
    mix::Matrix{Float64}
    unmix::Matrix{Float64}
end


function mixAndUnmix(L, K)
    # mixer matrix converts motor forces into total thrust and three moments
    # unmixer does the reverse, taking the thrust and moments and returning motor forces
    unmixer = [  [ 1  1  1  1]
               L*[ 1 -1 -1  1]
               L*[-1 -1  1  1]
               K*[ 1 -1  1 -1]]
    mixer = inv(unmixer)
    return (mixer, unmixer)
end

mixCF, unmixCF = mixAndUnmix(LCF, kappaCF)

CFConsts = QuadParams(massCF,
                      JMatCF,
                      LCF,
                      kappaCF,
                      maxMotorForceCF,
                      mixCF,
                      unmixCF)


