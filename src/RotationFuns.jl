# For reference see https://rotations.berkeley.edu/other-representations-of-a-rotation/

function skewMat(vec)
    # Promotes a vector in R^3 to a skew-symmetric matrix

    #M = [0 -vec[3] vec[2]
    #     vec[3] 0 -vec[1]
    #    -vec[2] vec[1] 0]
    M = reshape([0, vec[3], -vec[2],
                 -vec[3], 0, vec[1],
                 vec[2], -vec[1], 0], (3,3))
    return M
end



function quat2mat(q)
    # Turns a quaternion, represented as a length 4 array, into a rotation matrix
    # q = [q0, q1, q2, q3]
    # Careful with the subscripts/indices here.  q0 = q[1]!!
    q0, q1, q2, q3 = q
    #I3 = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    M1 = (2.0*q0^2 - 1.0)*I(3)
    @views M2 = 2.0 .* q[2:4]*q[2:4]'
    @views M3 = 2.0*q0*skewMat(q[2:4])
    return M1 + M2 + M3
end

function quat2euler(q)
    q0, q1, q2, q3 = q
    Φ = atan(2*(q0*q1 + q2*q3), q0^2 - q1^2 - q2^2 + q3^2)
    θ = asin(2*(q0*q2 - q3*q1))
    ψ = atan(2*(q0*q3 + q1*q2), q0^2 + q1^2 - q2^2 - q3^2)
    return rad2deg.([Φ, θ, ψ])
end

function euler2quat(Φ, θ, ψ)
    # ψ: yaw
    # θ: pitch
    # Φ: roll
    Φ, θ, ψ = deg2rad.([Φ, θ, ψ])
    q0 = cos(Φ/2)*cos(θ/2)*cos(ψ/2) + sin(Φ/2)*sin(θ/2)*sin(ψ/2)
    q1 = sin(Φ/2)*cos(θ/2)*cos(ψ/2) - cos(Φ/2)*sin(θ/2)*sin(ψ/2)
    q2 = cos(Φ/2)*sin(θ/2)*cos(ψ/2) + sin(Φ/2)*cos(θ/2)*sin(ψ/2)
    q3 = cos(Φ/2)*cos(θ/2)*sin(ψ/2) - sin(Φ/2)*sin(θ/2)*cos(ψ/2)
    return [q0, q1, q2, q3]
end

function euler2quat(eulerAngles)
    Φ, θ, ψ = eulerAngles
    return euler2quat(Φ, θ, ψ)
end
