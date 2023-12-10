using ControlSystems
using MatrixEquations

##### Controllers ##########################
function controllerNone()
    # x: state vector
    # p: ODE parameter dict
    # t: time
    # Just commands motors to do nothing (beyond nominal hover thrust).  Intended to be used as a dummy to get the model of the linear,
    # open-loop system from the nonlinear closed-loop version using AD.
    function nullController(x, p, t)
        mass = p[:quadConsts].mass
        hoverThrust = [mass*g 0 0 0]
        return vec([0 0 0 0] + hoverThrust)
    end
    return nullController
end

function makeLQR(f, xLin, p, tLin, Q, R, tSample=0.01)
    # It takes a nonlinear derivative function and linearizes it around xLin, then designs an optimal controller
    # based on Q and R
    A, B = getAandBMatrices(f, xLin, p, tLin)
    # Our original A and B matrices are for the continuous-time system, so we have to
    # convert to the discrete time equivalent and use that for the next step.
    sys = c2d(ss(A, B, I(size(A)[1]), 0), tSample)
    #return makeLQRk(sys.A, sys.B, Q, R, p, tSample)
    K = tryLQR(sys.A, sys.B, Q, R, lqr)
    function controller(x, p, t)
        mass = p[:quadConsts].mass
        r = p[:reference](x, p, t)
        hoverThrust = vec([mass*g 0 0 0])
        return @views -K*(x[1:13]-r) + hoverThrust
    end
    return controller
end  

function tryLQR(A, B, Q, R, LQRfun)
    # Repeatedly try to apply LQRfun, removing the least controllable state until success or it runs out of states. Fortunately the test for 
    # controllability is the same for continuous and discrete case
    # LQRfun can be either lqr or dlqr depending on whether the controller is intended to be continous or discrete
    # If states were removed, then the K matrix that is returned will have zeros in the rows and columns corresponding to those states
    # The controller will not even attempt to control those state variables.  But it is expected that the uncontrollability
    # will be the result of the quaternion representation being redunant, so this process shouldn't remove more than the one
    # (redunant) state.
    numStates = size(A)[1]
    # Reduced-state matrices
    ARed = A
    BRed = B
    QRed = Q
    RRed = R
    K = nothing
    removedStates = []
    while numStates > 0
        try
            #print(string("Solving ARE with " , numStates, " states...\n"))
            P, evals, evecs = ared(ARed, BRed, RRed, QRed)
            K = inv(RRed + BRed'*P*BRed)*BRed'*P*ARed
            #print("Success\n")
            break
        catch e
            #print("Failure\n")
            removedState, ARed, BRed, QRed, RRed = removeLeastControllableState(ARed, BRed, QRed, RRed)
            print(string("Removed state #", removedState, "\n"))
            numStates = numStates - 1
            push!(removedStates, removedState)
        end
    end
    if length(removedStates) > 0
        K = insertCols(K, removedStates)
    end
    return K
end

function removeLeastControllableState(A, B, Q, R)
# Determines the degree of controllability of each state, and removes the least controllable one
# The eigenvalues of the controllability gramian are the singular values of the controllability
# matrix.
    U, Σ, Vt = svd(ctrb(A, B))
    # The singular values are in decreasing size, and the last one is the least controllable direction
    # The corresponding direction in state space will be the last column of U.
    worstVec = U[:,end]
    # In terms of our original basis, the state that aligns most with the uncontrollable direction
    # is the one that contributes the largest component to our singular vector
    worstState = argmax(abs.(worstVec))
    ARed = deleteRowAndCol(A, worstState)
    BRed = deleteRow(B, worstState)
    QRed = deleteRowAndCol(Q, worstState)
    return (worstState, ARed, BRed, QRed, R)
end

############################################


############################################

function linearize(nonLinDFun, x0, p, t0)
    # Get the jacobian matrix corresponding to nonLinDFun
    # The syntax is a bit weird, but basically we first turn nonLinDFun into an anonymous function of 
    # one variable, then create a new anonymous function of that variable that actually computes the derivative.
    # Then we actually plug in our values of x0 and u0 to get the result at (x0, u0)
    dfun = x -> nonLinDFun(x, p, t0)
    dfdx = x -> jacobian(dfun, x)[1]
    jac = dfdx(x0)
    return jac
end

function getAandBMatrices(nonLinDFun, x0, p, t0)
# Get the A and B matrices for the system ẋ = Ax + Bu.
# This is slightly more specific than linearize() since there are other 
# states in the ODE state vector that shouldn't be in the linear representation
# of the system.
    jac = linearize(nonLinDFun, x0, p, t0)
    A = jac[1:13, 1:13]
    B = jac[1:13, 14:17]
    return (A, B)
end

function ctrbMat(A, B)
    # Computes the controllability for a linear system with ẋ = Ax + Bu

    # Interestingly, there is a connection with the Krylov subspace - the subspace spanned by
    # b, Ab, A^2b etc. for vector b
    # https://en.wikipedia.org/wiki/Krylov_subspace
    M, N = size(A)
    ctrlMat = B
    for m in 1:(M-1)
        ctrlMat = hcat(ctrlMat, A^m*B)
    end
    return ctrlMat
end

function isControllable(A, B)
    M, N = size(A)
    ctrl = ctrbMat(A, B)
    if rank(ctrl) == M
        return true
    else
        return false
    end
end

function uncontrolledStates(A, B)
    # Determine which of the states are not controllable

    # TODO: if the system turns out to be uncontrollable, maybe I should use SVD to transform it into a basis where that
    # uncontrollable subspace affects as few states as possible.  In the standard basis there's no guarantee that the 
    # subspace is orthogonal to the rest of the states but using SVD that should be possible
    M, N = size(A)
    ctrl = ctrbMat(A, B)
    U, Σ, Vt = svd(ctrl)
    # If the rank of the controllability matrix is less than N (the number of states/rows), then the SVD will have some singular values equal to zero.  
    # In that case, we can trace those through the transformation to figure out which states they correspond to.

    # Count nonzero elements of Σ.  They will be in descending order, so we can assume the zero elements are at the end.
    rank = sum(Σ .!= 0) # I'm comparing against exactly zero here, but I could imagine needing to set a tolerance at some point
    if rank == M
        # If the controllability matrix has rank equal to the number of states, then all states are controllable
        return []
    else
        # Collect all of the vectors that correspond to zero-valued singular values, and check each state to see if it
        # can be contstructed from them.  This entails solving a system of linear equations of the form:
        # e_i = a_1v_1 + a_2v_2 + ... a_(M-rank)v_(M-rank) for i in 1:M.  This seems like an overly conservative
        # approach, since it might somehow be the case that the uncontrolled subspace just cuts across the i'th state
        # but I'll have to address that later

        # The singular values represent the degree of controllability along different directions
        # The columns of U are the directions corresponding to each singular value

        # Find all the directions that are not controllable (σ=0):
        zeroVecs = U[:,(rank+1):end]
        states = []
        nVecs = size(zeroVecs)[2]
        for i=1:nVecs
            for state in findall(v_i -> v_i != 0, zeroVecs[:,i])
                push!(states, state)
            end
        end
        return states
    end
end

function microcontroller!(integrator)
    # Encapsulates all the code that might be run on the microcontroller.  Responsible for updating values of controller signal.
    # Will eventually be responsible for running state estimation and motion planning as well.
    p = integrator.p
    x = integrator.u # The state variable is referred to as u by the integrator
    t = integrator.t
    controller = p[:controller]
    uNew = controller(x, p, t)
    # Plug the result from the controller back into the integrator for the next timestep
    integrator.u[14:17] .= uNew
    return nothing
end