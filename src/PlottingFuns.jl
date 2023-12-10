# Functions for plotting various things

function plotAandB(A, B)
    # Plots the A and B matrices in a reasonable form
    stateLabels = ["x", "y", "z", "v_x", "v_y", "v_z", "q0", "q1", "q2", "q3", "ω_x", "ω_y", "ω_z"]
    controlLabels = ["ΔT", "n1", "n2", "n3"]
    Aplot = heatmap(x=stateLabels, y=stateLabels, z=A, colorscale="Picnic",zmin=-1, zmax=1, xaxis="x1", yaxis="y1", hoverongaps=false)
    Bplot = heatmap(x=controlLabels, y=stateLabels, z=B, colorscale="Picnic",zmin=-1, zmax=1, xaxis="x2", yaxis="y2", hoverongaps=false)
    ABplot = [Aplot, Bplot]
    p = plot(ABplot, Layout(yaxis1=attr(scaleanchor="x1",scaleratio=1, autorange="reversed"),
                            yaxis2=attr(scaleanchor="x2",scaleratio=1, autorange="reversed", anchor="x2"),
                            xaxis1_domain=[0, 0.65], xaxis2_domain=[0.8, 1]))
end

function imshow(arr, centered=false)
    # Convenience function to plot a matrix with entries ascending from the upper left
    m, n = size(arr)
    if centered
        p = plot(heatmap(x=1:m, y=1:n, z=arr, colorscale="Picnic", zmid=0), Layout(yaxis_autorange="reversed"))
    else
        p = plot(heatmap(x=1:m, y=1:n, z=arr, colorscale="Picnic"), Layout(yaxis_autorange="reversed"))
    end
    return p
end

function plotEigvals(arr)
    # Plots the eigenvalues of an array, typically corresponding to the zeros of a linear system
    vals = eigvals(arr)
    p = plot(scatter(x=real(vals), y=imag(vals),
             mode="markers",
             marker=attr(size=8, color=1:length(vals), colorscale="Viridis", showscale=false)))
    return p
end

function plotStates(sol)
    # Plots the state variables stored in sol, the solution to an ODE problem
    titles = ["x" "y" "z" "" "v_x" "v_y" "v_z" "" "q2" "q1" "q3" "q0" "ω_y" "ω_x" "ω_z" "" "n2" "n1" "n3" "tt"]

    p_x = scatter(x=sol.t, y=sol[1,:], mode="lines")
    p_y = scatter(x=sol.t, y=sol[2,:], mode="lines")
    p_z = scatter(x=sol.t, y=sol[3,:], mode="lines")

    p_vx = scatter(x=sol.t, y=sol[4,:], mode="lines")
    p_vy = scatter(x=sol.t, y=sol[5,:], mode="lines")
    p_vz = scatter(x=sol.t, y=sol[6,:], mode="lines")

    p_q0 = scatter(x=sol.t, y=sol[7,:], mode="lines")
    p_q1 = scatter(x=sol.t, y=sol[8,:], mode="lines")
    p_q2 = scatter(x=sol.t, y=sol[9,:], mode="lines")
    p_q3 = scatter(x=sol.t, y=sol[10,:], mode="lines")

    p_ωx = scatter(x=sol.t, y=sol[11,:], mode="lines")
    p_ωy = scatter(x=sol.t, y=sol[12,:], mode="lines")
    p_ωz = scatter(x=sol.t, y=sol[13,:], mode="lines")

    p_tt = scatter(x=sol.t, y=sol[14,:])
    p_n1 = scatter(x=sol.t, y=sol[15,:])
    p_n2 = scatter(x=sol.t, y=sol[16,:])
    p_n3 = scatter(x=sol.t, y=sol[17,:])
    

    p1 = make_subplots(rows=5, cols=4, subplot_titles=titles, shared_xaxes=true, shared_yaxes=true)
    # plot positions
    add_trace!(p1, p_x, row=1, col=1)
    add_trace!(p1, p_y, row=1, col=2)
    add_trace!(p1, p_z, row=1, col=3)
    # plot velocities
    add_trace!(p1, p_vx, row=2, col=1)
    add_trace!(p1, p_vy, row=2, col=2)
    add_trace!(p1, p_vz, row=2, col=3)
    # plot rotation quaternion
    add_trace!(p1, p_q2, row=3, col=1)
    add_trace!(p1, p_q1, row=3, col=2)
    add_trace!(p1, p_q3, row=3, col=3)
    add_trace!(p1, p_q0, row=3, col=4)
    # plot angular velocities
    add_trace!(p1, p_ωy, row=4, col=1)
    add_trace!(p1, p_ωx, row=4, col=2)
    add_trace!(p1, p_ωz, row=4, col=3)
    # plot controller inputs
    add_trace!(p1, p_n2, row=5, col=1)
    add_trace!(p1, p_n1, row=5, col=2)
    add_trace!(p1, p_n3, row=5, col=3)
    add_trace!(p1, p_tt, row=5, col=4)

    relayout!(p1, showlegend=false)
    p1
end

