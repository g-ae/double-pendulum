using DataFrames, CSV, Plots, LinearAlgebra, Statistics

global m1 = 0.1
global m2 = 0.1
global u = [0.1, 0.1, 0.1, 0.1]

# Import tracked CSV
dfa = DataFrame(CSV.File("mass_a_200.csv"))
dfb = DataFrame(CSV.File("mass_b_200.csv"))

function pendulum_dynamics(u)
    theta1 = u[1]
    theta2 = u[2]
    theta1p = u[3]
    theta2p = u[4]

    delta = theta1 - theta2
    sin_delta, cos_delta = sincos(delta)
    
    M11 = (m1 + m2) * l1
    M12 = m2 * l2 * cos_delta
    R1  = -m2 * l2 * theta2p^2 * sin_delta - (m1 + m2) * g * sin(theta1)
    
    M21 = l1 * cos_delta
    M22 = l2
    R2  = l1 * theta1p^2 * sin_delta - g * sin(theta2)

    # Cramer -> https://fr.wikipedia.org/wiki/R%C3%A8gle_de_Cramer  (éviter de créer trop de tableaux)
    det = M11 * M22 - M12 * M21
    acc1 = (M22 * R1 - M12 * R2) / det
    acc2 = (-M21 * R1 + M11 * R2) / det
    
    return (theta1p, theta2p, acc1, acc2)
end

# Fonction qui calcule tous les points
function rk4()
    # Listes pour stocker l'historique (pour le plot)
    x1 = Float64[]
    y1 = Float64[]
    x2 = Float64[]
    y2 = Float64[]

    k = Float64[]
    p = Float64[]
    m = Float64[]

    println("Calcul en cours...")
    # 4 secondes de simulation
    for i in 1:4000
        global u
        
        # Étape RK4
        k1 = pendulum_dynamics(u)
        k2 = pendulum_dynamics(u .+ (delta_t/2) .* k1)
        k3 = pendulum_dynamics(u .+ (delta_t/2) .* k2)
        k4 = pendulum_dynamics(u .+ delta_t .* k3)
        u = u .+ (delta_t/6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
        
        # Affichage
        theta1 = u[1]
        theta2 = u[2]
        
        push!(x1, l1 * sin(theta1))
        push!(y1, -l1 * cos(theta1))
        push!(x2, x1[end] + l2 * sin(theta2))
        push!(y2, y1[end] - l2 * cos(theta2))
        
        # Calcul énergie pour affichage
        T = 0.5*m1*l1^2*u[3]^2+0.5*m2*(l1^2*u[3]^2+l2^2*u[4]^2+2*l1*l2*u[3]*u[4]*cos(theta1-theta2))
        V = -(m1+m2)*g*l1*cos(theta1)-m2*g*l2*cos(theta2)
        push!(k, T)
        push!(p, V)
        push!(m, T + V)
    end
    return ([x1,y1,x2,y2], [k,p,m])
end

function save_mp4(position_data)
    x1, y1, x2, y2 = position_data

    frame_finale = length(dfa.x)
    VIDEO_DELTA_T = 0.01
    CALC_DELTA = trunc(Int, VIDEO_DELTA_T / delta_t)

    HISTORY = 400

    # Animation
    function plot_at(i::Int)
        Plots.plot(x2[max(1, i-HISTORY):i], y2[max(1, i-HISTORY):i], aspect_ratio = :equal, size=(600, 600), xlims=(-(l1+l2), (l1+l2)), ylims=(-(l1+l2), (l1+l2)), legend=false, axis=([], false))
                
        Plots.plot!([0, x1[i]], [0, y1[i]], color=:black, linewidth=2)
        Plots.plot!([x1[i], x2[i]], [y1[i], y2[i]], color=:black, linewidth=2)

        Plots.scatter!([x1[i]], [y1[i]], color=:red, markersize=5)
        Plots.scatter!([x2[i]], [y2[i]], color=:red, markersize=5)

        ir = div(i, CALC_DELTA) + 1
        
        # On ne dépasse pas la taille du fichier (frame_finale)
        if ir <= frame_finale
            ir_max = clamp(ir, 1, frame_finale)
            
            # On affiche tout l'historique jusqu'à l'instant t
            Plots.plot!(dfb.x[trunc(Int,max(1, ir_max - HISTORY / CALC_DELTA)):ir_max], dfb.y[trunc(Int,max(1,ir_max - HISTORY / CALC_DELTA)):ir_max], color=:blue, linewidth=2, label="Tracked B", alpha = 0.9)
            
            if ir < frame_finale
                Plots.scatter!(dfa.x[ir_max:ir_max], dfa.y[ir_max:ir_max], color=:green, markersize=3, alpha= 0.9)
                Plots.scatter!(dfb.x[ir_max:ir_max], dfb.y[ir_max:ir_max], color=:blue, markersize=3, alpha = 0.9)
            end
        end
    end
    println("Génération de l'animation...")
    anim = @animate for i in 1:length(x1)
        if i % 500 == 0
            println("Progression: ", i, "/", length(x1))
        end
        plot_at(i)
    end every 10
    gif(anim, "double_pendulum.mp4", fps = 100)
end

function save_position_images(position_data, starting_u)
    x1,y1,x2,y2 = position_data

    println("Sauvegarde des images de position...")

    frame_finale = length(dfa.x)
    VIDEO_DELTA_T = 0.01
    CALC_DELTA = trunc(Int, VIDEO_DELTA_T / delta_t)

    x1_affichage = x1[1:CALC_DELTA:frame_finale*CALC_DELTA]
    y1_affichage = y1[1:CALC_DELTA:frame_finale*CALC_DELTA]

    # Preparation des deux
    dfa_x_affichage = dfa.x[1:frame_finale]
    dfa_y_affichage = dfa.y[1:frame_finale]

    # Plot de position masse A
    m1r = round(m1, digits=5)
    m2r = round(m2, digits=5)
    Plots.plot(1:frame_finale, x1_affichage, label="Calculated", xlabel="frame", ylabel="position X", title="MASSE A, m1: $m1r kg, m2: $m2r kg")
    Plots.plot!(1:frame_finale, dfa_x_affichage, label="Tracked A")
    savefig("masse_a_x.png")

    ur = round.(starting_u, digits=5)
    Plots.plot(1:frame_finale, y1_affichage, label="Calculated", xlabel="frame", ylabel="position Y", title="MASSE A, $ur")
    Plots.plot!(1:frame_finale, dfa_y_affichage, label="Tracked A")

    savefig("masse_a_y.png")

    # Plot de position masse B
    Plots.plot(1:frame_finale, x2[1:CALC_DELTA:frame_finale*CALC_DELTA + 1 - CALC_DELTA], label="Calculated", xlabel="frame", ylabel="position X", title="MASSE B, m1: $m1r kg, m2: $m2r kg")
    Plots.plot!(1:frame_finale, dfb.x[1:frame_finale], label="Tracked B")
    savefig("masse_b_x.png")

    Plots.plot(1:frame_finale, y2[1:CALC_DELTA:frame_finale*CALC_DELTA + 1 - CALC_DELTA], label="Calculated", xlabel="frame", ylabel="position Y", title="MASSE B, $ur")
    Plots.plot!(1:frame_finale, dfb.y[1:frame_finale], label="Tracked B")
    savefig("masse_b_y.png")
end

function calculate_RMSE(position_data)
    x1, y1, x2, y2 = position_data

    # These constants are defined in the calling script (optimisation_pendule.jl)
    # If running this file standalone, they would need to be defined here.
    frame_finale = length(dfa.x)
    VIDEO_DELTA_T = 0.01
    CALC_DELTA = trunc(Int, VIDEO_DELTA_T / delta_t)

    # Downsample calculated data to match tracked data's framerate
    num_points = min(frame_finale, floor(Int, length(x1) / CALC_DELTA))
    
    indices = 1:CALC_DELTA:(num_points * CALC_DELTA)
    
    x1_calc = x1[indices]
    y1_calc = y1[indices]
    x2_calc = x2[indices]
    y2_calc = y2[indices]

    # Create arrays of vectors for comparison
    traj_calc_a = [[x, y] for (x, y) in zip(x1_calc, y1_calc)]
    traj_track_a = [[x, y] for (x, y) in zip(dfa.x[1:num_points], dfa.y[1:num_points])]
    
    traj_calc_b = [[x, y] for (x, y) in zip(x2_calc, y2_calc)]
    traj_track_b = [[x, y] for (x, y) in zip(dfb.x[1:num_points], dfb.y[1:num_points])]

    # Internal RMSE helper
    _rmse(v1, v2) = sqrt(mean(norm.(v1 .- v2) .^ 2))

    rmse_a = _rmse(traj_calc_a, traj_track_a)
    rmse_b = _rmse(traj_calc_b, traj_track_b)

    println("RMSE Masse A: ", round(rmse_a, digits=5))
    println("RMSE Masse B: ", round(rmse_b, digits=5))

    return rmse_a + rmse_b
end

function calc_energie(energie_data)
    k,p,m = energie_data

    println("Génération de l'image des énergies...")
    Plots.plot(1:iterations, k, label="Ec", ylabel="Energie", xlabel="Iterations ($(round(delta_t * iterations))/frame)")
    Plots.plot!(1:iterations, p, label="Ep")
    Plots.plot!(1:iterations, m, label="Em")
    
    savefig("cinetique.png")
    
    println("Énergie totale au point 1: ", m[1])
    println("Énergie totale à la fin: ", m[end])
    println("Équivaut à une perte de ", 100 - m[end] / m[1] * 100, " %")
end

function load_double_pendulum_mp4(masse1, masse2, t1, t2, t1p, t2p)
    println("\n-- Loading double pendulum data to mp4 --")
    global m1 = masse1
    global m2 = masse2
    starting_u = [t1,t2,t1p,t2p]
    global u = [t1,t2,t1p,t2p]
    pos_data, energie_data = rk4()
    save_mp4(pos_data)
    save_position_images(pos_data, starting_u)
    #calc_energie(energie_data)
    println("Calculated RMSE : ", calculate_RMSE(pos_data))
end

if abspath(PROGRAM_FILE) == @__FILE__
    load_double_pendulum_mp4(0.0306, 0.003592, 3.1728, 3.225, 0.0, 0.0)
end