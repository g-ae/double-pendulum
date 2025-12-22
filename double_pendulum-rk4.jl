using Plots, Statistics, LinearAlgebra

# Constants
const l1 = 0.09174
const l2 = 0.06933
const m1 = 0.020283
const m2 = 0.002083
const g = 9.81
const delta_t = 0.001
const iterations = 4000

const skip_animation = true

# const m1 = 0.041
# const m2 = 0.006

# x=2.084 mm, y=91.66 mm

# État du système : [theta1, theta2, theta1p, theta2p]
# Passage dans fonctions facilité
u = [pi + 0.02274, pi + 0.0951247, 0.0, 0.0]
# u = [pi/2 *3, pi/2 * 3, 0.0, 0.0]

# Cette fonction renvoie [d_theta1, d_theta2, d_theta1p, d_theta2p]
function pendulum_dynamics(u)
    theta1 = u[1]
    theta2 = u[2]
    theta1p = u[3]
    theta2p = u[4]
    
    delta = theta1 - theta2

    # résolution pour theta1pp et theta2pp
    M11 = (m1 + m2) * l1
    M12 = m2 * l2 * cos(delta)
    R1  = -m2 * l2 * theta2p^2 * sin(delta) - (m1 + m2) * g * sin(theta1)
    
    M21 = l1 * cos(delta)
    M22 = l2
    R2  = l1 * theta1p^2 * sin(delta) - g * sin(theta2)
    
    Matrice = [M11 M12; M21 M22]
    Resultats = [R1, R2]
    
    accels = Matrice \ Resultats # [theta1pp, theta2pp]

    # On retourne les dérivées :
    # La dérivée de la position est la vitesse
    # La dérivée de la vitesse est l'accélération
    return [theta1p, theta2p, accels[1], accels[2]]
end

# Listes pour stocker l'historique (pour le plot)
x1 = Float64[]
y1 = Float64[]
x2 = Float64[]
y2 = Float64[]

k = Float64[]
p = Float64[]
m = Float64[]

println("Calcul en cours...")
for i in 1:iterations
    global u
    
    # Étape RK4
    k1 = pendulum_dynamics(u)
    k2 = pendulum_dynamics(u + delta_t*k1/2)
    k3 = pendulum_dynamics(u + delta_t*k2/2)
    k4 = pendulum_dynamics(u + delta_t*k3)
    u += delta_t * (k1 + 2*k2 + 2*k3 + k4) / 6
    
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

function plot_at(i::Int)
    Plots.plot(x2[1:i], y2[1:i], aspect_ratio = :equal, size=(600, 600), xlims=(-(l1+l2), (l1+l2)), ylims=(-(l1+l2), (l1+l2)), legend=false, axis=([], false))
               
    Plots.plot!([0, x1[i]], [0, y1[i]], color=:black, linewidth=2)
    Plots.plot!([x1[i], x2[i]], [y1[i], y2[i]], color=:black, linewidth=2)
    
    Plots.scatter!([x1[i]], [y1[i]], color=:red, markersize=5)
    Plots.scatter!([x2[i]], [y2[i]], color=:red, markersize=5)
end

# Animation
if !skip_animation
    println("Génération de l'animation...")
    anim = @animate for i in 1:length(x1)
        if i % 100 == 0
            println("Progression: ", i, "/", length(x1))
        end
        plot_at(i)
    end every 10
    gif(anim, "double_pendulum-rk4.mp4", fps = 100)
end


# Energies
function calc_energie()
    println("Génération de l'image des énergies...")
    Plots.plot(1:iterations, k, label="Ec", ylabel="Energie", xlabel="Iterations ($(round(delta_t * iterations))/frame)")
    Plots.plot!(1:iterations, p, label="Ep")
    Plots.plot!(1:iterations, m, label="Em")
    
    savefig("cinetique.png")
    
    println("Énergie totale au point 1: ", m[1])
    println("Énergie totale à la fin: ", m[end])
    println("Équivaut à une perte de ", 100 - m[end] / m[1] * 100, " %")
end

# Two arrays of same length
function RMSE(arrVect1, arrVect2)
    if length(arrVect1) != length(arrVect2)
        throw("error, not same length")
    end
    return sqrt(mean(norm.(arrVect1 .- arrVect2).^2))
end

# nombre de frames -> 70 frames réelles
using DataFrames, CSV
const frame_finale = 70
const VIDEO_DELTA_T = 0.01
const CALC_DELTA = trunc(Int, VIDEO_DELTA_T / delta_t)
const DIFF = 10*CALC_DELTA

# Import CSV tracked positions
dfa = DataFrame(CSV.File("mass_a.csv"))
dfb = DataFrame(CSV.File("mass_b.csv"))

# Arrays de positions préparées
x1_affichage = x1[DIFF:CALC_DELTA:frame_finale*CALC_DELTA + DIFF-CALC_DELTA]
y1_affichage = y1[DIFF:CALC_DELTA:frame_finale*CALC_DELTA + DIFF-CALC_DELTA]

# Changer ce bool si on veut/pas les graphiques des positions
if true
    # Preparation des deux
    dfa_x_affichage = dfa.x[1:frame_finale]
    dfa_y_affichage = dfa.y[1:frame_finale]

    # X
    Plots.plot(1:frame_finale, x1_affichage, label="Calculated", xlabel="frame", ylabel="position X", title="MASSE A, m1: $m1, m2: $m2")
    Plots.plot!(1:frame_finale, dfa_x_affichage, label="Tracked A")

    savefig("masse_a_x.png")

    # Y

    Plots.plot(1:frame_finale, y1_affichage, label="Calculated", xlabel="frame", ylabel="position Y", title="MASSE A, m1: $m1, m2: $m2")
    Plots.plot!(1:frame_finale, dfa_y_affichage, label="Tracked A")

    savefig("masse_a_y.png")

    # Plot de position masse B

    Plots.plot(1:frame_finale, x2[DIFF:CALC_DELTA:frame_finale*CALC_DELTA + DIFF-CALC_DELTA], label="Calculated", xlabel="frame", ylabel="position X", title="MASSE B, m1: $m1, m2: $m2")
    Plots.plot!(1:frame_finale, dfb.x[1:frame_finale], label="Tracked B")
    savefig("masse_b_x.png")

    Plots.plot(1:frame_finale, y2[DIFF:CALC_DELTA:frame_finale*CALC_DELTA + DIFF-CALC_DELTA], label="Calculated", xlabel="frame", ylabel="position Y", title="MASSE B, m1: $m1, m2: $m2")
    Plots.plot!(1:frame_finale, dfb.y[1:frame_finale], label="Tracked B")
    savefig("masse_b_y.png")
end
# Masse A
# XY
traj_calc_a = [[x,y] for (x,y) in zip(x1_affichage, y1_affichage)]
traj_track_a = [[x, y] for (x, y) in zip(dfa.x[1:frame_finale], dfa.y[1:frame_finale])]
rmse_a = RMSE(traj_calc_a, traj_track_a)
println("RMSE masse A ($m1) : ", rmse_a)

# Masse B
x2_affichage = x2[DIFF:CALC_DELTA:frame_finale*CALC_DELTA + DIFF-CALC_DELTA]
y2_affichage = y2[DIFF:CALC_DELTA:frame_finale*CALC_DELTA + DIFF-CALC_DELTA]

traj_calc_b = [[x,y] for (x,y) in zip(x2_affichage, y2_affichage)]
traj_track_b = [[x, y] for (x, y) in zip(dfb.x[1:frame_finale], dfb.y[1:frame_finale])]
rmse_b = RMSE(traj_calc_b, traj_track_b)
println("RMSE masse B ($m2) : ", rmse_b)

# Total
println("RMSE Total (A + B) : ", rmse_a + rmse_b)