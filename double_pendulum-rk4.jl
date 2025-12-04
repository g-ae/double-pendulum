using Plots

# Constants
const l1 = 5.0
const l2 = 4.0
const m1 = 1.0
const m2 = 1.0
const g = 9.81
const delta_t = 0.005
const iterations = 20000

# État du système : [theta1, theta2, theta1p, theta2p]
# Passage dans fonctions facilité
u = [pi, pi + 0.2, 0.0, 0.0]

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

println("Calcul en cours...")
for i in 1:iterations
    global u
    
    theta1 = u[1]
    theta2 = u[2]
    
    push!(x1, l1 * sin(theta1))
    push!(y1, -l1 * cos(theta1))
    push!(x2, x1[end] + l2 * sin(theta2))
    push!(y2, y1[end] - l2 * cos(theta2))
    
    # Étape RK4
    k1 = pendulum_dynamics(u)
    k2 = pendulum_dynamics(u + delta_t*k1/2)
    k3 = pendulum_dynamics(u + delta_t*k2/2)
    k4 = pendulum_dynamics(u + delta_t*k3)
    u += (delta_t/6) * (k1 + 2*k2 + 2*k3 + k4)
end

# --- ANIMATION (Inchangé ou presque) ---
println("Génération de l'animation...")
anim = @animate for i in 1:length(x1)
    if i % 100 == 0
        println("Progression: ", i, "/", length(x1))
    end
    Plots.plot(x2[1:i], y2[1:i], aspect_ratio = :equal, size=(600, 600), xlims=(-(l1+l2), (l1+l2)), ylims=(-(l1+l2), (l1+l2)), legend=false, axis=([], false))
               
    Plots.plot!([0, x1[i]], [0, y1[i]], color=:black, linewidth=2)
    Plots.plot!([x1[i], x2[i]], [y1[i], y2[i]], color=:black, linewidth=2)
    
    Plots.scatter!([x1[i]], [y1[i]], color=:red, markersize=5)
    Plots.scatter!([x2[i]], [y2[i]], color=:red, markersize=5)
end every 10

gif(anim, "double_pendulum-rk4.gif", fps = 50)