using DataFrames, CSV, LinearAlgebra, Statistics, Optim

# Constantes
const l1 = 0.09174
const l2 = 0.06933
const g = 9.81
const delta_t = 0.001
const iterations = 4000

# Paramètres vidéo
const frame_finale = 70
const VIDEO_DELTA_T = 0.01
const CALC_DELTA = trunc(Int, VIDEO_DELTA_T / delta_t)
const DIFF = 10 * CALC_DELTA

dfa = DataFrame(CSV.File("mass_a.csv"))
dfb = DataFrame(CSV.File("mass_b.csv"))

#region Fonctions
function RMSE(arrVect1, arrVect2)
    if length(arrVect1) != length(arrVect2)
        throw("error, not same length")
    end
    return sqrt(mean(norm.(arrVect1 .- arrVect2).^2))
end

# Même fonction que dans double_pendulum-rk4 mais on spécifie la masse
function pendulum_dynamics(u, m1, m2)
    theta1 = u[1]
    theta2 = u[2]
    theta1p = u[3]
    theta2p = u[4]
    
    delta = theta1 - theta2

    M11 = (m1 + m2) * l1
    M12 = m2 * l2 * cos(delta)
    R1  = -m2 * l2 * theta2p^2 * sin(delta) - (m1 + m2) * g * sin(theta1)
    
    M21 = l1 * cos(delta)
    M22 = l2
    R2  = l1 * theta1p^2 * sin(delta) - g * sin(theta2)
    
    Matrice = [M11 M12; M21 M22]
    Resultats = [R1, R2]
    
    accels = Matrice \ Resultats 
    return [theta1p, theta2p, accels[1], accels[2]]
end

# Fonction de Coût pour l'optimiseur

function cost_function(params)
    # params contient [m1, m2, theta1_debut, theta2_debut]
    # println("testing avec paramètres: ", params)
    m1, m2, theta1_init, theta2_init = params
    
    u = [theta1_init, theta2_init, 0.0, 0.0]
    
    x1 = Float64[]
    y1 = Float64[]
    x2 = Float64[]
    y2 = Float64[]

    # Simulation RK4
    for i in 1:iterations
        k1 = pendulum_dynamics(u, m1, m2)
        k2 = pendulum_dynamics(u + delta_t*k1/2, m1, m2)
        k3 = pendulum_dynamics(u + delta_t*k2/2, m1, m2)
        k4 = pendulum_dynamics(u + delta_t*k3, m1, m2)
        u += delta_t * (k1 + 2*k2 + 2*k3 + k4) / 6
        
        theta1, theta2 = u[1], u[2]
        push!(x1, l1 * sin(theta1))
        push!(y1, -l1 * cos(theta1))
        push!(x2, x1[end] + l2 * sin(theta2))
        push!(y2, y1[end] - l2 * cos(theta2))
    end
    
    # RMSE masse A
    calc_a = [[x1[i], y1[i]] for i in DIFF:CALC_DELTA:frame_finale*CALC_DELTA + DIFF-CALC_DELTA]
    track_a = [[dfa.x[i], dfa.y[i]] for i in 1:frame_finale]
    rmse_a = RMSE(calc_a, track_a)
    
    # RMSE masse B
    calc_b = [[x2[i], y2[i]] for i in DIFF:CALC_DELTA:frame_finale*CALC_DELTA + DIFF-CALC_DELTA]
    track_b = [[dfb.x[i], dfb.y[i]] for i in 1:frame_finale]
    rmse_b = RMSE(calc_b, track_b)
    
    # retourne valeur qui doit être minimisée = RMSE TOTAL
    return rmse_a + rmse_b
end
#endregion

#region Optimisation

# [m1, m2, theta1_init, theta2_init]
initial_guess = [0.028, 0.0022, pi + 0.0227, pi + 0.0951]

# m1 entre 10g et 50g, m2 entre 0.1g et 6g, angles proches de pi
lower = [0.01, 0.0001, pi, pi]
upper = [0.05, 0.0060, pi + 0.2, pi + 0.2]

println("Optimisation en cours (gradient descent)...")
result = optimize(cost_function, lower, upper, initial_guess, Fminbox(GradientDescent()))

# Résultats
best = Optim.minimizer(result)
println("\nFIN DE L'OPTIMISATION")
println("m1: ", best[1])
println("m2: ", best[2])
println("theta1: ", best[3])
println("theta2: ", best[4])
println("RMSE minimal: ", Optim.minimum(result))

#endregion