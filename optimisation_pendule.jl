using DataFrames, CSV, LinearAlgebra, Statistics, Optim, LineSearches

# Constantes
const l1 = 0.09174
const l2 = 0.06933
const g = 9.81
const delta_t = 0.001
const iterations = 2000

# Paramètres vidéo
dfa = DataFrame(CSV.File("mass_a_200.csv"))
dfb = DataFrame(CSV.File("mass_b_200.csv"))

const frame_finale = length(dfa.x)
const VIDEO_DELTA_T = 0.01
const CALC_DELTA = trunc(Int, VIDEO_DELTA_T / delta_t)
const DIFF = 10 * CALC_DELTA

#region Fonctions
function RMSE(arrVect1, arrVect2)
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

    # Cramer -> https://fr.wikipedia.org/wiki/R%C3%A8gle_de_Cramer  (éviter de créer trop de tableaux)
    det = M11 * M22 - M12 * M21
    acc1 = (M22 * R1 - M12 * R2) / det
    acc2 = (-M21 * R1 + M11 * R2) / det
    
    return [theta1p, theta2p, acc1, acc2]
end

# Fonction de Coût pour l'optimiseur

function cost_function(params)
    # params contient [m1, m2, theta1_debut, theta2_debut, theta1p_debut, theta2p_debut]
    m1, m2, theta1_init, theta2_init, theta1p_init, theta2p_init = params
    
    u = [theta1_init, theta2_init, theta1p_init, theta2p_init]
    
    sum_sq_err_a = 0.0
    sum_sq_err_b = 0.0
    track_idx = 1

    # Simulation RK4
    for i in 1:iterations

        if (i - 1) % CALC_DELTA == 0
            t1, t2 = u[1], u[2]

            x1_val = l1 * sin(t1)
            y1_val = -l1 * cos(t1)
            x2_val = x1_val + l2 * sin(t2)
            y2_val = y1_val - l2 * cos(t2)

            sum_sq_err_a += (x1_val - dfa.x[track_idx])^2 + (y1_val - dfa.y[track_idx])^2
            sum_sq_err_b += (x2_val - dfb.x[track_idx])^2 + (y2_val - dfb.y[track_idx])^2

            track_idx += 1
        end

        k1 = pendulum_dynamics(u, m1, m2)
        k2 = pendulum_dynamics(u + delta_t*k1/2, m1, m2)
        k3 = pendulum_dynamics(u + delta_t*k2/2, m1, m2)
        k4 = pendulum_dynamics(u + delta_t*k3, m1, m2)
        u += delta_t * (k1 + 2*k2 + 2*k3 + k4) / 6
    end
        
    count = track_idx - 1
    return sqrt(sum_sq_err_a / count) + sqrt(sum_sq_err_b / count)
end
#endregion

#region Optimisation

# [m1, m2, theta1_init, theta2_init, theta1p, theta2p]
initial_guess = [0.0306, 0.003592, 3.08977, 3.35064, 1.09264, 0.082822] # valeur trouvée avant avec petit rmse
lower = [0.005, 0.0005, pi - 0.5, pi - 0.5, -1.5, -1.0]
upper = [0.04, 0.01, pi + 0.5, pi + 0.5, 1.5, 1.0]

println("Optimisation (LBFGS)...")

result = optimize(cost_function, lower, upper, initial_guess, Fminbox(LBFGS()), 
                  Optim.Options(show_trace=true, iterations=1000, show_every=5))

best = Optim.minimizer(result)

println("\n--- RÉSULTATS ---")
println("Masse 1: ", round(best[1], digits=5))
println("Masse 2: ", round(best[2], digits=5))
println("Theta 1: ", round(best[3], digits=5))
println("Theta 2: ", round(best[4], digits=5))
println("Theta 1p: ", round(best[5], digits=5))
println("Theta 2p: ", round(best[6], digits=5))
println("RMSE final: ", round(Optim.minimum(result), digits=5))

#endregion