using DataFrames, CSV, LinearAlgebra, Statistics, Optim, LineSearches

include("double_pendulum.jl")

# Constantes
const l1 = 0.09174
const l2 = 0.06933
const g = 9.81
const delta_t = 0.001
const iterations = 2000

# Paramètres vidéo
dfa_raw = DataFrame(CSV.File("mass_a_200.csv"))
dfb_raw = DataFrame(CSV.File("mass_b_200.csv"))

const track_a_x = Float64.(dfa_raw.x)
const track_a_y = Float64.(dfa_raw.y)
const track_b_x = Float64.(dfb_raw.x)
const track_b_y = Float64.(dfb_raw.y)

const frame_finale = length(track_a_x)
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

# Fonction de Coût pour l'optimiseur

function cost_function(params)
    # params contient [m1, m2, theta1_debut, theta2_debut, theta1p_debut, theta2p_debut]
    m1, m2, theta1p_init, theta2p_init = params

    if m1 <= 0 || m2 <= 0 || theta1p_init <= 0 || theta2p_init <= 0
        return Inf
    end
    
    # angles de base fixes
    u = (3.1728, 3.225, theta1p_init, theta2p_init)
    
    sum_sq_err_a = 0.0
    sum_sq_err_b = 0.0
    track_idx = 1

    # Simulation RK4
    for i in 1:iterations

        if (i - 1) % CALC_DELTA == 0
            t1, t2 = u[1], u[2]

            x1 = l1 * sin(t1)
            y1 = -l1 * cos(t1)
            x2 = x1 + l2 * sin(t2)
            y2 = y1 - l2 * cos(t2)

            sum_sq_err_a += (x1 - track_a_x[track_idx])^2 + (y1 - track_a_y[track_idx])^2
            sum_sq_err_b += (x2 - track_b_x[track_idx])^2 + (y2 - track_b_y[track_idx])^2

            track_idx += 1
        end

        k1 = pendulum_dynamics(u, m1, m2)
        k2 = pendulum_dynamics(u .+ (delta_t/2) .* k1, m1, m2)
        k3 = pendulum_dynamics(u .+ (delta_t/2) .* k2, m1, m2)
        k4 = pendulum_dynamics(u .+ delta_t .* k3, m1, m2)
        u = u .+ (delta_t/6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)

        # Stability check to avoid DomainError
        if !all(isfinite, u) || any(abs.(x) > 1e5 for x in u)
            return 1e9
        end
    end
        
    count = track_idx - 1
    return sqrt(sum_sq_err_a / count) + sqrt(sum_sq_err_b / count)
end
#endregion

#region Optimisation

# [m1, m2, theta1_init, theta2_init, theta1p, theta2p]
initial_guess = [0.0306, 0.003592, 0, 0.2] # valeur trouvée avant avec petit rmse
lower = [0.005, 0.001, 0, 0]
upper = [0.06, 0.01, 0.5, 1]

println("Optimisation...")

result = optimize(cost_function, lower, upper, ParticleSwarm(n_particles=10000), 
                  Optim.Options(show_trace=true, iterations=100, show_every=2))

best = Optim.minimizer(result)

println("\n--- RÉSULTATS ---")
println("Masse 1: ", round(best[1], digits=5))
println("Masse 2: ", round(best[2], digits=5))
println("Theta 1p: ", round(best[3], digits=5))
println("Theta 2p: ", round(best[4], digits=5))
println("RMSE final: ", round(Optim.minimum(result), digits=5))

#endregion

best = round.(best, digits=5)
load_double_pendulum_mp4(best[1], best[2], 3.1728, 3.225, best[3], best[4])