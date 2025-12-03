using Plots

# Constants
const l1 = 5
const l2 = 4
const m1 = 1
const m2 = 1
const g = 9.81
const delta_t = 0.01
const iterations = 10000

# Base values
theta1 = pi
theta1p = 2.0 # base speed
theta2 = 1
theta2p = 0.0 # base speed

# Position : (0,0) = base
x1 = [sin(theta1) * l1]
y1 = [-l1 * cos(theta1)]
x2 = [x1[end] + sin(theta2) * l2]
y2 = [y1[end] - l2 * cos(theta2)]

for i in 1:iterations
    global theta1, theta1p, theta2, theta2p
    
    # Calcul des lagragiens
    delta = theta1 - theta2

    # --- ÉQUATION 1 ---
    # Forme : M11 * a1 + M12 * a2 = R1
    # On garde les accélérations à gauche, on passe le reste à droite
    
    M11 = (m1 + m2) * l1
    M12 = m2 * l2 * cos(delta)
    # Attention aux signes quand tu passes de l'autre côté du =
    R1  = -m2 * l2 * theta2p^2 * sin(delta) - (m1 + m2) * g * sin(theta1)
    
    # --- ÉQUATION 2 ---
    # Forme : M21 * a1 + M22 * a2 = R2
    
    M21 = l1 * cos(delta)
    M22 = l2
    R2  = l1 * theta1p^2 * sin(delta) - g * sin(theta2)
    
    # --- RÉSOLUTION ---
    # On construit la matrice et le vecteur de droite
    Matrice = [M11 M12; 
            M21 M22]
    Resultats = [R1, R2]
    
    # L'opérateur \ résout le système (trouve a1 et a2) magiquement
    accels = Matrice \ Resultats
    
    # TODO: Euler-explicite -> perte à chaque pas de temps
    theta1pp = accels[1]
    theta2pp = accels[2]
    
    # Mise à jour vitesse
    theta1p += theta1pp * delta_t
    theta2p += theta2pp * delta_t
    
    # Mise à jour position
    theta1 += theta1p * delta_t
    theta2 += theta2p * delta_t
    
    """
    push!(x, cos(theta) * l)
    push!(y, sin(theta) * l)
    """
    push!(x1, l1 * sin(theta1))
    push!(y1, -l1 * cos(theta1))
    push!(x2, x1[end] + l2 * sin(theta2))
    push!(y2, y1[end] - l2 * cos(theta2))
end

anim = @animate for i in 1:length(x1)
    if i % 100 == 0
        println("Progression ", i, "/", length(x1))
    end
	Plots.plot(x2[1:i], y2[1:i], aspect_ratio = :equal, size=(600, 600), axis=([], false))
	Plots.plot!([0, x1[i]], [0, y1[i]], color=:black)
	Plots.plot!([x1[i], x2[i]], [y1[i], y2[i]], color=:black)
    Plots.scatter!(x1[i:i], y1[i:i], xlims=(-9,9), ylims=(-9,9), marker=true, color=:red, legend=false)
    Plots.scatter!(x2[i:i], y2[i:i], xlims=(-9,9), ylims=(-9,9), marker=true, color=:red, legend=false)
end every 10

gif(anim, "doublependulum.gif", fps = 50)