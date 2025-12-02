# Using lagrangian, found that theta.. = -g/l * sin(theta)

using Plots

# Constants
const l = 5
const g = 9.81
const delta_t = 0.01
const iterations = 800

# Base values
theta = 3.1415
thetap = 0.6 # base speed

# Position : (0,0) = base
x = [sin(theta) * l]
y = [-l * cos(theta)]

for i in 1:iterations
    global theta, thetap
    
    # Calcul accel
    thetapp = -(g/l) * sin(theta)
    
    # Mise à jour vitesse
    thetap += thetapp * delta_t
    
    # Mise à jour position
    theta += thetap * delta_t
    
    """
    push!(x, cos(theta) * l)
    push!(y, sin(theta) * l)
    """
    push!(x, l * sin(theta))
    push!(y, -l * cos(theta))
end

anim = @animate for i in 1:length(x)
    if i % 100 == 0
        println("Progression ", i, "/", length(x))
    end
	Plots.plot(x[1:i], y[1:i], aspect_ratio = :equal, size=(600, 600), axis=([], false))
	Plots.plot!([0, x[i]], [0, y[i]], color=:black)
    Plots.scatter!(x[i:i], y[i:i], xlims=(-6,6), ylims=(-6,6), marker=true, color=:red, legend=false)
end

gif(anim, "basicpendulum.gif", fps = 100)