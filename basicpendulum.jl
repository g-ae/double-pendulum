# Using lagrangian, found that theta.. = -g/l * sin(theta)

using Plots

# Constants
const l = 5
const g = 9.81
const delta_t = 0.001
const iterations = 5000

# Base values
theta = 1.0
thetap = 0.0 # base speed

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

anim = @animate for i in 1:20:length(x)
	Plots.plot(x[1:i], y[1:i])
    Plots.scatter!(x[i:i], y[i:i], xlims=(-6,6), ylims=(-6,6), marker=true, color=:red, legend=false)
end

gif(anim, "basic_pendulum.gif", fps = 30)