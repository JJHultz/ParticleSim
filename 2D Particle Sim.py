import numpy as np
import matplotlib.pyplot as plt
import random

# Parameters for the simulation
num_particles = 4 # Number of particles in the simulation
time_step = 0.001  # Time step for the simulation
gravity = -0.1  # Gravity strength (applies downward force)
damping = 0.7 # Cofficient for damping velocity after reversal upon collision
attractive_force = 4
sim_speed = 0.001
colis_force = 1
drag = 0.99

# Simulation boundaries
boundary_x = (-10, 10)
boundary_y = (0, 10)

# Particle class to keep track of position and velocity
class Particle:
    def __init__(self):
        # Initialize particle at random position within boundaries
        self.x = random.uniform(boundary_x[0] + 1, boundary_x[1] - 1)
        self.y = random.uniform(boundary_y[0] + 1, boundary_y[1] - 1)
        
        # Initialize velocity with random components
        self.vx = random.uniform(-1, 1)
        self.vy = random.uniform(-1, 1)
        self.radius = 0.4
        
    def apply_gravity(self):
        self.vy += gravity  # Gravity affects y-velocity

    def apply_drag(self):
        self.vx *= drag
        self.vy *= drag

    def update_position(self):
        # Update particle position based on current velocity
        self.x += self.vx * time_step
        self.y += self.vy * time_step
        
        # Check for boundary collisions and reflect velocity if needed
        if self.x <= boundary_x[0] or self.x >= boundary_x[1]:
            self.vx = -self.vx * damping  # Reverse x-velocity on x-boundary collision
        if self.y <= boundary_y[0] or self.y >= boundary_y[1]:
            self.vy = -self.vy * damping # Reverse y-velocity on y-boundary collision

# Initialize particles
particles = [Particle() for _ in range(num_particles)]

# Set up the plot
plt.ion()  # Interactive mode for live updating
fig, ax = plt.subplots()
ax.set_xlim(-11, 11)
ax.set_ylim(-11, 11)

# Simulation loop
for _ in range(1000):  # Number of steps
    plt.cla()  # Clear the plot for the next frame
    ax.set_xlim(-11,11)
    ax.set_ylim(-11,11)
    
    # Update each particle
    for particle in particles:
        # particle.apply_gravity()
        particle.apply_drag()

        for particle_two in particles: # to apply attractive/repulsive force
            if particle != particle_two: # If the particle we are comparing position to is not the current particle, continue
                dx = particle_two.x - particle.x
                dy = particle_two.y - particle.y
                diam = particle.radius + particle_two.radius

                if ((dx ** 2) + (dy ** 2)) < (diam) ** 2:
                    if abs(dx) < abs(dy):
                        particle.vx = particle.vx * -colis_force
                    else:
                        particle.vy = particle.vy * -colis_force
                else:
                    particle.vx += (particle_two.x - particle.x) * attractive_force
                    particle.vy += (particle_two.y - particle.y) * attractive_force
                # in the future it would be best to create a method to update these values rather than directly modify the attributes of each particle

    for particle in particles:
        particle.update_position()
        
        # Draw the particle
        ax.plot(particle.x, particle.y, 'bo')  # 'bo' for blue circles
    
    plt.pause(sim_speed)  # Pause for real-time effect

plt.ioff()  # Turn off interactive mode
plt.show()  # Keep the final frame open