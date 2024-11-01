import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random

# Parameters for the simulation
num_particles = 10  # Number of particles in the simulation
time_step = 0.05  # Time step for the simulation
gravity = -0.1  # Gravity strength (applies downward force)
sim_speed = 0.01 # Pause between each simulation frame
damping_force = 0.7 # dampens each collision


# Simulation boundaries
boundary_x = (-10, 10)
boundary_y = (0, 10)
boundary_z = (-10, 10)

# Particle class to keep track of position and velocity
class Particle:
    def __init__(self):
        # Initialize particle at random position within boundaries
        self.x = random.uniform(boundary_x[0] + 1, boundary_x[1] - 1)
        self.y = random.uniform(boundary_y[0] + 1, boundary_y[1] - 1)
        self.z = random.uniform(boundary_z[0] + 1, boundary_z[1] - 1)
        
        # Initialize velocity with random components
        self.vx = random.uniform(-1, 1)
        self.vy = random.uniform(-1, 1)
        self.vz = random.uniform(-1, 1)
        
        # Initialize the particle property(s)
        self.radius = 0.2

    def apply_gravity(self):
        self.vz += gravity  # Gravity affects z-velocity only

    def update_position(self):
        # Update particle position based on current velocity
        self.x += self.vx * time_step
        self.y += self.vy * time_step
        self.z += self.vz * time_step
        
        # Check for boundary collisions and reflect velocity if needed
        if self.x <= boundary_x[0] or self.x >= boundary_x[1]:
            self.vx = -self.vx * damping_force  # Reverse x-velocity on x-boundary collision
        if self.y <= boundary_y[0] or self.y >= boundary_y[1]:
            self.vy = -self.vy * damping_force  # Reverse y-velocity on y-boundary collision
        if self.z <= boundary_z[0] or self.z >= boundary_z[1]:
            self.vz = -self.vz * damping_force  # Reverse z-velocity on z-boundary collision

    def update_velocity(self, vX, vY, vZ):
        # This method will set the particle's new velocity
        self.vx = vX
        self.vy = vY
        self.vz = vZ

    def add_velocity(self, vX, vY, vZ):
        # This method will add to the particle's velocity
        self.vx += vX
        self.vy += vY
        self.vz += vZ

# Initialize particles
particles = [Particle() for _ in range(num_particles)]

# Set up the plot
plt.ion()  # Interactive mode for live updating
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(boundary_x)
ax.set_ylim(boundary_y)
ax.set_zlim(boundary_z)

# Simulation loop
for _ in range(1000):  # Number of steps
    ax.cla()  # Clear the plot for the next frame
    ax.set_xlim(boundary_x)
    ax.set_ylim(boundary_y)
    ax.set_zlim(boundary_z)
    
    # Update each particle's velocity first
    for particle in particles:
        # particle.apply_gravity()
        particle.update_position()
        

    # Update each particle
    for particle in particles:
        # particle.apply_gravity()
        particle.update_position()

        # Draw the particle
        ax.scatter(particle.x, particle.y, particle.z, color='b')  # Blue particles
    
    plt.pause(sim_speed)  # Pause for real-time effect

plt.ioff()  # Turn off interactive mode
plt.show()  # Keep the final frame open
