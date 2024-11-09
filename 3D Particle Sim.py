import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import matplotlib.colors as mcolors

# Parameters for the simulation
num_particles = 3  # Number of particles in the simulation
time_step = 0.05  # Time step for the simulation
gravity = 0  # Gravity strength (applies force towards gravity location)
gx = 0
gy = 0
gz = 5
sim_speed = 0.01 # Pause between each simulation frame
damping_force = 0.95 # dampens each collision
drag = 1 # 'air resistance'
collision_adjust = 0.1 # This variable changes how much the position of a particle is adjusted post colision to lessen the chance of two collision frames in a row
attractive_force = 0.3 # Particle attraction force


# Simulation boundaries
boundary_x = (-10, 10)
boundary_y = (-10, 10)
boundary_z = (0, 10)

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
        # self.color = random.choice(list(mcolors.CSS4_COLORS.values()))
        self.color = '#0000FF'
        # print(self.color)
        self.static = False

    def apply_gravity(self):
        self.vx -= (self.x - gx) * gravity
        self.vy -= (self.y - gy) * gravity
        self.vz -= (self.z - gz) * gravity

    def update_position(self):
        # Update particle position based on current velocity
        self.x += self.vx * time_step
        self.y += self.vy * time_step
        self.z += self.vz * time_step
        
        # Check for boundary collisions and reflect velocity if needed
        if self.x <= boundary_x[0]:
            self.vx = -self.vx * damping_force  # Reverse x-velocity on x-boundary collision
            self.x = boundary_x[0]
        if self.x >= boundary_x[1]:
            self.vx = -self.vx * damping_force  # Reverse x-velocity on x-boundary collision
            self.x = boundary_x[1]
        if self.y <= boundary_y[0]:
            self.vy = -self.vy * damping_force  # Reverse x-velocity on x-boundary collision
            self.y = boundary_y[0]
        if self.y >= boundary_y[1]:
            self.vy = -self.vy * damping_force  # Reverse x-velocity on x-boundary collision
            self.y = boundary_y[1]
        if self.z <= boundary_z[0]:
            self.vz = -self.vz * damping_force  # Reverse x-velocity on x-boundary collision
            self.z = boundary_z[0]
        if self.z >= boundary_z[1]:
            self.vz = -self.vz * damping_force  # Reverse x-velocity on x-boundary collision
            self.z = boundary_z[1]

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
        
    def apply_drag(self):
        self.vx = self.vx * drag
        self.vy = self.vy * drag
        self.vz = self.vz * drag

# Initialize particles
particles = [Particle() for _ in range(num_particles)]

# Set up the plot
plt.ion()  # Interactive mode for live updating
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(boundary_x)
ax.set_ylim(boundary_y)
ax.set_zlim(boundary_z)

# Playing around with static particles --------------------------------------------------------------------------------------------------------------------------
if False:
    particles[0].static = True
    particles[0].x = 0
    particles[0].y = 0
    particles[0].z = 5
    particles[0].vx = 1
    particles[0].vy = 1
    particles[0].vz = 1
    particles[0].radius = 2
    particles[0].color = '#FF1111'

# Orbit test ----------------------------------------------------------------------------------------------------------------------------------------------------
if False:
    particles[0].static = True
    particles[0].x = 0
    particles[0].y = 0
    particles[0].z = 5
    particles[0].vx = 2
    particles[0].vy = 2
    particles[0].vz = 2
    particles[0].radius = 8
    particles[0].color = '#FF1111'

    # particles[1].static = True
    particles[1].x = 0
    particles[1].y = 5
    particles[1].z = 5
    particles[1].vx = 0.0001
    particles[1].vy = 0.0001
    particles[1].vz = 3
    particles[1].radius = 0.1
    particles[1].color = '#1111FF'

# Simulation loop
for _ in range(1000):  # Number of steps
    ax.cla()  # Clear the plot for the next frame
    ax.set_xlim(boundary_x)
    ax.set_ylim(boundary_y)
    ax.set_zlim(boundary_z)
    
    # Update each particle's velocity first
    for particle in particles:
        # We'll have to run another for loop to check the positions of every other particle
        for particle_two in particles:
            if not particle.static and particle != particle_two:
                dx = particle_two.x - particle.x
                dy = particle_two.y - particle.y
                dz = particle_two.z - particle.z
                dist = (dx ** 2 + dy ** 2 + dz ** 2) ** 0.5
                radii = particle.radius + particle_two.radius

                if dist < radii:
                    if abs(dx) < abs(dy) and abs(dx) < abs(dz):
                        particle.update_velocity(-particle.vx * damping_force, particle.vy * damping_force, particle.vz * damping_force)
                        particle.x = particle.x + abs(particle.x)/particle.x * collision_adjust * (radii - dist)
                    
                    elif abs(dy) < abs(dx) and abs(dy) < abs(dz):
                        particle.update_velocity(particle.vx * damping_force, -particle.vy * damping_force, particle.vz * damping_force)
                        particle.y = particle.y + abs(particle.y)/particle.y * collision_adjust * (radii - dist)
                    
                    elif abs(dz) < abs(dx) and abs(dz) < abs(dy):
                        particle.update_velocity(particle.vx * damping_force, particle.vy * damping_force, -particle.vz * damping_force)
                        particle.z = particle.z + abs(particle.z)/particle.z * collision_adjust * (radii - dist)

                else:
                    particle.add_velocity(dx/(dist ** 2), dy/(dist ** 2), dz/(dist ** 2))
        

    # Update each particle
    for particle in particles:
        if not particle.static:
            particle.apply_drag()
            particle.apply_gravity()
        particle.update_position()

        # Draw the particle
        ax.scatter(particle.x, particle.y, particle.z, color=particle.color, s=particle.radius*100)  # Blue particles
    
    plt.pause(sim_speed)  # Pause for real-time effect

plt.ioff()  # Turn off interactive mode
plt.show()  # Keep the final frame open
