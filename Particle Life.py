import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import matplotlib.colors as mcolors

# Parameters for the simulation
blue_particles = 20
red_particles = 20
green_particles = 10
purple_particles = 10

red = 0
green = 1
blue = 2
purple = 3

# color attraction matrix
atmat = [[1, -1, -1, -1],
         [-1, 1, -1, -1],
         [-1, -1, 1, -1],
         [-1, -1, -1, 1]]

time_step = 0.001  # Time step for the simulation
gravity = 0  # Gravity strength (applies force towards gravity location)
gx = 0
gy = 0
gz = 5
sim_speed = 0.01 # Pause between each simulation frame
damping_force = 0.95 # dampens each collision
drag = 1 # 'air resistance'
collision_adjust = 0.1 # This variable changes how much the position of a particle is adjusted post colision to lessen the chance of two collision frames in a row
attractive_force = 0.3 # Particle attraction force
search_radius = 3 # radius at which a particle will implement the color attraction matrix

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
        
        # Initialize velocity with random components
        self.vx = 0 # random.uniform(-1, 1)
        self.vy = 0 # random.uniform(-1, 1)
        
        # Initialize the particle property(s)
        self.radius = 0.2
        # self.color = random.choice(list(mcolors.CSS4_COLORS.values()))
        self.color = '#000000'
        self.color_num = 0
        # print(self.color)
        self.static = False

    def apply_gravity(self):
        self.vx -= (self.x - gx) * gravity
        self.vy -= (self.y - gy) * gravity

    def update_position(self):
        # Update particle position based on current velocity
        self.x += self.vx * time_step
        self.y += self.vy * time_step
        
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

    def update_velocity(self, vX, vY):
        # This method will set the particle's new velocity
        self.vx = vX
        self.vy = vY

    def add_velocity(self, vX, vY):
        # This method will add to the particle's velocity
        self.vx += vX
        self.vy += vY
        
    def apply_drag(self):
        self.vx = self.vx * drag
        self.vy = self.vy * drag

# Initialize particles
# Initialize all particles
particles = [Particle() for _ in range(red_particles + blue_particles + green_particles + purple_particles)]

# Assign colors for red
for i in range(0, red_particles):
    particles[i].color = '#FF0000'
    particles[i].color_num = 0

# Then assign green
for i in range(red_particles, red_particles + green_particles):
    particles[i].color = '#00FF00'
    particles[i].color_num = 1

# Then assign blue
for i in range(red_particles + green_particles, red_particles + green_particles + blue_particles):
    particles[i].color = '#0000FF'
    particles[i].color_num = 2

# Finally assign purple
for i in range(red_particles + green_particles + blue_particles,
               red_particles + green_particles + blue_particles + purple_particles):
    particles[i].color = '#FF00FF'
    particles[i].color_num = 3



# Set up the plot
plt.ion()  # Interactive mode for live updating
fig, ax = plt.subplots()
ax.set_xlim(-11, 11)
ax.set_ylim(-11, 11)

# Simulation loop
for _ in range(2000):  # Number of steps
    plt.cla()  # Clear the plot for the next frame
    ax.set_xlim(-11,11)
    ax.set_ylim(-11,11)
    
    # Update each particle
    for particle in particles:
        # particle.apply_gravity()
        # particle.apply_drag()

        for particle_two in particles: # to apply attractive/repulsive force
            if particle != particle_two: # If the particle we are comparing position to is not the current particle, continue
                distance_x = particle_two.x - particle.x
                distance_y = particle_two.y - particle.y
                distance = ((distance_x) ** 2 + (distance_y) ** 2) ** 0.5

                if distance < particle.radius + particle_two.radius:
                    if distance_x < distance_y:
                        particle.update_velocity((particle.vx * damping_force), (-particle.vy * damping_force))
                    else:
                        particle.update_velocity((-particle.vx * damping_force), (particle.vy * damping_force))

                if distance < search_radius:
                    particle.add_velocity((atmat[particle.color_num][particle_two.color_num] * (particle_two.x - particle.x) / (distance ** 2)), (atmat[particle.color_num][particle_two.color_num] * (particle_two.y - particle.y) / (distance ** 2)))

                # in the future it would be best to create a method to update these values rather than directly modify the attributes of each particle

    for particle in particles:
        particle.update_position()
        
        # Draw the particle
        ax.plot(particle.x, particle.y, color=particle.color, marker='o', linestyle='')
    
    plt.pause(sim_speed)  # Pause for real-time effect

plt.ioff()  # Turn off interactive mode
plt.show()  # Keep the final frame open
