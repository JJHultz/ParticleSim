import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import matplotlib.colors as mcolors

# Parameters for the simulation
blue_particles =   100
red_particles =    100
green_particles =  100
purple_particles = 100

RED = 0
GREEN = 1
BLUE = 2
PURPLE = 3

# (CHANGED) dictionary for pairwise interactions
# Key: (color1, color2)
# Value: a float representing how strongly color1 is attracted (+) or repelled (-) by color2
color_interactions = {
    (RED,   RED):        1.0,
    (RED,   GREEN):     -1.0,
    (RED,   BLUE):      -1.0,
    (RED,   PURPLE):    -1.0,
    (GREEN,   RED):      1.0, # Green is attracted to Red 
    (GREEN,   GREEN):   -1.0,
    (GREEN,   BLUE):    -1.0,
    (GREEN,   PURPLE):   1.0,
    (BLUE,   RED):       1.0,
    (BLUE,   GREEN):    -1.0,
    (BLUE,   BLUE):     -1.0,
    (BLUE,   PURPLE):    1.0,
    (PURPLE,   RED):     0.0,
    (PURPLE,   GREEN):  -1.0,
    (PURPLE,   BLUE):    0.0,
    (PURPLE,   PURPLE): -1.0
}

# color_interactions = {
#     (RED,   RED):        1.0,
#     (RED,   GREEN):      1.0,
#     (RED,   BLUE):       1.0,   # red -> blue is attracted
#     (RED,   PURPLE):     1.0,
#     (GREEN,   RED):      1.0,
#     (GREEN,   GREEN):    -1.0,
#     (GREEN,   BLUE):     1.0,
#     (GREEN,   PURPLE):   1.0,
#     (BLUE,   RED):       1.0,
#     (BLUE,   GREEN):     0.0,
#     (BLUE,   BLUE):      1.0,
#     (BLUE,   PURPLE):    1.0,
#     (PURPLE,   RED):     0.0,
#     (PURPLE,   GREEN):   1.0,
#     (PURPLE,   BLUE):    0.0,
#     (PURPLE,   PURPLE):  1.0
# }

def get_interaction_force(c1, c2):
    return color_interactions.get((c1, c2), 0)

time_step = 0.05  # Time step for the simulation
gravity = 0  # Gravity strength (applies force towards gravity location)
gx = 0
gy = 0
gz = 5
sim_speed = 0.00001 # Pause between each simulation frame
damping_force = 0.95 # dampens each collision
drag = 0.95 # 'air resistance'
collision_adjust = 0.1 # This variable changes how much the position of a particle is adjusted post colision to lessen the chance of two collision frames in a row
attractive_force = 0.3 # Particle attraction force
search_radius = 1 # radius at which a particle will implement the color attraction matrix

# Simulation boundaries
boundary_x = (-20, 20)
boundary_y = (-20, 20)

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
fig.set_facecolor('lightgray')

epsilon = 0.05    # Softening term for inverse-square
restitution = 0.5 # Collision coefficient of restitution
max_speed = 5.0   # Cap on velocity

# Simulation loop
for _ in range(2000):  # Number of steps
    plt.cla()  # Clear the plot for the next frame
    ax.set_xlim(-21,21)
    ax.set_ylim(-21,21)
    
    # We can do a double loop so we only handle each pair once
    for i, p1 in enumerate(particles):
        # p1.apply_gravity()
        # p1.apply_drag()

        for j in range(i+1, len(particles)):
            p2 = particles[j]

            dx = p2.x - p1.x
            dy = p2.y - p1.y
            dist = (dx*dx + dy*dy)**0.5

            # --- Collision detection/response ---
            if dist < p1.radius + p2.radius:
                # (CHANGED) Use a normal impulse approach
                if dist == 0:
                    # If two particles happen to have the exact same position,
                    # skip or nudge them apart to avoid division by zero
                    dist = p1.radius + p2.radius
                    dx = p1.radius
                    dy = p1.radius

                nx = dx / dist
                ny = dy / dist
                # Relative velocity along the normal
                rel_vx = p1.vx - p2.vx
                rel_vy = p1.vy - p2.vy
                rel_vel_normal = rel_vx * nx + rel_vy * ny

                # If they're moving toward each other, apply an impulse
                if rel_vel_normal < 0:
                    # j is the magnitude of impulse
                    j = -(1.0 + restitution) * rel_vel_normal
                    # If you had masses, you'd divide by sum(1/m_i)
                    # For simplicity, assume equal mass => half to each
                    jx = j * nx * 0.5
                    jy = j * ny * 0.5

                    # Apply impulse
                    p1.vx += jx
                    p1.vy += jy
                    p2.vx -= jx
                    p2.vy -= jy

                # (CHANGED) Position Correction
                overlap = (p1.radius + p2.radius) - dist
                if overlap > 0:
                    correction = overlap / 2.0
                    p1.x -= correction * nx
                    p1.y -= correction * ny
                    p2.x += correction * nx
                    p2.y += correction * ny


            # (CHANGED) color-based attraction/repulsion
            if dist < search_radius and dist > 0:
                force_coef = get_interaction_force(p1.color_num, p2.color_num)
                force_mag = force_coef / (dist ** 2 + epsilon)

                # apply force to p1
                fx = force_mag * dx
                fy = force_mag * dy
                p1.add_velocity(fx, fy)

                # If you want equal & opposite reaction:
                # to keep momentum balanced
                force_coef = get_interaction_force(p2.color_num, p1.color_num)
                force_mag = force_coef / (dist ** 2 + 2)
                fx2 = force_mag * dx
                fy2 = force_mag * dy
                p2.add_velocity(fx2, fy2)

    # Now update positions (and optionally clamp velocity)
    for p in particles:
        # (CHANGED) Optionally clamp velocity
        speed = (p.vx*p.vx + p.vy*p.vy)**0.5
        if speed > max_speed:
            scale = max_speed / speed
            p.vx *= scale
            p.vy *= scale

        p.apply_drag()
        p.update_position()

        # Plot
        ax.plot(p.x, p.y, color=p.color, marker='o', linestyle='')

    plt.pause(sim_speed)

plt.ioff()
plt.show()