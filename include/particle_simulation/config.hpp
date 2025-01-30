#ifndef PARTICLE_SIMULATION_CONFIG_HPP
#define PARTICLE_SIMULATION_CONFIG_HPP

#define WINDOW_WIDTH 900.0
#define WINDOW_HEIGHT 900.0
#define CIRCLE_SEGMENTS 10

#define MIN_GEN_DT 0.01
#define FRAME_DT 0.01
#define OPT_ITER 3

#define PARTICLE_RADIUS 2.0        // Radius of each particle
#define DAMPING_FACTOR 0.9      // Factor to reduce velocity after collision with walls
#define RESTITUTION_COEFFICIENT 1.0  // Bounciness factor for particle collisions 1.0: complete colision
#define VISCOSITY_COEFFICIENT 10.0  // 0~0.1 pure water, ~99 maximum
#define GRAVITY 30              // Downward acceleration force

#define INTERACT_RADIUS 20


#endif // PARTICLE_SIMULATION_CONFIG_HPP
