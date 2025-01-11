#include "particle_simulation/particle.hpp"


// Constructor
Particle::Particle(float d_x, float d_y, float d_vx, float d_vy, float d_ax, float d_ay, int id)
    : pos(d_x, d_y), vel(d_vx, d_vy), acc(d_ax, d_ay), id_(id) {}

// Update the particle's position based on its velocity and acceleration
void Particle::Update(float d_time_step) {
    // Update velocity based on acceleration
    vel.x() += acc.x() * d_time_step;
    vel.y() += acc.y() * d_time_step;

    // Update position based on velocity
    pos.x() += vel.x() * d_time_step;
    pos.y() += vel.y() * d_time_step;
}
