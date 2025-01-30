#ifndef PARTICLE_SIMULATION_PARTICLE_MANAGER_HPP
#define PARTICLE_SIMULATION_PARTICLE_MANAGER_HPP

#include "particle_simulation/particle.hpp"
#include "particle_simulation/config.hpp"
#include <vector>
#include <utility>
#include <chrono>



class ParticleManager {
public:
    void SetBounds(float d_min_x, float d_max_x, float d_min_y, float d_max_y);
    void AddParticle(float d_x, float d_y, float d_vx, float d_vy, float d_ax, float d_ay, int id);
    void InteractParticle(float x, float y, float vx, float vy, float radius);
    void UpdateParticleIteration(float d_time_step, int iteration = 1);
    void UpdateParticles(float d_time_step);

    // Get all particle positions
    std::vector<Particle> GetParticles() const;

    int GetParticleCount() const;

private:
    std::pair<int, int> GetGridIndex(const float &x, const float &y) const;
    void AssignParticlesToGrid();

    std::vector<Particle> particles_; // List of particles
    float grid_cell_size_ = PARTICLE_RADIUS * 1.7; // Example cell size
    int grid_width_, grid_height_;
    std::vector<std::vector<std::vector<int>>> grid_; // 3D vector: [grid_y][grid_x][particle_indices]

    float d_min_x_, d_max_x_; // X bounds
    float d_min_y_, d_max_y_; // Y bounds
};
#endif // PARTICLE_SIMULATION_PARTICLE_MANAGER_HPP

