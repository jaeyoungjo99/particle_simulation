#include "particle_simulation/particle_manager.hpp"


// Set the bounds within which particles can move
void ParticleManager::SetBounds(float d_min_x, float d_max_x, float d_min_y, float d_max_y) {
    d_min_x_ = d_min_x;
    d_max_x_ = d_max_x;
    d_min_y_ = d_min_y;
    d_max_y_ = d_max_y;

    grid_height_ = static_cast<int>((d_max_y_ - d_min_y_) / grid_cell_size_);
    grid_width_ = static_cast<int>((d_max_x_ - d_min_x_) / grid_cell_size_);

    grid_.resize(grid_height_, std::vector<std::vector<int>>(grid_width_));
}

// Add a new particle to the manager
void ParticleManager::AddParticle(float d_x, float d_y, float d_vx, float d_vy, float d_ax, float d_ay, int id) {
    std::cout << "Add Particle: (" << d_x << ", " << d_y << ")" << std::endl;
    particles_.emplace_back(d_x, d_y, d_vx, d_vy, d_ax, d_ay, id);
}

void ParticleManager::InteractParticle(float x, float y, float vx, float vy, float radius){
    auto [grid_x, grid_y] = GetGridIndex(x,y);

    float min_distance = radius + PARTICLE_RADIUS;

    int i_search_index = static_cast<int>((radius) / grid_cell_size_);

    for (int i_dy = -i_search_index; i_dy <= i_search_index; ++i_dy) {
        for (int i_dx = -i_search_index; i_dx <= i_search_index; ++i_dx) {
            int neighbor_x = grid_x + i_dx;
            int neighbor_y = grid_y + i_dy;

            if (neighbor_x >= 0 && neighbor_x < grid_width_ && neighbor_y >= 0 && neighbor_y < grid_height_) {
                for (int j : grid_[neighbor_y][neighbor_x]) {
                    auto& p2 = particles_[j];

                    float dx = p2.pos.x() - x;
                    float dy = p2.pos.y() - y;
                    float distance = std::sqrt(dx * dx + dy * dy);
                    
                    if (distance < min_distance) {
                        // Handle collision response
                        float nx = dx / distance;
                        float ny = dy / distance;
                        float relative_velocity = (p2.vel.x() - vx) * nx + 
                                                    (p2.vel.y() - vy) * ny;

                        if (relative_velocity < 0) {
                            float m1 = radius*radius; // p1.mass;
                            float m2 = PARTICLE_RADIUS * PARTICLE_RADIUS; // .mass;

                            // Calculate impulse considering different masses
                            float impulse = (1 + RESTITUTION_COEFFICIENT) * relative_velocity / (1/m1 + 1/m2);

                            float impulse_nx = impulse * nx;
                            float impulse_ny = impulse * ny;

                            // Update velocities with impulse
                            p2.vel.noalias() -= Eigen::Vector2f(impulse_nx / m2, impulse_ny / m2);
                        }

                        // Adjust positions to prevent overlap
                        float overlap = min_distance - distance;

                        float correction_x = overlap * nx;
                        float correction_y = overlap * ny;

                        // Use noalias() to optimize the vector addition
                        float velocity_factor = 1.0f; // Adjust this factor to control the influence on velocity
                        p2.pos.noalias() += Eigen::Vector2f(correction_x, correction_y);
                        p2.vel.noalias() += Eigen::Vector2f(correction_x, correction_y) * velocity_factor;

                    }
                }
            }
        }
    }
}

std::pair<int, int> ParticleManager::GetGridIndex(const float &x, const float &y) const {
    int grid_x = static_cast<int>((x - d_min_x_) / grid_cell_size_);
    int grid_y = static_cast<int>((y - d_min_y_) / grid_cell_size_);
    return {grid_x, grid_y};
}

void ParticleManager::AssignParticlesToGrid() {
    // Clear the grid
    for (auto& row : grid_) {
        for (auto& cell : row) {
            cell.clear();
        }
    }

    // Assign particles to grid
    for (int i = 0; i < particles_.size(); ++i) {
        // Reset checked flag
        particles_[i].resetMutualId();

        const auto& p = particles_[i];
        auto [grid_x, grid_y] = GetGridIndex(p.pos.x(), p.pos.y());
        if (grid_x >= 0 && grid_x < grid_width_ && grid_y >= 0 && grid_y < grid_height_) {
            grid_[grid_y][grid_x].push_back(i);
        }
    }
}

void ParticleManager::UpdateParticleIteration(float d_time_step, int iteration){
    for (int i = 0; i < iteration; ++i) {
        if(i == 0){
            UpdateParticles(d_time_step);
        }
        else{
            UpdateParticles(0.0);
        }
    }
}

// Update all particles
void ParticleManager::UpdateParticles(float d_time_step) {

    AssignParticlesToGrid();

    float min_distance = PARTICLE_RADIUS * 2.0;

    for (auto& p1 : particles_) {
        p1.Update(d_time_step);

        // Check for boundary collision and reverse velocity if needed
        if (p1.pos.x() < d_min_x_ + PARTICLE_RADIUS) {
            p1.pos.x() = d_min_x_ + PARTICLE_RADIUS;
            p1.vel.x() = -p1.vel.x() * DAMPING_FACTOR;
        }
        if (p1.pos.x() > d_max_x_ - PARTICLE_RADIUS) {
            p1.pos.x() = d_max_x_ - PARTICLE_RADIUS;
            p1.vel.x() = -p1.vel.x() * DAMPING_FACTOR;
        }
        if (p1.pos.y() < d_min_y_ + PARTICLE_RADIUS) {
            p1.pos.y() = d_min_y_ + PARTICLE_RADIUS;
            p1.vel.y() = -p1.vel.y() * DAMPING_FACTOR;
        }
        if (p1.pos.y() > d_max_y_ - PARTICLE_RADIUS) {
            p1.pos.y() = d_max_y_ - PARTICLE_RADIUS;
            p1.vel.y() = -p1.vel.y() * DAMPING_FACTOR;
        }


        auto [grid_x, grid_y] = GetGridIndex(p1.pos.x(), p1.pos.y());

        // Check this cell and adjacent cells
        for (int i_dy = -1; i_dy <= 1; ++i_dy) {
            for (int i_dx = -1; i_dx <= 1; ++i_dx) {
                int neighbor_x = grid_x + i_dx;
                int neighbor_y = grid_y + i_dy;

                if (neighbor_x >= 0 && neighbor_x < grid_width_ && neighbor_y >= 0 && neighbor_y < grid_height_) {
                    for (int j : grid_[neighbor_y][neighbor_x]) {
                        auto& p2 = particles_[j];

                        // Skip self
                        if (p1.id() == p2.id()) continue;

                        // Skip mutually checked particles
                        if (p1.mutualId() == p2.id() || p2.mutualId() == p1.id()) continue;

                        float dx = p2.pos.x() - p1.pos.x();
                        float dy = p2.pos.y() - p1.pos.y();
                        float distance = std::sqrt(dx * dx + dy * dy);
                        

                        if (distance < min_distance) {
                            p1.setMutualId(p2.id());
                            p2.setMutualId(p1.id());

                            // Handle collision response
                            float nx = dx / distance;
                            float ny = dy / distance;
                            float relative_velocity = (p2.vel.x() - p1.vel.x()) * nx + 
                                                      (p2.vel.y() - p1.vel.y()) * ny;

                            if (relative_velocity < 0) {
                                float impulse = (1 + RESTITUTION_COEFFICIENT) * relative_velocity / (1 + 1);  // Assuming equal mass

                                float impulse_nx = impulse * nx;
                                float impulse_ny = impulse * ny;

                                // Use noalias() to optimize the vector addition
                                p1.vel.noalias() += Eigen::Vector2f(impulse_nx, impulse_ny);
                                p2.vel.noalias() -= Eigen::Vector2f(impulse_nx, impulse_ny);
                            }

                            // Adjust positions to prevent overlap
                            float overlap = min_distance - distance;
                            float correction_factor = 0.5; // Adjust this factor as needed

                            float correction_x = overlap * nx * correction_factor;
                            float correction_y = overlap * ny * correction_factor;

                            // Use noalias() to optimize the vector addition
                            p1.pos.noalias() -= Eigen::Vector2f(correction_x, correction_y);
                            p2.pos.noalias() += Eigen::Vector2f(correction_x, correction_y);

                        }

                        // Calculate velocity difference for viscosity
                        Eigen::Vector2f velocity_difference = p2.vel - p1.vel;

                        // Apply viscosity force
                        Eigen::Vector2f viscosity_force = VISCOSITY_COEFFICIENT * velocity_difference;

                        // Update velocities with viscosity using noalias()
                        p1.vel.noalias() += viscosity_force * d_time_step;
                        p2.vel.noalias() -= viscosity_force * d_time_step;
                    }
                }
            }
        }
    }
}

std::vector<Particle> ParticleManager::GetParticles() const {
    return particles_;
}

int ParticleManager::GetParticleCount() const {
    return particles_.size();
}