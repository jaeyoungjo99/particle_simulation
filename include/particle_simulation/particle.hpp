#ifndef PARTICLE_SIMULATION_PARTICLE_HPP
#define PARTICLE_SIMULATION_PARTICLE_HPP

#include <iostream>
#include <cmath>
#include <Eigen/Dense>

#include "particle_simulation/config.hpp"

class Particle {
public:
    // Constructor
    Particle(float d_x, float d_y, float d_vx, float d_vy, float d_ax, float d_ay, int id);

    // Update the particle's position based on its velocity and acceleration
    void Update(float d_time_step);

    int id() const { return id_; }
    int mutualId() const { return i_mutual_id_; }
    void setMutualId(int id) { i_mutual_id_ = id; }
    void resetMutualId() { i_mutual_id_ = -1; }

    bool isChecked() const { return checked_; }
    void setChecked(bool checked) { checked_ = checked; }

    Eigen::Vector2f pos;
    Eigen::Vector2f vel;
    Eigen::Vector2f acc;

private:

    int id_;
    int i_mutual_id_ = -1;

    bool checked_ = false;
    
};
#endif // PARTICLE_SIMULATION_PARTICLE_HPP

