#include <eigen3/Eigen/Dense>

namespace smc {

class GimbalObserver {
public:
    GimbalObserver(double rotational_inertia, double damping_coefficient);

    void update_measurements(double measured_angle, double measured_rate);

    void set_reference_trajectory(
        double target_angle, double target_rate, double target_acceleration);

    const Eigen::Vector2d& current_state_vector() const;

    const Eigen::Vector2d& desired_state_vector() const;

    const Eigen::Vector2d& desired_state_vector_dot() const;

private:
    void estimate_current_state();

    void update_reference_vectors();

    Eigen::Vector2d x_current_;
    Eigen::Vector2d x_desired_;
    Eigen::Vector2d x_desired_dot_;

    double measured_angle_;
    double measured_rate_;

    double target_angle_;
    double target_rate_;
    double target_acceleration_;

    double J_;
    double B_;
};
} // namespace smc