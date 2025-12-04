#include <cstdlib>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>

namespace smc {
class SMCCalculator {
public:
    SMCCalculator(double c, double k, double phi)
        : c_(c)
        , k_(k)
        , phi_(phi) {
        C_transpose_ << c_, 1.0;
    }

    double update_control_torque(
        const Eigen::Vector2d& x, const Eigen::Vector2d& x_d, const Eigen::Vector2d& dot_x_d,
        const Eigen::Matrix2d& A, const Eigen::Vector2d& B) {

        auto s = calc_sliding_surface(x, x_d);
        auto u_eq = calc_equal_control(x, dot_x_d, A, B);
        auto u_sw = calc_switch_control(s, B);
        return u_eq + u_sw;
    }

private:
    double calc_sliding_surface(const Eigen::Vector2d& x, const Eigen::Vector2d& x_d) const {
        // 计算滑模面 s = C^T * e
        Eigen::Vector2d e = x - x_d;
        return C_transpose_ * e;
    }

    double calc_equal_control(
        const Eigen::Vector2d& x, const Eigen::Vector2d& dot_x_d, const Eigen::Matrix2d& A,
        const Eigen::Vector2d& B) const {
        // 计算等效控制 u_eq = (C^T * B)^-1 * C^T * (dot_x_d - A * x)
        double term_eq = (C_transpose_ * (dot_x_d - A * x));
        double s = 1.0 / (C_transpose_ * B) * term_eq;
        return s;
    }

    double calc_switch_control(double sliding_surface, const Eigen::Vector2d& B) const {
        // 计算切换控制 u_sw = -(C^T * B)^-1 * K * sat(s/Phi)
        double sat_val = sat(sliding_surface, phi_);
        double u_eq = -1.0 / (C_transpose_ * B) * k_ * sat_val;
        return u_eq;
    }

    static double sgn(double s) {
        if (s == 0) {
            return 0;
        } else if (s > 0) {
            return 1;
        } else {
            return -1;
        }
    }

    static double sat(double s, double phi) {
        if (std::abs(s) < phi) {
            return s / phi;
        } else {
            return sgn(s);
        }
    }

    double c_;
    double k_;
    double phi_;
    Eigen::RowVector2d C_transpose_;
};
} // namespace smc