#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
namespace smc {

class MonoAxisModel {
public:
    MonoAxisModel(double J, double B)
        : rotational_inertia_(J)
        , damping_coefficient_(B) {}

    Eigen::Matrix2d system_matrix() const {
        Eigen::Matrix2d A;
        A << 0, 1, 0, -damping_coefficient_ / rotational_inertia_;
        return A;
    }

    Eigen::Vector2d input_matrix() const { return Eigen::Vector2d{0, 1.0 / rotational_inertia_}; }

    double inertia() const { return rotational_inertia_; }
    double damping() const { return damping_coefficient_; }

private:
    double rotational_inertia_;
    double damping_coefficient_;
};

class TwoAxisGimbalModel {
public:
    TwoAxisGimbalModel(double Jy, double By, double Jp, double Bp)
        : rotational_inertia_y_(Jy)
        , damping_coefficient_y_(By)
        , rotational_inertia_p_(Jp)
        , damping_coefficient_p_(Bp) {}

    Eigen::Matrix4d system_matrix() const {
        double Jy = rotational_inertia_y_;
        double By = damping_coefficient_y_;
        double Jp = rotational_inertia_p_;
        double Bp = damping_coefficient_p_;

        Eigen::Matrix4d A;
        A << 0, 1, 0, 0,       //
            0, -By / Jy, 0, 0, //
            0, 0, 0, 1,        //
            0, 0, 0, -Bp / Jp; //

        return A;
    }

    Eigen::Matrix<double, 4, 2> input_matrix() const {
        double Jy = rotational_inertia_y_;
        double Jp = rotational_inertia_p_;

        Eigen::Matrix<double, 4, 2> B;
        B << 0, 0, 1.0 / Jy, 0, 0, 0, 0, 1.0 / Jp;

        return B;
    }

    double inertia_yaw() const { return rotational_inertia_y_; }
    double inertia_pitch() const { return rotational_inertia_p_; }
    double damping_yaw() const { return damping_coefficient_y_; }
    double damping_pitch() const { return damping_coefficient_p_; }

private:
    double rotational_inertia_y_;
    double damping_coefficient_y_;
    double rotational_inertia_p_;
    double damping_coefficient_p_;
};

} // namespace smc