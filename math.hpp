#ifndef MATH_HPP
#define MATH_HPP

#include <Eigen/Core>
#include <cassert>

namespace math {

constexpr double PRECISION = 1.0e-14;

template <typename DerivedLHS, typename DerivedRHS>
bool isEqual(const Eigen::DenseBase<DerivedLHS> &lhs,
             const Eigen::DenseBase<DerivedRHS> &rhs) {
  return lhs.isApprox(rhs, PRECISION);
}

inline bool isZero(double val) {
  return std::abs(val) < PRECISION;
}

inline double triangleArea(
    const Eigen::VectorXd &a,
    const Eigen::VectorXd &b,
    const Eigen::VectorXd &c)
{
  assert(a.size() == b.size());
  assert(b.size() == c.size());
  if (a.size() == 2) {
    Eigen::Vector2d A = b;
    A -= a;
    Eigen::Vector2d B = c;
    B -= a;
    return 0.5 * (A(0) * B(1) - A(1) * B(0));
  } else {
    assert(a.size() == 3);
    Eigen::Vector3d A = b;
    A -= a;
    Eigen::Vector3d B = c;
    B -= a;
    return 0.5 * A.cross(B).norm();
  }
}

} // namespace math
#endif /* end of include guard */
