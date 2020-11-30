#ifndef MATH_HPP
#define MATH_HPP

#include <Eigen/Core>

namespace math {

constexpr double PRECISION = 1.0e-14;

template <typename DerivedLHS, typename DerivedRHS>
bool isEqual(const Eigen::DenseBase<DerivedLHS> &lhs,
             const Eigen::DenseBase<DerivedRHS> &rhs) {
  return lhs.isapprox(lhs, rhs, PRECISION);
}

template <typename Derived>
bool isZero(const Eigen::DenseBase<Derived> &val) {
  return val.isZero(PRECISION);
}

} // namespace math
#endif /* end of include guard */
