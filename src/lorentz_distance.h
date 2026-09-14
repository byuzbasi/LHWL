#ifndef LHWL_LORENTZ_DISTANCE_H
#define LHWL_LORENTZ_DISTANCE_H

#include <algorithm>
#include <cmath>
#include <limits>

// Stable unit-hyperboloid distance. The direct Lorentz contraction is fast,
// but it loses precision for nearby points with large radial coordinates.
inline double lhwl_lorentz_distance(double z0, double z1, double z2,
                                    double w0, double w1, double w2) {
  const long double lz0 = z0;
  const long double lz1 = z1;
  const long double lz2 = z2;
  const long double lw0 = w0;
  const long double lw1 = w1;
  const long double lw2 = w2;

  const long double arg = lz0 * lw0 - lz1 * lw1 - lz2 * lw2;
  const long double scale = std::max(
    static_cast<long double>(1.0),
    std::fabs(lz0 * lw0) + std::fabs(lz1 * lw1) + std::fabs(lz2 * lw2)
  );
  const long double guard =
    64.0L * std::numeric_limits<long double>::epsilon() * scale;

  if (std::isfinite(arg) && arg >= 1.0L && arg - 1.0L > guard) {
    return static_cast<double>(std::acosh(arg));
  }

  const long double sz = std::hypot(lz1, lz2);
  const long double sw = std::hypot(lw1, lw2);
  const long double rz = std::asinh(sz);
  const long double rw = std::asinh(sw);
  long double half_chord_sq = 0.0L;
  if (sz > 0.0L && sw > 0.0L) {
    const long double du1 = lz1 / sz - lw1 / sw;
    const long double du2 = lz2 / sz - lw2 / sw;
    half_chord_sq = 0.25L * (du1 * du1 + du2 * du2);
  }

  const long double radial = std::sinh(0.5L * (rz - rw));
  const long double sinh_half_d_sq =
    radial * radial + sz * sw * half_chord_sq;
  if (!std::isfinite(sinh_half_d_sq) || sinh_half_d_sq < 0.0L) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return static_cast<double>(
    2.0L * std::asinh(std::sqrt(sinh_half_d_sq))
  );
}

#endif
