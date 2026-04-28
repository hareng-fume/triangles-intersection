#include "triangle_intersect.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>

namespace _Details
{
  constexpr double EPS = 1e-8;

  inline int Sign(double i_val, double i_eps = EPS)
  {
    return (i_val > i_eps) - (i_val < -i_eps);
  }

  inline bool Equal(double i_val1, double i_val2, double i_eps = EPS)
  {
    return std::abs(i_val1 - i_val2) < i_eps;
  }

  template <std::size_t NDim>
  struct Point {std::array<double, NDim> coords;};

  using Point2D = Point<2>;
  using Point3D = Point<3>;

  template <std::size_t NDim>
  Point<NDim> operator-(const Point<NDim>& i_p1, const Point<NDim>& i_p2)
  {
    Point<NDim> res{};
    for(std::size_t i = 0; i < NDim; ++i)
      res.coords[i] = i_p1.coords[i] - i_p2.coords[i];
    return res;
  }

  template <std::size_t NDim>
  Point<NDim> operator/(const Point<NDim>& i_p, double i_denom)
  {
    Point<NDim> res{};
    if(std::abs(i_denom) > EPS)
    {
      auto inv = 1.0 / i_denom;
      for(std::size_t i = 0; i < NDim; ++i)
        res.coords[i] = i_p.coords[i] * inv;
    }
    return res;
  }

  template <std::size_t NDim>
  inline bool Equal(const Point<NDim>& i_p1, const Point<NDim>& i_p2)
  {
    for(std::size_t i = 0; i < NDim; ++i)
      if(!Equal(i_p1.coords[i], i_p2.coords[i]))
        return false;
    return true;
  }

  template <std::size_t NDim>
  using Vector = Point<NDim>;
  using Vector2D = Vector<2>;
  using Vector3D = Vector<3>;

  template <std::size_t NDim>
  inline double LengthSquared(const Vector<NDim>& i_vec)
  {
    auto sq_sum{0.0};
    for(std::size_t i = 0; i < NDim; ++i)
      sq_sum += i_vec.coords[i]*i_vec.coords[i];
	return sq_sum;
  }

  template <std::size_t NDim>
  inline double Length(const Vector<NDim>& i_vec)
  {
    auto len_squared = LengthSquared(i_vec);
    return std::sqrt(len_squared);
  }

  template <std::size_t NDim>
  inline Vector<NDim> Normalize(const Vector<NDim>& i_vec)
  {
    if (auto len = Length(i_vec); len > EPS)
      return i_vec / len;
    return {};
  }

  template <std::size_t NDim>
  inline bool IsUnit(const Vector<NDim>& i_vec)
  {
    auto len = Length(i_vec);
    return Equal(len, 1.0);
  }

  template <std::size_t NDim>
  inline double DotProduct(const Vector<NDim>& i_v1, const Vector<NDim>& i_v2)
  {
    return std::inner_product(i_v1.coords.begin(), i_v1.coords.end(), i_v2.coords.begin(), 0.0);
  }

  template <std::size_t NDim>
  inline Vector<NDim> CrossProduct(const Vector<NDim>& /*i_vec1*/, const Vector<NDim>& /*i_vec1*/)
  {
    static_assert(NDim == 3, "Cross product is only defined for 3D");
  }

  template <>
  Vector3D CrossProduct<3>(const Vector3D& i_vec1, const Vector3D& i_vec2)
  {
    return {
      i_vec1.coords[1] * i_vec2.coords[2] - i_vec1.coords[2] * i_vec2.coords[1],
      i_vec1.coords[2] * i_vec2.coords[0] - i_vec1.coords[0] * i_vec2.coords[2],
      i_vec1.coords[0] * i_vec2.coords[1] - i_vec1.coords[1] * i_vec2.coords[0]};
  }

  inline Vector2D To2D(const Vector3D& i_vec, const Vector3D& i_e1, const Vector3D& i_e2)
  {
    return {DotProduct(i_vec, i_e1), DotProduct(i_vec, i_e2)};
  }

  inline Vector2D Perp(const Vector2D& i_vec)
  {
    return {-i_vec.coords[1], i_vec.coords[0]};
  }

  // distance from a point to a plane defined by a normal vector and a point on the plane
  inline double SignedDistance(const Vector3D& i_n, const Point3D& i_p0, const Point3D& i_p)
  {
    assert(IsUnit(i_n));
    return DotProduct(i_n, i_p - i_p0);
  }

  template <std::size_t NPoints, std::size_t NDim>
  struct ConvexHull
  {
  public:
    using point_type = Point<NDim>;
    using vector_type = Point<NDim>;

  public:
    static constexpr std::size_t NCoords = NPoints*NDim;

    explicit ConvexHull(const double* ip_coords)
    {
      for(std::size_t i = 0; i < NPoints; ++i)
        std::copy_n(ip_coords + i * NDim, NDim, m_points[i].coords.data());
    }

    explicit ConvexHull(std::initializer_list<double> i_il)
    {
      assert(i_il.size() == NCoords);
      for(std::size_t i = 0; i < NPoints; ++i)
        std::copy_n(i_il.begin()+i*NDim, NDim, m_points[i].coords.data());
    }

    explicit ConvexHull(std::initializer_list<point_type> i_il)
    {
      assert(i_il.size() == NPoints);
      std::copy(i_il.begin(), i_il.end(), m_points.begin());
    }

    const point_type& operator[](std::size_t i_idx) const
    {
      return m_points[i_idx];
    }

    const std::array<point_type, NPoints>& Points() const
    {
      return m_points;
    }

  private:
    std::array<point_type, NPoints> m_points;
  };

  using Triangle2D = ConvexHull<3, 2>;
  using Triangle3D = ConvexHull<3, 3>;

  typename Triangle3D::vector_type Normal(const Triangle3D& i_triangle)
  {
    const auto& p1 = i_triangle[0];
    const auto& p2 = i_triangle[1];
    const auto& p3 = i_triangle[2];

    auto n = CrossProduct(p2 - p1, p3 - p1);
    return Normalize(n);
  }

  inline bool IsZeroArea(const Triangle3D& i_tr)
  {
		auto v1 = i_tr[1] - i_tr[0];
		auto v2 = i_tr[2] - i_tr[0];

		auto n = CrossProduct(v1, v2);
		return Equal(LengthSquared(n), 0.0);
	}

  struct Interval {double min, max;};
  inline bool Overlap(const Interval& i_i1, const Interval& i_i2)
  {
    // bool i1_valid = std::abs(i_i1.max - i_i1.min) > EPS;
    // bool i2_valid = std::abs(i_i2.max - i_i2.min) > EPS;

    // if (!i1_valid || !i2_valid)
    //   return false;

    return !(i_i1.min > i_i2.max || i_i2.min > i_i1.max);
  }

  // convex polygon intersection in 2D using edge normals as separating axes
  bool SAT(const Triangle2D& i_tr1, const Triangle2D& i_tr2)
  {
    auto project = [](const Triangle2D& i_tr, const Vector2D& i_axis) -> Interval
    {
      double min = DotProduct(i_axis, i_tr[0]), max = min;
      for(int i = 1; i < 3; ++i)
      {
        auto projection = DotProduct(i_axis, i_tr[i]);
        min = std::min(min, projection);
        max = std::max(max, projection);
      }
      return {min, max};
    };

    auto get_axes = [](const Triangle2D& i_tr)
    {
      return std::array{
          Perp(i_tr[1] - i_tr[0])
        , Perp(i_tr[2] - i_tr[1])
        , Perp(i_tr[0] - i_tr[2])};
    };

    auto axes1 = get_axes(i_tr1);
    auto axes2 = get_axes(i_tr2);

    for(const auto& axis: {axes1[0], axes1[1], axes1[2], axes2[0], axes2[1], axes2[2]})
    {
      auto interval1 = project(i_tr1, axis);
      auto interval2 = project(i_tr2, axis);

      if(!Overlap(interval1, interval2))
        return false;
    }

    return true;
  }

} // namespace _Details

namespace triangle_3d {

  bool intersect(double t1[9], double t2[9])
  {
    _Details::Triangle3D tr1{t1};
    _Details::Triangle3D tr2{t2};

    auto n1 = Normal(tr1);
    auto n2 = Normal(tr2);

    if(Length(CrossProduct(n1, n2)) > _Details::EPS) // planes intersect
    {
      // define a plane of tr1 -> Plane_tr1: n1, tr1[0]
      // compute distances of tr2 points to the plane of tr1
      std::array<double, 3> dist_tr2;
      for(int i = 0; i < 3; ++i)
        dist_tr2[i] = SignedDistance(n1, tr1[0], tr2[i]);

      // early reject
      if((dist_tr2[0] >  _Details::EPS && dist_tr2[1] >  _Details::EPS && dist_tr2[2] >  _Details::EPS) ||
         (dist_tr2[0] < -_Details::EPS && dist_tr2[1] < -_Details::EPS && dist_tr2[2] < -_Details::EPS))
        return false;

      // define a plane of tr2 -> Plane_tr2: n2, tr2[0]
      // compute distances of tr1 points to the plane of tr2
      std::array<double, 3> dist_tr1;
      for(int i = 0; i < 3; ++i)
        dist_tr1[i] = SignedDistance(n2, tr2[0], tr1[i]);

      // early reject
      if((dist_tr1[0] >  _Details::EPS && dist_tr1[1] >  _Details::EPS && dist_tr1[2] >  _Details::EPS) ||
         (dist_tr1[0] < -_Details::EPS && dist_tr1[1] < -_Details::EPS && dist_tr1[2] < -_Details::EPS))
        return false;

      // compute intersection line of two planes (direction)
      auto intersection_dir = Normalize(CrossProduct(n1, n2));

      // project triangles onto this line (reduces everything to 1D interval overlap)
      auto project = [&intersection_dir](const _Details::Triangle3D& i_tr, const std::array<double, 3>& i_dist) -> _Details::Interval
      {
        std::vector<double> proj;
        for (int i = 0; i < 3; ++i)
        {
          int j = (i+1)%3;
          if(_Details::Sign(i_dist[i])*_Details::Sign(i_dist[j]) < 0) // edge inersects the plane
          {
            auto t = i_dist[i]/(i_dist[i]-i_dist[j]);
            _Details::Point3D p{
              i_tr[i].coords[0] + t*(i_tr[j].coords[0]-i_tr[i].coords[0]),
              i_tr[i].coords[1] + t*(i_tr[j].coords[1]-i_tr[i].coords[1]),
              i_tr[i].coords[2] + t*(i_tr[j].coords[2]-i_tr[i].coords[2])};
            proj.push_back(DotProduct(p, intersection_dir));
          }
          else if (std::abs(i_dist[i]) < _Details::EPS) // vertex on the plane
          {
            proj.push_back(DotProduct(i_tr[i], intersection_dir));
          }
        }

        if(proj.size() < 2) // e.g. one vertex is exactly on the plane
          return {0.0, 0.0};

        return {std::min(proj[0], proj[1]), std::max(proj[0], proj[1])};
      };

      auto interval1 = project(tr1, dist_tr1);
      auto interval2 = project(tr2, dist_tr2);

      // check interval overlap (if intervals overlap, triangles intersect)
      return Overlap(interval1, interval2);
    }
    else
    {
      if(!_Details::Equal(DotProduct(n1, tr2[0] - tr1[0]), 0.0)) // planes are parallel, tr1 and tr2 don't intersect
        return false;

      // tr1 and tr2 are coplanar
      auto e1 = Normalize(tr1[1] - tr1[0]);      // first orthonormal basis
      auto e2 = Normalize(CrossProduct(n1, e1)); // second orthonormal basis

      // translate tr1 and tr2 from 3D to 2D
      _Details::Triangle2D tr1_2d{ To2D(tr1[0], e1, e2), To2D(tr1[1], e1, e2), To2D(tr1[2], e1, e2) };
      _Details::Triangle2D tr2_2d{ To2D(tr2[0], e1, e2), To2D(tr2[1], e1, e2), To2D(tr2[2], e1, e2) };

      return SAT(tr1_2d, tr2_2d);
    }
  }
} // namespace triangle_3d
