/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#ifndef _BASES_HPP_WAS_INCLUDED_
    #error "You must NOT include this file directly. Include bases.hpp"
#endif

#ifndef _BASES_ALL_HPP_
#define _BASES_ALL_HPP_

#include <vector>

#include "common/eigen.hpp"

namespace disk {




template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_scalar_basis<Mesh<T,3,Storage>, typename Mesh<T,3,Storage>::cell>
    : public priv::monomial_basis_bones<3>
{
    typedef Mesh<T,3,Storage>                       mesh_type;
    typedef typename mesh_type::scalar_type         scalar_type;
    typedef typename mesh_type::cell                cell_type;
    typedef priv::monomial_basis_bones<3>           base;

public:

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    eval_functions(const mesh_type& msh, const cell_type& cl,
                   const point<T,3>& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? max_degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto bar = barycenter(msh, cl);
        auto h = diameter(msh, cl);

        auto ep = (pt - bar)/h;

        Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> ret;
        ret.resize( eval_range.size(), 1 );

#ifdef POWER_CACHE
        power_cache<scalar_type, 3> pow(ep, this->max_degree()+1);
#endif
        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());
        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;

#ifdef POWER_CACHE
            auto vx = pow.x(m[0]);
            auto vy = pow.y(m[1]);
            auto vz = pow.z(m[2]);
#else
            auto vx = iexp_pow(ep.x(), m[0]);
            auto vy = iexp_pow(ep.y(), m[1]);
            auto vz = iexp_pow(ep.z(), m[2]);
#endif

            ret(i++) = vx * vy * vz;
        }

        return ret;
    }

    Eigen::Matrix<scalar_type, Eigen::Dynamic, 3>
    eval_gradients(const mesh_type& msh, const cell_type& cl,
                   const point<T,3>& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? max_degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto bar = barycenter(msh, cl);
        auto h = diameter(msh, cl);
        auto ih = 1./h;

        auto ep = (pt - bar)/h;

        Eigen::Matrix<scalar_type, Eigen::Dynamic, 3> ret;
        ret.resize( eval_range.size(), 3 );

#ifdef POWER_CACHE
        power_cache<scalar_type, 3> pow(ep, this->max_degree()+1);
#endif

        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());
        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;

#ifdef POWER_CACHE
            auto px = pow.x(m[0]);
            auto py = pow.y(m[1]);
            auto pz = pow.z(m[2]);

            auto dx = (m[0] == 0) ? 0 : m[0]*ih*pow.x(m[0]-1);
            auto dy = (m[1] == 0) ? 0 : m[1]*ih*pow.y(m[1]-1);
            auto dz = (m[2] == 0) ? 0 : m[2]*ih*pow.z(m[2]-1);
#else
            auto px = iexp_pow(ep.x(), m[0]);
            auto py = iexp_pow(ep.y(), m[1]);
            auto pz = iexp_pow(ep.z(), m[2]);

            auto dx = (m[0] == 0) ? 0 : m[0]*ih*iexp_pow(ep.x(), m[0]-1);
            auto dy = (m[1] == 0) ? 0 : m[1]*ih*iexp_pow(ep.y(), m[1]-1);
            auto dz = (m[2] == 0) ? 0 : m[2]*ih*iexp_pow(ep.z(), m[2]-1);
#endif

            ret(i,0) = dx * py * pz;
            ret(i,1) = dy * px * pz;
            ret(i,2) = dz * px * py;

            i++;
        }

        return ret;
    }
};


template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
point<T, 2>
map_point(const Mesh<T,3,Storage>& msh,
          const typename Mesh<T,3,Storage>::face& fc,
          const point<T,3>& pt)
{
    auto pts = points(msh, fc);
    auto diam = diameter(msh, fc);
    auto bar = barycenter(msh, fc);

    static_vector<T,3> v0;
    static_vector<T,3> v1;

    bool ok = false;

    size_t npts = pts.size();
    for (size_t i = 0; i < npts; i++)
    {
        size_t i0, i1;
        i0 = (i+1)%npts;
        i1 = (i-1)%npts;
        v0 = (pts[i0] - pts[i]).to_vector();
        v1 = (pts[i1] - pts[i]).to_vector();

        static_vector<T,3> v0n = v0/v0.norm();
        static_vector<T,3> v1n = v1/v1.norm();

        if ( v0n.dot(v1n) < 0.99 ) // we want at least 8 degrees angle
        {
            ok = true;
            break;
        }
    }

    if (!ok)
        throw std::invalid_argument("Degenerate polyhedron, cannot proceed");

    static_vector<T,3> e0 = v0 / v0.norm();
    static_vector<T,3> e1 = v1 - (v1.dot(v0) * v0) / (v0.dot(v0));
    e1 = e1 / e1.norm();

    static_vector<T,3> v = (pt-bar).to_vector();

    auto eta = v.dot(e0);
    auto xi = v.dot(e1);

    return point<T,2>({eta, xi});
}


template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_scalar_basis<Mesh<T,3,Storage>, typename Mesh<T,3,Storage>::face>
    : public priv::monomial_basis_bones<2>
{
    typedef Mesh<T,3,Storage>                   mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;
    typedef typename mesh_type::face            face_type;

    typedef priv::monomial_basis_bones<2>   base;

public:

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    eval_functions(const mesh_type& msh, const face_type& fc,
                   const point<T,3>& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? max_degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto ep = map_point(msh, fc, pt);

        Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> ret;
        ret.resize( eval_range.size(), 1 );

#ifdef POWER_CACHE
        power_cache<scalar_type, 2> pow(ep, this->max_degree()+1);
#endif

        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());
        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;

#ifdef POWER_CACHE
            auto vx = pow.x(m[0]);
            auto vy = pow.y(m[1]);
#else
            auto vx = iexp_pow(ep.x(), m[0]);
            auto vy = iexp_pow(ep.y(), m[1]);
#endif

            ret(i++) = vx * vy;
        }

        return ret;
    }
};





template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_scalar_basis<Mesh<T,2,Storage>, typename Mesh<T,2,Storage>::cell>
    : public priv::monomial_basis_bones<2>
{
    typedef Mesh<T,2,Storage>                       mesh_type;
    typedef typename mesh_type::scalar_type         scalar_type;
    typedef typename mesh_type::cell                cell_type;
    typedef typename mesh_type::point_type          point_type;
    typedef priv::monomial_basis_bones<2>           base;

public:

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    eval_functions(const mesh_type& msh, const cell_type& cl,
                   const point_type& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? max_degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto bar = barycenter(msh, cl);
        auto h = diameter(msh, cl);

        auto ep = (pt - bar)/h;

        Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> ret;
        ret.resize( eval_range.size(), 1 );

#ifdef POWER_CACHE
        power_cache<scalar_type, 2> pow(ep, this->max_degree()+1);
#endif

        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());
        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;
#ifdef POWER_CACHE
            auto vx = pow.x(m[0]);
            auto vy = pow.y(m[1]);
#else
            auto vx = iexp_pow(ep.x(), m[0]);
            auto vy = iexp_pow(ep.y(), m[1]);
#endif
            ret(i++) = vx * vy;
        }

        return ret;
    }

    Eigen::Matrix<scalar_type, Eigen::Dynamic, 2>
    eval_gradients(const mesh_type& msh, const cell_type& cl,
                   const point_type& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? max_degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto bar = barycenter(msh, cl);
        auto h = diameter(msh, cl);
        auto ih = 1./h;

        auto ep = (pt - bar)/h;

        Eigen::Matrix<scalar_type, Eigen::Dynamic, 2> ret;
        ret.resize( eval_range.size(), 2 );

#ifdef POWER_CACHE
        power_cache<scalar_type, 2> pow(ep, this->max_degree()+1);
#endif
        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());
        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;

#ifdef POWER_CACHE
            auto px = pow.x(m[0]);
            auto py = pow.y(m[1]);

            auto dx = (m[0] == 0) ? 0 : m[0]*ih*pow.x(m[0]-1);
            auto dy = (m[1] == 0) ? 0 : m[1]*ih*pow.y(m[1]-1);
#else
            auto px = iexp_pow(ep.x(), m[0]);
            auto py = iexp_pow(ep.y(), m[1]);

            auto dx = (m[0] == 0) ? 0 : m[0]*ih*iexp_pow(ep.x(), m[0]-1);
            auto dy = (m[1] == 0) ? 0 : m[1]*ih*iexp_pow(ep.y(), m[1]-1);
#endif
            ret(i,0) = dx * py;
            ret(i,1) = dy * px;
            i++;
        }

        return ret;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_scalar_basis<Mesh<T,2,Storage>, typename Mesh<T,2,Storage>::face>
    : public priv::monomial_basis_bones<1>
{
    typedef Mesh<T,2,Storage>                       mesh_type;
    typedef typename mesh_type::scalar_type         scalar_type;
    typedef typename mesh_type::point_type          point_type;
    typedef typename mesh_type::face                face_type;
    typedef priv::monomial_basis_bones<1>           base;

public:

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    eval_functions(const mesh_type& msh, const face_type& fc,
                   const point_type& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? max_degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto pts = points(msh, fc);
        auto bar = barycenter(msh, fc);
        auto h = diameter(msh, fc);
        auto v = (pts[1] - pts[0]).to_vector();
        auto t = (pt - bar).to_vector();
        T dot = v.dot(t);
        auto ep = point<T, 1>({dot/(h*h)});

        Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> ret;
        ret.resize( eval_range.size(), 1 );

        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());
        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;
            auto vx = iexp_pow(ep.x(), m[0]);
            ret(i++) = vx;
        }

        return ret;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_scalar_basis<Mesh<T,1,Storage>, typename Mesh<T,1,Storage>::cell>
    : public priv::monomial_basis_bones<1>
{
    typedef Mesh<T,1,Storage>                       mesh_type;
    typedef typename mesh_type::scalar_type         scalar_type;
    typedef typename mesh_type::cell                cell_type;
    typedef typename mesh_type::point_type          point_type;
    typedef priv::monomial_basis_bones<1>           base;

public:
    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    eval_functions(const mesh_type& msh, const cell_type& cl,
                   const point_type& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? max_degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto bar = barycenter(msh, cl);
        auto h = diameter(msh, cl);

        auto ep = (pt - bar)/h;

        Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> ret;
        ret.resize( eval_range.size(), 1 );

        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());
        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;
            auto vx = iexp_pow(ep.x(), m[0]);
            ret(i++) = vx;
        }

        return ret;
    }

    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    eval_gradients(const mesh_type& msh, const cell_type& cl,
                   const point_type& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? max_degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto bar = barycenter(msh, cl);
        auto h = diameter(msh, cl);

        auto ep = (pt - bar)/h;

        Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> ret;
        ret.resize( eval_range.size(), 1 );

        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());
        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;

            auto dx = (m[0] == 0) ? 0 : (m[0]/h)*iexp_pow(ep.x(), m[0]-1);

            ret(i++) = dx;
        }

        return ret;
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_scalar_basis<Mesh<T,1,Storage>, typename Mesh<T,1,Storage>::face>
    : public priv::monomial_basis_bones<0>
{
    typedef Mesh<T,1,Storage>                           mesh_type;
    typedef typename mesh_type::scalar_type         scalar_type;
    typedef typename mesh_type::point_type              point_type;
    typedef typename mesh_type::face                    face_type;
    typedef priv::monomial_basis_bones<0>               base;

public:

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    eval_functions(const mesh_type& msh, const face_type& fc, const point_type& pt) const
    {
        Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> ret(1);
        ret(0) = 1;
        return ret;
    }
};


template<template<typename, size_t, typename> class Mesh,
                  typename T, typename Storage>
std::vector<point<T,2>>
make_test_points(const Mesh<T,2,Storage>& msh,
                 const typename Mesh<T,2,Storage>::cell& cl)
{
    std::vector<point<T,2>> test_points;
    auto pts = points(msh, cl);
    auto bar = barycenter(msh, cl);

    test_points.insert(test_points.begin(), pts.begin(), pts.end()); //vertices
    test_points.push_back(bar); //barycenter

    //split in triangles and take barycenter
    for (size_t i = 0; i < pts.size(); i++)
        test_points.push_back( (pts[i] + pts[(i+1) % pts.size()] + bar) / 3. );

    return test_points;
}

template<template<typename, size_t, typename> class Mesh,
                  typename T, typename Storage>
std::vector<point<T,2>>
make_test_points(const Mesh<T,2,Storage>& msh,
                 const typename Mesh<T,2,Storage>::face& fc)
{
    std::vector<point<T,2>> test_points(5);
    auto pts = points(msh, fc);
    assert(pts.size() == 2); /* we are in 2D and we are dealing with a face,
                                which is an edge */

    test_points[0] = pts[0];
    test_points[1] = pts[1];
    test_points[2] = (pts[0] + pts[1]) / 2.;
    test_points[3] = (test_points[2] + pts[0]) / 2.;
    test_points[4] = (test_points[2] + pts[1]) / 2.;

    return test_points;
}

template<template<typename, size_t, typename> class Mesh,
                  typename T, typename Storage>
std::vector<point<T,1>>
make_test_points(const Mesh<T,1,Storage>& msh,
                 const typename Mesh<T,1,Storage>::cell& cl)
{
    std::vector<point<T,1>> test_points;
    auto pts = points(msh, cl);
    auto bar = barycenter(msh, cl);

    test_points.insert(test_points.begin(), pts.begin(), pts.end()); //vertices
    test_points.push_back(bar); //barycenter
    test_points.push_back((pts[0] + bar) / 2.);
    test_points.push_back((pts[1] + bar) / 2.);

    return test_points;
}


// vectorial basis function

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_vector_basis<Mesh<T,3,Storage>, typename Mesh<T,3,Storage>::cell>
: public priv::monomial_basis_bones<3,3>
{
   typedef Mesh<T,3,Storage>                       mesh_type;
   typedef typename mesh_type::scalar_type         scalar_type;
   typedef typename mesh_type::cell                cell_type;
   typedef priv::monomial_basis_bones<3,3>           base;

public:
   typedef static_vector<T,3>              function_value_type;
   typedef static_matrix<T,3,3>            gradient_value_type;

   scaled_monomial_vector_basis()
   : base(1)
   {}

   scaled_monomial_vector_basis(size_t degree)
      : base(degree)
      {}


   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const cell_type& cl, const point<T,3>& pt) const
   {
      auto bar = barycenter(msh, cl);
      auto h = measure(msh, cl);

      auto ep = (pt - bar)/h;

      std::vector<function_value_type> ret;
      ret.reserve( this->size() );

      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;
         auto vx = iexp_pow(ep.x(), m[0]);
         auto vy = iexp_pow(ep.y(), m[1]);
         auto vz = iexp_pow(ep.z(), m[2]);
         auto val = vx * vy * vz;
         ret.push_back( static_vector<T,3>({val,   0,   0}) );
         ret.push_back( static_vector<T,3>({  0, val,   0}) );
         ret.push_back( static_vector<T,3>({  0,   0, val}) );
      }

      return ret;
   }

   std::vector<gradient_value_type>
   eval_sgradients(const mesh_type& msh, const cell_type& cl, const point<T,3>& pt) const
   {
      auto bar = barycenter(msh, cl);
      auto h = measure(msh, cl);

      auto ep = (pt - bar)/h;

      std::vector<gradient_value_type> ret;
      ret.reserve( this->size() );

      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;

         auto px = iexp_pow(ep.x(), m[0]);
         auto py = iexp_pow(ep.y(), m[1]);
         auto pz = iexp_pow(ep.z(), m[2]);

         auto dx = (m[0] == 0) ? 0 : (m[0]/h)*iexp_pow(ep.x(), m[0]-1);
         auto dy = (m[1] == 0) ? 0 : (m[1]/h)*iexp_pow(ep.y(), m[1]-1);
         auto dz = (m[2] == 0) ? 0 : (m[2]/h)*iexp_pow(ep.z(), m[2]-1);

         auto val1 = dx * py * pz;
         auto val2 = px * dy * pz;
         auto val3 = px * py * dz;

         gradient_value_type sg;
         sg = gradient_value_type::Zero();
         sg(0,0) = val1;
         sg(0,1) = val2;
         sg(0,2) = val3;
         ret.push_back(0.5*(sg + sg.transpose()));

         sg = gradient_value_type::Zero();
         sg(1,0) = val1;
         sg(1,1) = val2;
         sg(1,2) = val3;
         ret.push_back(0.5*(sg + sg.transpose()));

         sg = gradient_value_type::Zero();
         sg(2,0) = val1;
         sg(2,1) = val2;
         sg(2,2) = val3;
         ret.push_back(0.5*(sg + sg.transpose()));
      }

      return ret;
   }



   std::vector<gradient_value_type>
   eval_ssgradients(const mesh_type& msh, const cell_type& cl, const point<T,3>& pt) const
   {
      auto bar = barycenter(msh, cl);
      auto h = measure(msh, cl);

      auto ep = (pt - bar)/h;

      std::vector<gradient_value_type> ret;
      ret.reserve( this->size() );

      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;

         auto px = iexp_pow(ep.x(), m[0]);
         auto py = iexp_pow(ep.y(), m[1]);
         auto pz = iexp_pow(ep.z(), m[2]);

         auto dx = (m[0] == 0) ? 0 : (m[0]/h)*iexp_pow(ep.x(), m[0]-1);
         auto dy = (m[1] == 0) ? 0 : (m[1]/h)*iexp_pow(ep.y(), m[1]-1);
         auto dz = (m[2] == 0) ? 0 : (m[2]/h)*iexp_pow(ep.z(), m[2]-1);

         auto val1 = dx * py * pz;
         auto val2 = px * dy * pz;
         auto val3 = px * py * dz;

         gradient_value_type sg;
         sg = gradient_value_type::Zero();
         sg(0,0) = val1;
         sg(0,1) = val2;
         sg(0,2) = val3;
         ret.push_back(0.5*(sg - sg.transpose()));

         sg = gradient_value_type::Zero();
         sg(1,0) = val1;
         sg(1,1) = val2;
         sg(1,2) = val3;
         ret.push_back(0.5*(sg - sg.transpose()));

         sg = gradient_value_type::Zero();
         sg(2,0) = val1;
         sg(2,1) = val2;
         sg(2,2) = val3;
         ret.push_back(0.5*(sg - sg.transpose()));
      }

      return ret;
   }

std::vector<gradient_value_type>
eval_ssgradients_const(const mesh_type& msh, const cell_type& cl) const
{

   std::vector<gradient_value_type> ret;
   ret.reserve(3);

   gradient_value_type ssg = gradient_value_type::Zero();
   ssg(0,1) = 1.0;
   ssg(1,0) = -1.0;
   ret.push_back(ssg);

   ssg = gradient_value_type::Zero();
   ssg(0,2) = 1.0;
   ssg(2,0) = -1.0;
   ret.push_back(ssg);


   ssg = gradient_value_type::Zero();
   ssg(1,2) = 1.0;
   ssg(2,1) = -1.0;
   ret.push_back(ssg);

   return ret;
}


   std::vector<gradient_value_type>
   eval_gradients(const mesh_type& msh, const cell_type& cl, const point<T,3>& pt) const
   {
      auto bar = barycenter(msh, cl);
      auto h = measure(msh, cl);

      auto ep = (pt - bar)/h;

      std::vector<gradient_value_type> ret;
      ret.reserve( this->size() );

      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;

         auto px = iexp_pow(ep.x(), m[0]);
         auto py = iexp_pow(ep.y(), m[1]);
         auto pz = iexp_pow(ep.z(), m[2]);

         auto dx = (m[0] == 0) ? 0 : (m[0]/h)*iexp_pow(ep.x(), m[0]-1);
         auto dy = (m[1] == 0) ? 0 : (m[1]/h)*iexp_pow(ep.y(), m[1]-1);
         auto dz = (m[2] == 0) ? 0 : (m[2]/h)*iexp_pow(ep.z(), m[2]-1);

         auto val1 = dx * py * pz;
         auto val2 = px * dy * pz;
         auto val3 = px * py * dz;

         gradient_value_type sg;
         sg = gradient_value_type::Zero();
         sg(0,0) = val1;
         sg(0,1) = val2;
         sg(0,2) = val3;
         ret.push_back(sg);

         sg = gradient_value_type::Zero();
         sg(1,0) = val1;
         sg(1,1) = val2;
         sg(1,2) = val3;
         ret.push_back(sg);

         sg = gradient_value_type::Zero();
         sg(2,0) = val1;
         sg(2,1) = val2;
         sg(2,2) = val3;
         ret.push_back(sg);
      }

      return ret;
   }
};



template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_vector_basis<Mesh<T,3,Storage>, typename Mesh<T,3,Storage>::face>
: public priv::monomial_basis_bones<2,3>
{
   typedef Mesh<T,3,Storage>                   mesh_type;
   typedef typename mesh_type::scalar_type     scalar_type;
   typedef typename mesh_type::face            face_type;

   typedef priv::monomial_basis_bones<2,3>   base;

public:
   typedef static_vector<T,3>                      function_value_type;

   scaled_monomial_vector_basis()
      : base(1)
      {}

   scaled_monomial_vector_basis(size_t degree)
      : base(degree)
      {}

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const face_type& fc, const point<T,3>& pt) const
   {
      //auto bar = barycenter(msh, fc);
      //auto h = measure(msh, fc);

      auto ep = map_point(msh, fc, pt);

      std::vector<function_value_type> ret;
      ret.reserve( this->size() );

      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;
         auto vx = iexp_pow(ep.x(), m[0]);
         auto vy = iexp_pow(ep.y(), m[1]);

         auto val = vx * vy;

         ret.push_back( static_vector<T,3>({val,   0,   0}) );
         ret.push_back( static_vector<T,3>({  0, val,   0}) );
         ret.push_back( static_vector<T,3>({  0,   0, val}) );
      }

      return ret;
   }
};




template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_vector_basis<Mesh<T,2,Storage>, typename Mesh<T,2,Storage>::cell>
: public priv::monomial_basis_bones<2,2>
{
   typedef Mesh<T,2,Storage>                       mesh_type;
   typedef typename mesh_type::scalar_type         scalar_type;
   typedef typename mesh_type::cell                cell_type;
   typedef typename mesh_type::point_type          point_type;
   typedef priv::monomial_basis_bones<2,2>           base;

public:
   typedef static_vector<T,2>              function_value_type;
   typedef static_matrix<T,2,2>            gradient_value_type;

   scaled_monomial_vector_basis()
      : base(1)
      {}

   scaled_monomial_vector_basis(size_t degree)
      : base(degree)
      {}

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const cell_type& cl, const point<T,2>& pt) const
   {
      auto bar = barycenter(msh, cl);
      auto h = measure(msh, cl);

      auto ep = (pt - bar)/h;

      std::vector<function_value_type> ret;
      ret.reserve( this->size() );

      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;
         auto vx = iexp_pow(ep.x(), m[0]);
         auto vy = iexp_pow(ep.y(), m[1]);
         auto val = vx * vy;
         ret.push_back( static_vector<T,2>({val,   0}) );
         ret.push_back( static_vector<T,2>({  0, val}) );
      }

      return ret;
   }

   std::vector<gradient_value_type>
   eval_sgradients(const mesh_type& msh, const cell_type& cl, const point<T,2>& pt) const
   {
      auto bar = barycenter(msh, cl);
      auto h = measure(msh, cl);

      auto ep = (pt - bar)/h;

      std::vector<gradient_value_type> ret;
      ret.reserve( this->size() );

      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;

         auto px = iexp_pow(ep.x(), m[0]);
         auto py = iexp_pow(ep.y(), m[1]);

         auto dx = (m[0] == 0) ? 0 : (m[0]/h)*iexp_pow(ep.x(), m[0]-1);
         auto dy = (m[1] == 0) ? 0 : (m[1]/h)*iexp_pow(ep.y(), m[1]-1);

         auto val1 = dx * py;
         auto val2 = px * dy;

         gradient_value_type sg;
         sg = gradient_value_type::Zero();
         sg(0,0) = val1 ;
         sg(0,1) = val2 ;
         ret.push_back(0.5*(sg + sg.transpose()));

         sg = gradient_value_type::Zero();
         sg(1,0) = val1 ;
         sg(1,1) = val2 ;
         ret.push_back(0.5*(sg + sg.transpose()));

      }

      return ret;
   }

   std::vector<gradient_value_type>
   eval_ssgradients(const mesh_type& msh, const cell_type& cl, const point<T,2>& pt) const
   {
      auto bar = barycenter(msh, cl);
      auto h = measure(msh, cl);

      auto ep = (pt - bar)/h;

      std::vector<gradient_value_type> ret;
      ret.reserve( this->size() );

      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;

         auto px = iexp_pow(ep.x(), m[0]);
         auto py = iexp_pow(ep.y(), m[1]);

         auto dx = (m[0] == 0) ? 0 : (m[0]/h)*iexp_pow(ep.x(), m[0]-1);
         auto dy = (m[1] == 0) ? 0 : (m[1]/h)*iexp_pow(ep.y(), m[1]-1);

         auto val1 = dx * py;
         auto val2 = px * dy;

         gradient_value_type sg;
         sg = gradient_value_type::Zero();
         sg(0,0) = val1 ;
         sg(0,1) = val2 ;
         ret.push_back(0.5*(sg - sg.transpose()));

         sg = gradient_value_type::Zero();
         sg(1,0) = val1 ;
         sg(1,1) = val2 ;
         ret.push_back(0.5*(sg - sg.transpose()));

      }

      return ret;
   }

   std::vector<gradient_value_type>
   eval_gradients(const mesh_type& msh, const cell_type& cl, const point<T,2>& pt) const
   {
      auto bar = barycenter(msh, cl);
      auto h = measure(msh, cl);

      auto ep = (pt - bar)/h;

      std::vector<gradient_value_type> ret;
      ret.reserve( this->size() );

      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;

         auto px = iexp_pow(ep.x(), m[0]);
         auto py = iexp_pow(ep.y(), m[1]);

         auto dx = (m[0] == 0) ? 0 : (m[0]/h)*iexp_pow(ep.x(), m[0]-1);
         auto dy = (m[1] == 0) ? 0 : (m[1]/h)*iexp_pow(ep.y(), m[1]-1);

         auto val1 = dx * py;
         auto val2 = px * dy;

         gradient_value_type sg;
         sg = gradient_value_type::Zero();
         sg(0,0) = val1 ;
         sg(0,1) = val2 ;
         ret.push_back(sg);

         sg = gradient_value_type::Zero();
         sg(1,0) = val1 ;
         sg(1,1) = val2 ;
         ret.push_back(sg);

      }

      return ret;
   }

std::vector<gradient_value_type>
eval_ssgradients_const(const mesh_type& msh, const cell_type& cl) const
{

    std::vector<gradient_value_type> ret;
    ret.reserve(1);

    gradient_value_type ssg = gradient_value_type::Zero();
    ssg(0,1) = 1.0;
    ssg(1,0) = -1.0;
    ret.push_back(ssg);

    return ret;
}
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_vector_basis<Mesh<T,2,Storage>, typename Mesh<T,2,Storage>::face>
: public priv::monomial_basis_bones<1,2>
{
   typedef Mesh<T,2,Storage>                       mesh_type;
   typedef typename mesh_type::scalar_type         scalar_type;
   typedef typename mesh_type::point_type          point_type;
   typedef typename mesh_type::face                face_type;
   typedef priv::monomial_basis_bones<1,2>           base;

public:
   typedef static_vector<T,2>                      function_value_type;

   scaled_monomial_vector_basis()
      : base(1)
      {}

   scaled_monomial_vector_basis(size_t degree)
      : base(degree)
      {}

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const face_type& fc, const point<T,2>& pt) const
   {
      auto pts = points(msh, fc);
      auto bar = barycenter(msh, fc);
      auto h = diameter(msh, fc);
      auto v = (pts[1] - pts[0]).to_vector();
      auto t = (pt - bar).to_vector();
      T dot = v.dot(t);
      auto ep = point<T, 1>({dot/(h*h)});

      std::vector<function_value_type> ret;
      ret.reserve( this->size() );

      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;
         auto vx = iexp_pow(ep.x(), m[0]);
         //std::cout << ep.x() <<" " << m[0] << " " << vx << '\n';

         auto val = vx ;

         ret.push_back( static_vector<T,2>({val,   0}) );
         ret.push_back( static_vector<T,2>({  0, val}) );
      }

      return ret;
   }
};



template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_vector_basis<Mesh<T,1,Storage>, typename Mesh<T,1,Storage>::cell>
: public priv::monomial_basis_bones<1>
{
   typedef Mesh<T,1,Storage>                       mesh_type;
   typedef typename mesh_type::scalar_type         scalar_type;
   typedef typename mesh_type::cell                cell_type;
   typedef typename mesh_type::point_type          point_type;
   typedef priv::monomial_basis_bones<1>           base;

public:
   typedef T              function_value_type;
   typedef T           gradient_value_type;

   scaled_monomial_vector_basis()
      : base(1)
      {}

   scaled_monomial_vector_basis(size_t degree)
      : base(degree)
      {}

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const cell_type& cl, const point_type& pt) const
   {
      auto bar = barycenter(msh, cl);
      auto h = diameter(msh, cl);

      auto ep = (pt - bar)/h;

      std::vector<function_value_type> ret;
      ret.reserve( this->size() );

      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;
         auto vx = iexp_pow(ep.x(), m[0]);
         function_value_type val;
         val = vx;

         ret.push_back(val);
      }

      return ret;
   }

   std::vector<gradient_value_type>
   eval_gradients(const mesh_type& msh, const cell_type& cl, const point_type& pt) const
   {
      auto bar = barycenter(msh, cl);
      auto h = diameter(msh, cl);

      auto ep = (pt - bar)/h;

      std::vector<gradient_value_type> ret;
      ret.reserve( this->size() );

      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;

         auto dx = (m[0] == 0) ? 0 : (m[0]/h)*iexp_pow(ep.x(), m[0]-1);

         gradient_value_type sg;
         sg = dx ;
         ret.push_back(sg);
      }

      return ret;
   }

   std::vector<gradient_value_type>
   eval_sgradients(const mesh_type& msh, const cell_type& cl, const point_type& pt) const
   {
      auto bar = barycenter(msh, cl);
      auto h = diameter(msh, cl);

      auto ep = (pt - bar)/h;

      std::vector<gradient_value_type> ret;
      ret.reserve( this->size() );

      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;

         auto dx = (m[0] == 0) ? 0 : (m[0]/h)*iexp_pow(ep.x(), m[0]-1);

         gradient_value_type sg;
         sg = dx ;
         ret.push_back(sg);
      }

      return ret;
   }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_vector_basis<Mesh<T,1,Storage>, typename Mesh<T,1,Storage>::face>
: public priv::monomial_basis_bones<0>
{
   typedef Mesh<T,1,Storage>                           mesh_type;
   typedef typename mesh_type::scalar_type         scalar_type;
   typedef typename mesh_type::point_type              point_type;
   typedef typename mesh_type::face                    face_type;
   typedef priv::monomial_basis_bones<0>               base;

public:
   typedef T              function_value_type;

   scaled_monomial_vector_basis()
      : base(1)
      {}

   scaled_monomial_vector_basis(size_t degree)
      : base(degree)
      {}

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const face_type& fc, const point_type& pt) const
   {
      std::vector<function_value_type> ret;
      ret.reserve( this->size() );

      assert(this->size() == 1);

      function_value_type val;
      val = T{1};

      ret.push_back(val);
      return ret;
   }
};


/////////////////////////
// matrix basis function
//////////////////////////

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_matrix_basis<Mesh<T,3,Storage>, typename Mesh<T,3,Storage>::cell>
: public priv::monomial_basis_bones<3,3>
{
   typedef Mesh<T,3,Storage>                       mesh_type;
   typedef typename mesh_type::scalar_type         scalar_type;
   typedef typename mesh_type::cell                cell_type;
   typedef priv::monomial_basis_bones<3,3>           base;

public:

   typedef static_matrix<T,3,3>            function_value_type;

   scaled_monomial_matrix_basis()
   : base(1)
   {}

   scaled_monomial_matrix_basis(size_t degree)
   : base(degree)
   {}

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const cell_type& cl, const point<T,3>& pt) const
   {
      auto bar = barycenter(msh, cl);
      auto h = diameter(msh, cl);

      auto ep = (pt - bar)/h;

      std::vector<function_value_type> ret;
      ret.reserve( 9*this->size());


      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;

         auto vx = iexp_pow(ep.x(), m[0]);
         auto vy = iexp_pow(ep.y(), m[1]);
         auto vz = iexp_pow(ep.z(), m[2]);

         auto val = vx * vy * vz;

         function_value_type fc;

         for( size_t i = 0; i < 3; i++){
            for( size_t j = 0; j < 3; j++){
               fc = function_value_type::Zero();
               fc(i,j) = val;
               ret.push_back(fc);
            }
         }
      }

      return ret;
   }
};


template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_matrix_basis<Mesh<T,3,Storage>, typename Mesh<T,3,Storage>::face>
: public priv::monomial_basis_bones<2,3>
{
   typedef Mesh<T,3,Storage>                   mesh_type;
   typedef typename mesh_type::scalar_type     scalar_type;
   typedef typename mesh_type::face            face_type;

   typedef priv::monomial_basis_bones<2,3>   base;

public:
   typedef static_matrix<T,3,3>                  function_value_type;

   scaled_monomial_matrix_basis()
   : base(1)
   {}

   scaled_monomial_matrix_basis(size_t degree)
   : base(degree)
   {}

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const face_type& fc, const point<T,3>& pt) const
   {
      auto ep = map_point(msh, fc, pt);

      std::vector<function_value_type> ret;
      ret.reserve(9* this->size());


      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;

         auto vx = iexp_pow(ep.x(), m[0]);
         auto vy = iexp_pow(ep.y(), m[1]);

         auto val = vx * vy;

         function_value_type fc;

         for( size_t i = 0; i < 3; i++){
            for( size_t j = 0; j < 3; j++){
               fc = function_value_type::Zero();
               fc(i,j) = val;
               ret.push_back(fc);
            }
         }
      }

      return ret;
   }
};





template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_matrix_basis<Mesh<T,2,Storage>, typename Mesh<T,2,Storage>::cell>
: public priv::monomial_basis_bones<2,2>
{
   typedef Mesh<T,2,Storage>                       mesh_type;
   typedef typename mesh_type::scalar_type         scalar_type;
   typedef typename mesh_type::cell                cell_type;
   typedef typename mesh_type::point_type          point_type;
   typedef priv::monomial_basis_bones<2,2>           base;

public:
   typedef static_matrix<T,2,2>                      function_value_type;

   scaled_monomial_matrix_basis()
   : base(1)
   {}

   scaled_monomial_matrix_basis(size_t degree)
   : base(degree)
   {}

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const cell_type& cl, const point_type& pt) const
   {
      auto bar = barycenter(msh, cl);
      auto h = diameter(msh, cl);

      auto ep = (pt - bar)/h;

      std::vector<function_value_type> ret;
      ret.reserve( 4*this->size());

      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;

         auto vx = iexp_pow(ep.x(), m[0]);
         auto vy = iexp_pow(ep.y(), m[1]);

         auto val = vx * vy;

         function_value_type fc;

         for( size_t i = 0; i < 2; i++){
            for( size_t j = 0; j < 2; j++){
               fc = function_value_type::Zero();
               fc(i,j) = val;
               ret.push_back(fc);
            }
         }
      }

      return ret;
   }

};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_matrix_basis<Mesh<T,2,Storage>, typename Mesh<T,2,Storage>::face>
: public priv::monomial_basis_bones<1,2>
{
   typedef Mesh<T,2,Storage>                       mesh_type;
   typedef typename mesh_type::scalar_type         scalar_type;
   typedef typename mesh_type::point_type          point_type;
   typedef typename mesh_type::face                face_type;
   typedef priv::monomial_basis_bones<1,2>           base;

public:
   typedef static_matrix<T,2,2>                      function_value_type;

   scaled_monomial_matrix_basis()
   : base(1)
   {}

   scaled_monomial_matrix_basis(size_t degree)
   : base(degree)
   {}

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const face_type& fc, const point_type& pt) const
   {
      auto pts = points(msh, fc);
      auto bar = barycenter(msh, fc);
      auto h = diameter(msh, fc);
      auto v = (pts[1] - pts[0]).to_vector();
      auto t = (pt - bar).to_vector();
      T dot = v.dot(t);
      auto ep = point<T, 1>({dot/(h*h)});

      std::vector<function_value_type> ret;
      ret.reserve(4 * this->size());


      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;
         auto vx = iexp_pow(ep.x(), m[0]);
         function_value_type fc;

         for( size_t i = 0; i < 2; i++){
            for( size_t j = 0; j < 2; j++){
               fc = function_value_type::Zero();
               fc(i,j) = vx;
               ret.push_back(fc);
            }
         }
      }

      return ret;
   }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_matrix_basis<Mesh<T,1,Storage>, typename Mesh<T,1,Storage>::cell>
: public priv::monomial_basis_bones<1>
{
   typedef Mesh<T,1,Storage>                       mesh_type;
   typedef typename mesh_type::scalar_type         scalar_type;
   typedef typename mesh_type::cell                cell_type;
   typedef typename mesh_type::point_type          point_type;
   typedef priv::monomial_basis_bones<1>           base;

public:
   typedef T              function_value_type;

   scaled_monomial_matrix_basis()
   : base(1)
   {}

   scaled_monomial_matrix_basis(size_t degree)
   : base(degree)
   {}

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const cell_type& cl,
                  const point_type& pt,
                  size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
   {
      maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? max_degree() : maxdeg;
      auto eval_range = range(mindeg, maxdeg);

      auto bar = barycenter(msh, cl);
      auto h = diameter(msh, cl);

      auto ep = (pt - bar)/h;

      std::vector<function_value_type> ret;
      ret.reserve( eval_range.size());

      auto begin = this->monomials_begin();
      std::advance(begin, eval_range.min());
      auto end = this->monomials_begin();
      std::advance(end, eval_range.max());
      for (auto itor = begin; itor != end; itor++)
      {
         auto m = *itor;
         auto vx = iexp_pow(ep.x(), m[0]);

         ret.push_back(vx);
      }

      return ret;
   }

};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_matrix_basis<Mesh<T,1,Storage>, typename Mesh<T,1,Storage>::face>
: public priv::monomial_basis_bones<0>
{
   typedef Mesh<T,1,Storage>                           mesh_type;
   typedef typename mesh_type::scalar_type         scalar_type;
   typedef typename mesh_type::point_type              point_type;
   typedef typename mesh_type::face                    face_type;
   typedef priv::monomial_basis_bones<0>               base;

public:
   typedef T              function_value_type;

   scaled_monomial_matrix_basis()
   : base(1)
   {}

   scaled_monomial_matrix_basis(size_t degree)
   : base(degree)
   {}

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const face_type& fc, const point_type& pt) const
   {
      assert(this->size() == 1);
      std::vector<function_value_type> ret(1, T{1});

      return ret;
   }
};
} // namespace disk


#endif /* _BASES_ALL_HPP_ */
