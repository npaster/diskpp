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
    typedef T                                       function_value_type;
    typedef static_vector<T,3>                      gradient_value_type;

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    scaled_monomial_scalar_basis(size_t degree, size_t computed_degree)
        : base(degree, computed_degree)
    {}

    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    eval_functions(const mesh_type& msh, const cell_type& cl,
                   const point<T,3>& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto bar = barycenter(msh, cl);
        auto h = diameter(msh, cl);

        auto ep = (pt - bar)/h;

        Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> ret;
        ret.resize( eval_range.size(), 1 );

#ifdef POWER_CACHE
        power_cache<scalar_type, 3> pow(ep, this->computed_degree()+1);
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
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto bar = barycenter(msh, cl);
        auto h = diameter(msh, cl);
        auto ih = 1./h;

        auto ep = (pt - bar)/h;

        Eigen::Matrix<scalar_type, Eigen::Dynamic, 3> ret;
        ret.resize( eval_range.size(), 3 );

#ifdef POWER_CACHE
        power_cache<scalar_type, 3> pow(ep, this->computed_degree()+1);
#endif

        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());
        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;
            gradient_value_type grad;

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
    typedef T                               function_value_type;

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    scaled_monomial_scalar_basis(size_t degree, size_t computed_degree)
        : base(degree, computed_degree)
    {}

    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    eval_functions(const mesh_type& msh, const face_type& fc,
                   const point<T,3>& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto ep = map_point(msh, fc, pt);

        Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> ret;
        ret.resize( eval_range.size(), 1 );

#ifdef POWER_CACHE
        power_cache<scalar_type, 2> pow(ep, this->computed_degree()+1);
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

#if 0
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class cell_basis<Mesh<T,2,Storage>>
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::cell_type       cell_type;
    typedef typename mesh_type::point_type      point_type;
    typedef typename mesh_type::scalar_type     scalar_type;

    point_type                      bar;
    scalar_type                     h, invh;

    std::shared_ptr<monomial_generator<2>>    monomials;

public:
    basis(std::shared_ptr<monomial_generator<2>> mg, const Mesh& msh, const Mesh::cell_type& cl)
        : monomials(mg)
    {
        bar     = barycenter(msh, cl);
        h       = diameter(msh, cl)/2.0;
        invh    = 1./h;
    }

    dynamic_vector<scalar_type>
    eval_functions(const point_type& pt) const
    {
        dynamic_vector<scalar_type> ret;
        ret.resize( monomials->size(), 1 );

        auto eval_point = (pt - bar)*invh;

        size_t i = 0;
        auto begin = monomials->begin();
        auto end = monomials->end();
        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;

            auto vx = iexp_pow(eval_point.x(), m[0]);
            auto vy = iexp_pow(eval_point.y(), m[1]);

            ret(i++) = vx * vy;
        }

        return ret;
    }
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class hho_space<Mesh<T,2,Storage>>
{
    typedef Mesh<T,2,Storage>       mesh_type;

    monomial_generator<2>           cell_mg;
    monomial_generator<1>           face_mg;

public:
    hho_space(const mesh_type& msh, size_t degree)
        : cell_mg(degree), face_mg(degree)
    {}

    hho_space(const mesh_type& msh, size_t cell_degree, size_t face_degree)
        : cell_mg(cell_degree), face_mg(face_degree)
    {}
};


#endif




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
    typedef T                                       function_value_type;
    typedef static_vector<T,2>                      gradient_value_type;

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    scaled_monomial_scalar_basis(size_t degree, size_t computed_degree)
        : base(degree, computed_degree)
    {}

    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    eval_functions(const mesh_type& msh, const cell_type& cl,
                   const point_type& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto bar = barycenter(msh, cl);
        auto h = diameter(msh, cl);

        auto ep = (pt - bar)/h;

        Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> ret;
        ret.resize( eval_range.size(), 1 );

#ifdef POWER_CACHE
        power_cache<scalar_type, 2> pow(ep, this->computed_degree()+1);
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
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto bar = barycenter(msh, cl);
        auto h = diameter(msh, cl);
        auto ih = 1./h;

        auto ep = (pt - bar)/h;

        Eigen::Matrix<scalar_type, Eigen::Dynamic, 2> ret;
        ret.resize( eval_range.size(), 2 );

#ifdef POWER_CACHE
        power_cache<scalar_type, 2> pow(ep, this->computed_degree()+1);
#endif
        size_t i = 0;
        auto begin = this->monomials_begin();
        std::advance(begin, eval_range.min());
        auto end = this->monomials_begin();
        std::advance(end, eval_range.max());
        for (auto itor = begin; itor != end; itor++)
        {
            auto m = *itor;
            gradient_value_type grad;

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
    typedef T                                       function_value_type;

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    scaled_monomial_scalar_basis(size_t degree, size_t computed_degree)
        : base(degree, computed_degree)
    {}

    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    eval_functions(const mesh_type& msh, const face_type& fc,
                   const point_type& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? degree() : maxdeg;
        auto eval_range = range(mindeg, maxdeg);

        auto pts = points(msh, fc);
        auto bar = barycenter(msh, fc);
        //auto h = (pts[1] - pts[0]).to_vector().norm()/2.0;
        auto h = diameter(msh, fc);
        auto v = (pts[1] - pts[0]).to_vector();
        //auto v = (bar - pts[0]).to_vector();
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

#if 0
template<typename T>
class legendre_bones
{
public:
    T legendre_p0(T x)  { return 1.0; }
    T legendre_p1(T x)  { return x; }
    T legendre_p2(T x)  { return 0.5*(3*x*x - 1); }
    T legendre_p3(T x)  { return 0.5*(5*x*x*x - 3*x); }
    T legendre_p4(T x)  { return 0.125*(35*x*x*x*x - 30*x*x + 3); }

    T legendre_d0(T x)  { return 0.0; }
    T legendre_d1(T x)  { return 1.0; }
    T legendre_d2(T x)  { return 3*x; }
    T legendre_d3(T x)  { return 0.5*(15*x*x - 3); }
    T legendre_d4(T x)  { return 0.125*(140*x*x*x - 60*x); }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class legendre_scalar_basis<Mesh<T,1,Storage>, typename Mesh<T,1,Storage>::cell>
    : public legendre_bones
{

};
#endif

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
    typedef T                                       function_value_type;
    typedef T                                       gradient_value_type;

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    scaled_monomial_scalar_basis(size_t degree, size_t computed_degree)
        : base(degree, computed_degree)
    {}

    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    eval_functions(const mesh_type& msh, const cell_type& cl,
                   const point_type& pt,
                   size_t mindeg = 0, size_t maxdeg = VERY_HIGH_DEGREE) const
    {
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? degree() : maxdeg;
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
        maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? degree() : maxdeg;
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
    typedef T                                           function_value_type;

    scaled_monomial_scalar_basis()
        : base(1)
    {}

    scaled_monomial_scalar_basis(size_t degree)
        : base(degree)
    {}

    scaled_monomial_scalar_basis(size_t degree, size_t computed_degree)
        : base(degree, computed_degree)
    {}

    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    eval_functions(const mesh_type& msh, const face_type& fc, const point_type& pt) const
    {
        Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> ret(1);
        ret(0) = 1;
        return ret;
    }
};


// vectorial basis functions

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_vector_basis<Mesh<T,3,Storage>, typename Mesh<T,3,Storage>::cell>
: public priv::monomial_basis_bones<3,3>
{
   typedef Mesh<T,3,Storage>                       mesh_type;
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
   typedef typename mesh_type::cell                cell_type;
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

};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class scaled_monomial_vector_basis<Mesh<T,2,Storage>, typename Mesh<T,2,Storage>::face>
: public priv::monomial_basis_bones<1,2>
{
   typedef Mesh<T,2,Storage>                       mesh_type;
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
   typedef typename mesh_type::cell                cell_type;
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
   eval_functions(const mesh_type& msh, const cell_type& cl, const point<T,1>& pt) const
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
   eval_gradients(const mesh_type& msh, const cell_type& cl, const point<T,1>& pt) const
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
   eval_sgradients(const mesh_type& msh, const cell_type& cl, const point<T,1>& pt) const
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
   eval_functions(const mesh_type& msh, const face_type& fc, const point<T,1>& pt) const
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


   //overloaded function


   size_t degree_index(const size_t degree) const {
      return 3 * base::degree_index(degree);
   }

   size_t size() const {
      return 3 * base::size();
   }

   size_t computed_size() const {
      return 3 * base::computed_size();
   }

   dof_range range(const size_t min_degree, const size_t max_degree) const
   {
      const dof_range vector_range = base::range(min_degree, max_degree);

      return dof_range(3*vector_range.min(), 3*vector_range.max());
   }

   dof_range range() const {
      return dof_range(0, size());
   }

   dof_range computed_range() const {
      return dof_range(0, computed_size());
   }

   //evaluation

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const cell_type& cl, const point<T,3>& pt) const
   {
      auto bar = barycenter(msh, cl);
      auto h = diameter(msh, cl);

      auto ep = (pt - bar)/h;

      std::vector<function_value_type> ret;
      ret.reserve( 3*this->size());


      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;

         auto vx = iexp_pow(ep.x(), m[0]);
         auto vy = iexp_pow(ep.y(), m[1]);
         auto vz = iexp_pow(ep.z(), m[2]);

         auto val = vx * vy * vz;

         function_value_type fc;

         for( size_t j = 0; j < 3; j++){
            for( size_t i = 0; i < 3; i++){
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

   //overloaded function


   size_t degree_index(const size_t degree) const {
      return 3 * base::degree_index(degree);
   }

   size_t size() const {
      return 3 * base::size();
   }

   size_t computed_size() const {
      return 3 * base::computed_size();
   }

   dof_range range(const size_t min_degree, const size_t max_degree) const
   {
      const dof_range vector_range = base::range(min_degree, max_degree);

      return dof_range(3*vector_range.min(), 3*vector_range.max());
   }

   dof_range range() const {
      return dof_range(0, size());
   }

   dof_range computed_range() const {
      return dof_range(0, computed_size());
   }

   //evaluation

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const face_type& fc, const point<T,3>& pt) const
   {
      auto ep = map_point(msh, fc, pt);

      std::vector<function_value_type> ret;
      ret.reserve(3* this->size());


      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;

         auto vx = iexp_pow(ep.x(), m[0]);
         auto vy = iexp_pow(ep.y(), m[1]);

         auto val = vx * vy;

         function_value_type fc;

         for( size_t j = 0; j < 3; j++){
            for( size_t i = 0; i < 3; i++){
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

   //overloaded function


   size_t degree_index(const size_t degree) const {
      return 2 * base::degree_index(degree);
   }

   size_t size() const {
      return 2 * base::size();
   }

   size_t computed_size() const {
      return 2 * base::computed_size();
   }

   dof_range range(const size_t min_degree, const size_t max_degree) const
   {
      const dof_range vector_range = base::range(min_degree, max_degree);

      return dof_range(2*vector_range.min(), 2*vector_range.max());
   }

   dof_range range() const {
      return dof_range(0, size());
   }

   dof_range computed_range() const {
      return dof_range(0, computed_size());
   }

   //evaluation

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const cell_type& cl, const point_type& pt) const
   {
      auto bar = barycenter(msh, cl);
      auto h = diameter(msh, cl);

      auto ep = (pt - bar)/h;

      std::vector<function_value_type> ret;
      ret.reserve( 2*this->size());


      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;

         auto vx = iexp_pow(ep.x(), m[0]);
         auto vy = iexp_pow(ep.y(), m[1]);

         auto val = vx * vy;

         function_value_type fc;

         for( size_t j = 0; j < 2; j++){
            for( size_t i = 0; i < 2; i++){
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

   //overloaded function


   size_t degree_index(const size_t degree) const {
      return 2 * base::degree_index(degree);
   }

   size_t size() const {
      return 2 * base::size();
   }

   size_t computed_size() const {
      return 2 * base::computed_size();
   }

   dof_range range(const size_t min_degree, const size_t max_degree) const
   {
      const dof_range vector_range = base::range(min_degree, max_degree);

      return dof_range(2*vector_range.min(), 2*vector_range.max());
   }

   dof_range range() const {
      return dof_range(0, size());
   }

   dof_range computed_range() const {
      return dof_range(0, computed_size());
   }

   //evaluation

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
      ret.reserve(2 * this->size());


      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         auto m = *itor;
         auto vx = iexp_pow(ep.x(), m[0]);
         function_value_type fc;

         for( size_t j = 0; j < 2; j++){
            for( size_t i = 0; i < 2; i++){
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
      maxdeg = (maxdeg == VERY_HIGH_DEGREE) ? degree() : maxdeg;
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






/////////////////////////
//// RT elements ////////
/////////////////////////

// RT^{k+1}(T; R^d) = P^k(T;R^d) + x * P^{k,H}(T;R)

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class Raviart_Thomas_vector_basis<Mesh<T,3,Storage>, typename Mesh<T,3,Storage>::cell>
: public priv::monomial_basis_bones<3,3>
{
   typedef Mesh<T,3,Storage>                       mesh_type;
   typedef typename mesh_type::cell                cell_type;
   typedef priv::monomial_basis_bones<3,3>           base;
   typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    scaled_monomial_scalar_cell_basis_type;
   typedef disk::scaled_monomial_vector_basis<mesh_type, cell_type>    scaled_monomial_vector_cell_basis_type;
   scaled_monomial_scalar_cell_basis_type smscb;
   scaled_monomial_vector_cell_basis_type smvcb;

public:
   typedef static_vector<T,3>              function_value_type;
   typedef static_matrix<T,3,3>            gradient_value_type;

   Raviart_Thomas_vector_basis()
   : base(1)
   {
      smscb = scaled_monomial_scalar_cell_basis_type(this->degree()-1);
      smvcb = scaled_monomial_vector_cell_basis_type(this->degree()-1);
   }

   Raviart_Thomas_vector_basis(size_t degree)
   : base(degree)
   {
      smscb = scaled_monomial_scalar_cell_basis_type(this->degree()-1);
      smvcb = scaled_monomial_vector_cell_basis_type(this->degree()-1);
   }

   //overloaded function

   size_t size() const {
      const size_t degree = this->degree() - 1;
      return (degree+1) * (degree+2) * (degree+4)/2;
   }

   size_t computed_size() const {
      const size_t computed_degree = this->computed_degree() - 1;
      return (computed_degree+1) * (computed_degree+2) * (computed_degree+4)/2;
   }

   dof_range range(const size_t min_degree, const size_t max_degree) const
   {
      const size_t degree = this->degree() - 1;
      if(max_degree <= degree){
         return base::range(min_degree, max_degree);
      }
      else if(max_degree == this->degree()){
         if(min_degree < this->degree()){
            const dof_range vector_range = base::range(min_degree, max_degree);
            return dof_range(vector_range.min(), (degree+1) * (degree+2) * (degree+4)/2);
         }
         else{
            return dof_range(base::size()+1, (degree+1) * (degree+2) * (degree+4)/2);
         }
      }
      else{
         return base::range(min_degree, max_degree);
      }
   }

   dof_range range() const {
      return dof_range(0, size());
   }

   dof_range computed_range() const {
      return dof_range(0, computed_size());
   }

   //evaluation


   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const cell_type& cl, const point<T,3>& pt) const
   {
      std::vector<function_value_type> ret;
      const size_t poly_degree = (this->degree() - 1);
      ret.reserve( (poly_degree+1) * (poly_degree+2) * (poly_degree+4)/2 );

      // basis functions of P^{poly_degree}_d(T; R^d)
      const auto smvcb_phi = smvcb.eval_functions(msh, cl, pt);

      for (size_t i = 0; i < smvcb_phi.size(); i++){
         ret.push_back(smvcb_phi[i]);
      }

      // basis functions of P^{poly_degree},H_d(T; R)
      const auto smscb_phi = smscb.eval_functions(msh, cl, pt, poly_degree, poly_degree);
      const size_t dim_PHk = (poly_degree + 1) * (poly_degree + 2) / 2;
      assert( smscb_phi.rows() == dim_PHk);

      // basis functions of x*P^k,H_d(T; R)
      for(size_t i = 0; i < dim_PHk; i++){
         const auto vx = pt.x() * smscb_phi(i);
         const auto vy = pt.y() * smscb_phi(i);
         const auto vz = pt.z() * smscb_phi(i);
         ret.push_back( static_vector<T,3>({vx,  vy,  vz}) );
      }

      return ret;
   }
};



template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class Raviart_Thomas_vector_basis<Mesh<T,3,Storage>, typename Mesh<T,3,Storage>::face>
: public priv::monomial_basis_bones<2,3>
{
   typedef Mesh<T,3,Storage>                   mesh_type;
   typedef typename mesh_type::face            face_type;

   typedef priv::monomial_basis_bones<2,3>   base;
   typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    scaled_monomial_scalar_face_basis_type;
   typedef disk::scaled_monomial_vector_basis<mesh_type, face_type>    scaled_monomial_vector_face_basis_type;
   scaled_monomial_scalar_face_basis_type smsfb;
   scaled_monomial_vector_face_basis_type smvfb;

public:
   typedef static_vector<T,3>                      function_value_type;

   Raviart_Thomas_vector_basis()
   : base(1)
   {
      smsfb = scaled_monomial_scalar_face_basis_type(this->degree()-1);
      smvfb = scaled_monomial_vector_face_basis_type(this->degree()-1);
   }

   Raviart_Thomas_vector_basis(size_t degree)
   : base(degree)
   {
      smsfb = scaled_monomial_scalar_face_basis_type(this->degree()-1);
      smvfb = scaled_monomial_vector_face_basis_type(this->degree()-1);
   }

   //overloaded function

   size_t size() const {
      const size_t degree = this->degree() - 1;
      return (degree+1) * (degree+2) * (degree+4)/2;
   }

   size_t computed_size() const {
      const size_t computed_degree = this->computed_degree() - 1;
      return (computed_degree * (computed_degree+2) * (computed_degree+4)/2);
   }

   dof_range range(const size_t min_degree, const size_t max_degree) const
   {
      const size_t degree = this->degree() - 1;
      if(max_degree <= degree){
         return base::range(min_degree, max_degree);
      }
      else if(max_degree == this->degree()){
         if(min_degree < this->degree()){
            const dof_range vector_range = base::range(min_degree, max_degree);
            return dof_range(vector_range.min(), (degree+1) * (degree+2) * (degree+4)/2);
         }
         else{
            return dof_range(base::size()+1, (degree+1) * (degree+2) * (degree+4)/2);
         }
      }
      else{
         return base::range(min_degree, max_degree);
      }
   }

   dof_range range() const {
      return dof_range(0, size());
   }

   dof_range computed_range() const {
      return dof_range(0, computed_size());
   }

   //evaluation

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const face_type& fc, const point<T,3>& pt) const
   {
      std::vector<function_value_type> ret;
      const size_t poly_degree = (this->degree() - 1);
      ret.reserve( (poly_degree+1) * (poly_degree+2) * (poly_degree+4)/2 );

      // basis functions of P^{poly_degree}_d(T; R^d)
      const auto smvfb_phi = smvfb.eval_functions(msh, fc, pt);

      for (size_t i = 0; i < smvfb_phi.size(); i++){
         ret.push_back(smvfb_phi[i]);
      }

      // basis functions of P^k,H_d(T; R)
      const auto smsfb_phi = smsfb.eval_functions(msh, fc, pt, poly_degree, poly_degree);
      const size_t dim_PHk = (poly_degree + 1) * (poly_degree + 2) / 2;
      assert( smsfb_phi.rows() == dim_PHk);

      // basis functions of x*P^k,H_d(T; R)
      for(size_t i = 0; i < dim_PHk; i++){
         const auto vx = pt.x() * smsfb_phi(i);
         const auto vy = pt.y() * smsfb_phi(i);
         const auto vz = pt.z() * smsfb_phi(i);
         ret.push_back( static_vector<T,3>({vx,  vy,  vz}) );
      }

      return ret;
   }
};


template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class Raviart_Thomas_vector_basis<Mesh<T,2,Storage>, typename Mesh<T,2,Storage>::cell>
: public priv::monomial_basis_bones<2,2>
{
   typedef Mesh<T,2,Storage>                       mesh_type;
   typedef typename mesh_type::cell                cell_type;
   typedef priv::monomial_basis_bones<2,2>           base;

   typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    scaled_monomial_scalar_cell_basis_type;
   typedef disk::scaled_monomial_vector_basis<mesh_type, cell_type>    scaled_monomial_vector_cell_basis_type;
   scaled_monomial_scalar_cell_basis_type smscb;
   scaled_monomial_vector_cell_basis_type smvcb;

public:
   typedef static_vector<T,2>              function_value_type;
   typedef static_matrix<T,2,2>            gradient_value_type;

   Raviart_Thomas_vector_basis()
   : base(1)
   {
      smscb = scaled_monomial_scalar_cell_basis_type(this->degree()-1);
      smvcb = scaled_monomial_vector_cell_basis_type(this->degree()-1);
   }

   Raviart_Thomas_vector_basis(size_t degree)
   : base(degree)
   {
      smscb = scaled_monomial_scalar_cell_basis_type(this->degree()-1);
      smvcb = scaled_monomial_vector_cell_basis_type(this->degree()-1);
   }

   //overloaded function

   size_t size() const {
      const size_t degree = this->degree() - 1;
      return (degree+1) * (degree+3);
   }

   size_t computed_size() const {
      const size_t computed_degree = this->computed_degree() - 1;
      return (computed_degree+1) * (computed_degree+3);
   }

   dof_range range(const size_t min_degree, const size_t max_degree) const
   {
      const size_t degree = this->degree() - 1;
      if(max_degree <= degree){
         return base::range(min_degree, max_degree);
      }
      else if(max_degree == this->degree()){
         if(min_degree < this->degree()){
            const dof_range vector_range = base::range(min_degree, max_degree);
            return dof_range(vector_range.min(), (degree+1) * (degree+3));
         }
         else{
            return dof_range(base::size()+1, (degree+1) * (degree+3));
         }
      }
      else{
         return base::range(min_degree, max_degree);
      }
   }

   dof_range range() const {
      return dof_range(0, size());
   }

   dof_range computed_range() const {
      return dof_range(0, computed_size());
   }

   //evaluation

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const cell_type& cl, const point<T,2>& pt) const
   {
      std::vector<function_value_type> ret;
      const size_t poly_degree = this->degree()-1;
      ret.reserve( (poly_degree+1) * (poly_degree+3));


      // basis functions of P^k_d(T; R^d)
      const auto smvcb_phi = smvcb.eval_functions(msh, cl, pt);

      for (size_t i = 0; i < smvcb_phi.size(); i++){
         ret.push_back(smvcb_phi[i]);
      }

      // basis functions of P^k,H_d(T; R)
      const auto smscb_phi = smscb.eval_functions(msh, cl, pt, poly_degree, poly_degree);
      const size_t dim_PHk = (poly_degree + 1);
      assert( smscb_phi.rows() == dim_PHk);

      // basis functions of x*P^k,H_d(T; R)
      for(size_t i = 0; i < dim_PHk; i++){
         const auto vx = pt.x() * smscb_phi(i);
         const auto vy = pt.y() * smscb_phi(i);
         ret.push_back( static_vector<T,2>({vx,  vy}) );
      }

      return ret;
   }

};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class Raviart_Thomas_vector_basis<Mesh<T,2,Storage>, typename Mesh<T,2,Storage>::face>
: public priv::monomial_basis_bones<1,2>
{
   typedef Mesh<T,2,Storage>                       mesh_type;
   typedef typename mesh_type::face                face_type;
   typedef priv::monomial_basis_bones<1,2>           base;
   typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    scaled_monomial_scalar_face_basis_type;
   typedef disk::scaled_monomial_vector_basis<mesh_type, face_type>    scaled_monomial_vector_face_basis_type;
   scaled_monomial_scalar_face_basis_type smsfb;
   scaled_monomial_vector_face_basis_type smvfb;

public:
   typedef static_vector<T,2>                      function_value_type;

   Raviart_Thomas_vector_basis()
   : base(1)
   {
      smsfb = scaled_monomial_scalar_face_basis_type(this->degree()-1);
      smvfb = scaled_monomial_vector_face_basis_type(this->degree()-1);
   }

   Raviart_Thomas_vector_basis(size_t degree)
   : base(degree)
   {
      smsfb = scaled_monomial_scalar_face_basis_type(this->degree()-1);
      smvfb = scaled_monomial_vector_face_basis_type(this->degree()-1);
   }


   //overloaded function

   size_t size() const {
      const size_t degree = this->degree() - 1;
      return (degree+1) * (degree+3);
   }

   size_t computed_size() const {
      const size_t computed_degree = this->computed_degree() - 1;
      return (computed_degree+1) * (computed_degree+3);
   }

   dof_range range(const size_t min_degree, const size_t max_degree) const
   {
      const size_t degree = this->degree() - 1;
      if(max_degree <= degree){
         return base::range(min_degree, max_degree);
      }
      else if(max_degree == this->degree()){
         if(min_degree < this->degree()){
            const dof_range vector_range = base::range(min_degree, max_degree);
            return dof_range(vector_range.min(), (degree+1) * (degree+3));
         }
         else{
            return dof_range(base::size()+1, (degree+1) * (degree+3));
         }
      }
      else{
         return base::range(min_degree, max_degree);
      }
   }

   dof_range range() const {
      return dof_range(0, size());
   }

   dof_range computed_range() const {
      return dof_range(0, computed_size());
   }

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const face_type& fc, const point<T,2>& pt) const
   {
      std::vector<function_value_type> ret;
      const size_t poly_degree = this->degree()-1;
      ret.reserve( (poly_degree+1) * (poly_degree+3));


      // basis functions of P^k_d(T; R^d)
      const auto smvfb_phi = smvfb.eval_functions(msh, fc, pt);

      for (size_t i = 0; i < smvfb_phi.size(); i++){
         ret.push_back(smvfb_phi[i]);
      }

      // basis functions of P^k,H_d(T; R)
      const auto smsfb_phi = smsfb.eval_functions(msh, fc, pt, poly_degree, poly_degree);
      const size_t dim_PHk = poly_degree + 1;
      assert( smsfb_phi.rows() == dim_PHk);

      // basis functions of x*P^k,H_d(T; R)
      for(size_t i = 0; i < dim_PHk; i++){
         const auto vx = pt.x() * smsfb_phi(i);
         const auto vy = pt.y() * smsfb_phi(i);
         ret.push_back( static_vector<T,2>({vx,  vy}) );
      }

      return ret;
   }
};



template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class Raviart_Thomas_matrix_basis<Mesh<T,3,Storage>, typename Mesh<T,3,Storage>::cell>
: public priv::monomial_basis_bones<3,3>
{
   typedef Mesh<T,3,Storage>                       mesh_type;
   typedef typename mesh_type::cell                cell_type;
   typedef priv::monomial_basis_bones<3,3>           base;
   typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    scaled_monomial_scalar_cell_basis_type;

   scaled_monomial_scalar_cell_basis_type smscb;

public:
   typedef static_matrix<T,3,3>              function_value_type;

   Raviart_Thomas_matrix_basis()
   : base(1)
   {
      smscb = scaled_monomial_scalar_cell_basis_type(this->degree());
   }

   Raviart_Thomas_matrix_basis(size_t degree)
   : base(degree)
   {
      smscb = scaled_monomial_scalar_cell_basis_type(this->degree());
   }

   //overloaded function

   size_t size() const {
      return 3*(this->degree()+1) * (this->degree()+2) * (this->degree()+4)/2;
   }

   size_t computed_size() const {
      return 3*(this->computed_degree()+1) * (this->computed_degree()+2) * (this->computed_degree()+4)/2;
   }

   dof_range range(const size_t min_degree, const size_t max_degree) const
   {
      const dof_range vector_range = base::range(min_degree, max_degree);

      if(max_degree == computed_degree())
         return dof_range(3*vector_range.min(), 3*(vector_range.max() + (max_degree + 1) * (max_degree + 2) / 2));
      else
         return dof_range(3*vector_range.min(), 3*vector_range.max());
   }

   dof_range range() const {
      return dof_range(0, size());
   }

   dof_range computed_range() const {
      return dof_range(0, computed_size());
   }

   //evaluation

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const cell_type& cl, const point<T,3>& pt) const
   {
      const auto bar = barycenter(msh, cl);
      const auto h = measure(msh, cl);

      const auto ep = (pt - bar)/h;

      std::vector<function_value_type> ret;
      ret.reserve( 3*(this->degree()+1) * (this->degree()+2) * (this->degree()+4)/2 );

      // basis functions of P^k_d(T; R^d)
      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         const auto m = *itor;
         const auto vx = iexp_pow(ep.x(), m[0]);
         const auto vy = iexp_pow(ep.y(), m[1]);
         const auto vz = iexp_pow(ep.z(), m[2]);
         const auto val = vx * vy * vz;

         function_value_type fc;

         for( size_t j = 0; j < 3; j++){
            for( size_t i = 0; i < 3; i++){
               fc = function_value_type::Zero();
               fc(i,j) = val;
               ret.push_back(fc);
            }
         }
      }

      // basis functions of P^k,H_d(T; R)
      const auto smscb_phi = smscb.eval_functions(msh, cl, pt, this->degree(), this->degree());
      const size_t dim_PHk = (this->degree() + 1) * (this->degree() + 2) / 2;
      assert( smscb_phi.rows() == dim_PHk);

      // basis functions of P^k,H_d(T; R) * X
      for(size_t i = 0; i < dim_PHk; i++){
         const auto vx = pt.x() * smscb_phi(i);
         const auto vy = pt.y() * smscb_phi(i);
         const auto vz = pt.z() * smscb_phi(i);

         function_value_type fc;

         for(size_t k = 0; k < 3; k++){
            fc = function_value_type::Zero();

            fc(k,0) = vx;
            fc(k,1) = vy;
            fc(k,2) = vz;

            ret.push_back(fc);
         }
      }

      return ret;
   }
};



template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class Raviart_Thomas_matrix_basis<Mesh<T,3,Storage>, typename Mesh<T,3,Storage>::face>
: public priv::monomial_basis_bones<2,3>
{
   typedef Mesh<T,3,Storage>                   mesh_type;
   typedef typename mesh_type::face            face_type;

   typedef priv::monomial_basis_bones<2,3>   base;
   typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    scaled_monomial_scalar_face_basis_type;

   scaled_monomial_scalar_face_basis_type smsfb;

public:
   typedef static_matrix<T,3,3>                      function_value_type;

   Raviart_Thomas_matrix_basis()
   : base(1)
   {
      smsfb = scaled_monomial_scalar_face_basis_type(this->degree());
   }

   Raviart_Thomas_matrix_basis(size_t degree)
   : base(degree)
   {
      smsfb = scaled_monomial_scalar_face_basis_type(this->degree());
   }

   //overloaded function

   size_t size() const {
      return 3*(this->degree()+1) * (this->degree()+2) * (this->degree()+4)/2;
   }

   size_t computed_size() const {
      return 3*(this->computed_degree()+1) * (this->computed_degree()+2) * (this->computed_degree()+4)/2;
   }

   dof_range range(const size_t min_degree, const size_t max_degree) const
   {
      const dof_range vector_range = base::range(min_degree, max_degree);

      if(max_degree == computed_degree())
         return dof_range(3*vector_range.min(), 3*(vector_range.max() + (max_degree + 1) * (max_degree + 2) / 2));
      else
         return dof_range(3*vector_range.min(), 3*vector_range.max());
   }

   dof_range range() const {
      return dof_range(0, size());
   }

   dof_range computed_range() const {
      return dof_range(0, computed_size());
   }

   //evaluation

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const face_type& fc, const point<T,3>& pt) const
   {
      const auto ep = map_point(msh, fc, pt);

      std::vector<function_value_type> ret;
      ret.reserve(3*(this->degree()+1) * (this->degree()+2) * (this->degree()+4)/2 );

      // basis functions of P^k_d(T; R^d)
      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         const auto m = *itor;
         const auto vx = iexp_pow(ep.x(), m[0]);
         const auto vy = iexp_pow(ep.y(), m[1]);

         const auto val = vx * vy;

         function_value_type fc;

         for( size_t j = 0; j < 3; j++){
            for( size_t i = 0; i < 3; i++){
               fc = function_value_type::Zero();
               fc(i,j) = val;
               ret.push_back(fc);
            }
         }
      }

      // basis functions of P^k,H_d(T; R)
      const auto smsfb_phi = smsfb.eval_functions(msh, fc, pt, this->degree(), this->degree());
      const size_t dim_PHk = (this->degree() + 1) * (this->degree() + 2) / 2;
      assert( smsfb_phi.rows() == dim_PHk);

      // basis functions of P^k,H_d(T; R) * x
      for(size_t i = 0; i < dim_PHk; i++){
         const auto vx = pt.x() * smsfb_phi(i);
         const auto vy = pt.y() * smsfb_phi(i);
         const auto vz = pt.z() * smsfb_phi(i);

         function_value_type fc;

         for(size_t k = 0; k < 3; k++){
            fc = function_value_type::Zero();

            fc(k,0) = vx;
            fc(k,1) = vy;
            fc(k,2) = vz;

            ret.push_back(fc);
         }
      }

      return ret;
   }
};


template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class Raviart_Thomas_matrix_basis<Mesh<T,2,Storage>, typename Mesh<T,2,Storage>::cell>
: public priv::monomial_basis_bones<2,2>
{
   typedef Mesh<T,2,Storage>                       mesh_type;
   typedef typename mesh_type::cell                cell_type;
   typedef priv::monomial_basis_bones<2,2>           base;

   typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    scaled_monomial_scalar_cell_basis_type;

   scaled_monomial_scalar_cell_basis_type smscb;

public:
   typedef static_matrix<T,2,2>              function_value_type;

   Raviart_Thomas_matrix_basis()
   : base(1)
   {
      smscb = scaled_monomial_scalar_cell_basis_type(this->degree());
   }

   Raviart_Thomas_matrix_basis(size_t degree)
   : base(degree)
   {
      smscb = scaled_monomial_scalar_cell_basis_type(this->degree());
   }

   //overloaded function

   size_t size() const {
      return 2*(this->degree()+1) * (this->degree()+3);
   }

   size_t computed_size() const {
      return 2*(this->computed_degree()+1) * (this->computed_degree()+3);
   }

   dof_range range(const size_t min_degree, const size_t max_degree) const
   {
      const dof_range vector_range = base::range(min_degree, max_degree);

      if(max_degree == computed_degree())
         return dof_range(2*vector_range.min(), 2*(vector_range.max() + (max_degree + 1)));
      else
         return dof_range(2*vector_range.min(), 2*vector_range.max());
   }

   dof_range range() const {
      return dof_range(0, size());
   }

   dof_range computed_range() const {
      return dof_range(0, computed_size());
   }

   //evaluation

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const cell_type& cl, const point<T,2>& pt) const
   {
      auto bar = barycenter(msh, cl);
      auto h = measure(msh, cl);

      auto ep = (pt - bar)/h;

      std::vector<function_value_type> ret;
      ret.reserve( 2*(this->degree()+1) * (this->degree()+3));


      // basis functions of P^k_d(T; R^d)
      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         const auto m = *itor;
         const auto vx = iexp_pow(ep.x(), m[0]);
         const auto vy = iexp_pow(ep.y(), m[1]);
         const auto val = vx * vy;

         function_value_type fc;

         for( size_t j = 0; j < 2; j++){
            for( size_t i = 0; i < 2; i++){
               fc = function_value_type::Zero();
               fc(i,j) = val;
               ret.push_back(fc);
            }
         }
      }

      // basis functions of P^k,H_d(T; R)
      const auto smscb_phi = smscb.eval_functions(msh, cl, pt, this->degree(), this->degree());
      const size_t dim_PHk = (this->degree() + 1);
      assert( smscb_phi.rows() == dim_PHk);

      // basis functions of x*P^k,H_d(T; R)
      for(size_t i = 0; i < dim_PHk; i++){
         const auto vx = pt.x() * smscb_phi(i);
         const auto vy = pt.y() * smscb_phi(i);

         function_value_type fc;

         for(size_t k = 0; k < 2; k++){
            fc = function_value_type::Zero();

            fc(k,0) = vx;
            fc(k,1) = vy;

            ret.push_back(fc);
         }
      }

      return ret;
   }

};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class Raviart_Thomas_matrix_basis<Mesh<T,2,Storage>, typename Mesh<T,2,Storage>::face>
: public priv::monomial_basis_bones<1,2>
{
   typedef Mesh<T,2,Storage>                       mesh_type;
   typedef typename mesh_type::face                face_type;
   typedef priv::monomial_basis_bones<1,2>           base;
   typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    scaled_monomial_scalar_face_basis_type;

   scaled_monomial_scalar_face_basis_type smsfb;

public:
   typedef static_matrix<T,2,2>                      function_value_type;

   Raviart_Thomas_matrix_basis()
   : base(1)
   {
      smsfb = scaled_monomial_scalar_face_basis_type(this->degree());
   }

   Raviart_Thomas_matrix_basis(size_t degree)
   : base(degree)
   {
      smsfb = scaled_monomial_scalar_face_basis_type(this->degree());
   }


   //overloaded function

   size_t size() const {
      return 2*(this->degree()+1) * (this->degree()+3);
   }

   size_t computed_size() const {
      return 2*(this->computed_degree()+1) * (this->computed_degree()+3);
   }

   dof_range range(const size_t min_degree, const size_t max_degree) const
   {
      const dof_range vector_range = base::range(min_degree, max_degree);

      if(max_degree == computed_degree())
         return dof_range(2*vector_range.min(), 2*(vector_range.max() + (max_degree + 1)));
      else
         return dof_range(2*vector_range.min(), 2*vector_range.max());
   }

   dof_range range() const {
      return dof_range(0, size());
   }

   dof_range computed_range() const {
      return dof_range(0, computed_size());
   }

   //evaluation

   std::vector<function_value_type>
   eval_functions(const mesh_type& msh, const face_type& fc, const point<T,2>& pt) const
   {
      const auto pts = points(msh, fc);
      const auto bar = barycenter(msh, fc);
      const auto h = diameter(msh, fc);
      const auto v = (pts[1] - pts[0]).to_vector();
      const auto t = (pt - bar).to_vector();
      const T dot = v.dot(t);
      const auto ep = point<T, 1>({dot/(h*h)});

      std::vector<function_value_type> ret;
      ret.reserve(2*(this->degree()+1) * (this->degree()+3));

      // basis functions of P^k_d(T; R^d)
      for (auto itor = this->monomials_begin(); itor != this->monomials_end(); itor++)
      {
         const auto m = *itor;
         const auto vx = iexp_pow(ep.x(), m[0]);

         const auto val = vx ;

         function_value_type fc;

         for( size_t j = 0; j < 2; j++){
            for( size_t i = 0; i < 2; i++){
               fc = function_value_type::Zero();
               fc(i,j) = val;
               ret.push_back(fc);
            }
         }
      }

      // basis functions of P^k,H_d(T; R)
      const auto smsfb_phi = smsfb.eval_functions(msh, fc, pt, this->degree(), this->degree());
      const size_t dim_PHk = this->degree() + 1;
      assert( smsfb_phi.rows() == dim_PHk);

      // basis functions of x*P^k,H_d(T; R)
      for(size_t i = 0; i < dim_PHk; i++){
         const auto vx = pt.x() * smsfb_phi(i);
         const auto vy = pt.y() * smsfb_phi(i);

         function_value_type fc;

         for(size_t k = 0; k < 2; k++){
            fc = function_value_type::Zero();

            fc(k,0) = vx;
            fc(k,1) = vy;

            ret.push_back(fc);
         }
      }

      return ret;
   }
};

} // namespace disk


#endif /* _BASES_ALL_HPP_ */
