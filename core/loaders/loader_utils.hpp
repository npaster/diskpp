/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */

#pragma once

#include <algorithm>
#include <cctype>
#include <locale>
#include <thread>
#include <vector>

namespace disk {

namespace priv {

template<typename T>
void
sort_uniq(std::vector<T>& v)
{
    std::sort(v.begin(), v.end());
    auto uniq_iter = std::unique(v.begin(), v.end());
    v.erase(uniq_iter, v.end());
}

} //namespace priv

// trim from start (in place)
inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
inline void
trim(std::string& s)
{
    rtrim(s);
    ltrim(s);
}

// trim from start (copying)
inline std::string
ltrim_copy(std::string s)
{
    ltrim(s);
    return s;
}

// trim from end (copying)
inline std::string
rtrim_copy(std::string s)
{
    rtrim(s);
    return s;
}

// trim from both ends (copying)
inline std::string
trim_copy(std::string s)
{
    trim(s);
    return s;
}

} //namespace disk

#define THREADED
#ifdef THREADED
    #define THREAD(name, body) std::thread name([&]{body})
    #define WAIT_THREAD(name) name.join()
#else
    #define THREAD(name, body) {body}
    #define WAIT_THREAD(name)
#endif
