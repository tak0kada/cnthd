#pragma once

#include <unordered_map>
#include <functional>
#include <boost/math/constants/constants.hpp>


namespace cnthd
{

using real_t = double;
static constexpr real_t pi = boost::math::constants::pi<real_t>();

struct pair_hash
{
    template <typename T0, typename T1>
    std::size_t operator() (const std::pair<T0, T1>& x) const
    {
        return std::hash<T0>{}(x.first) ^ std::hash<T1>{}(x.second);
    }
};

using unordered_map = std::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t, pair_hash>;


bool has_key(const unordered_map& map, const std::pair<std::size_t, std::size_t>& k)
{
    return map.find(k) != map.end();
}

} // end of namespace cnthd
