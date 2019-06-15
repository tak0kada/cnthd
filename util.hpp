#pragma once

#include <array>
#include <vector>
#include <unordered_map>
#include <utility>
#include <functional>
#include <iostream>
#include <boost/math/constants/constants.hpp>


namespace cnthd
{

using real_t = double;
constexpr real_t pi = boost::math::constants::pi<real_t>();

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

std::ostream& operator<<(std::ostream& os, const unordered_map& umap)
{
    for (const auto& pair: umap)
    {
        os << "{";
        const auto& key{std::get<0>(pair)};
        const auto& elem{std::get<1>(pair)};

        os << "{" << std::get<0>(key) << "-> " << std::get<1>(key) << "}: " << elem << "\n";
    }
    os << "\b";
    return os;
}


template <typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr)
{
    if (arr.empty())
    {
        os << "{}";
    }
    else{
        os << "{";
        for (const auto& a: arr)
        {
            os << a << ", ";
        }
        os << "\b\b}";
    }
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
    if (vec.empty())
    {
        os << "{}";
    }
    else{
        os << "{";
        for (const auto& v: vec)
        {
            os << v << ", ";
        }
        os << "\b\b}";
    }
    return os;
}

} // end of namespace cnthd
