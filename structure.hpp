#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <utility>
#include <tuple>
#include <iomanip>
#include "util.hpp"
// #include "vector.hpp"


namespace cnthd {

//-----------------------------------------------------------------------------
// foward declaration
//-----------------------------------------------------------------------------

struct HalfEdge;
struct Vertex;
struct Edge;
struct HalfEdge;
struct Face;


//-----------------------------------------------------------------------------
// definition
//-----------------------------------------------------------------------------

struct HalfEdge
{
    std::size_t idx;
    Vertex* from;
    Vertex* to;
    HalfEdge* prev;
    HalfEdge* next;
    HalfEdge* pair;
    Edge* edge;
    Face* face;

    HalfEdge(
        const std::size_t& idx,
        Vertex * const from, Vertex * const to,
        HalfEdge * const prev, HalfEdge * const next, HalfEdge * const pair,
        Edge * const edge, Face * const face)
    : idx{idx}, from{from}, to{to}, prev{prev}, next{next}, pair{pair}, edge{edge}, face{face}
    {}
};

struct Vertex
{
    std::size_t idx;
    std::array<real_t, 3> p;
    HalfEdge* he_out;

    Vertex(std::size_t idx, real_t x, real_t y, real_t z, HalfEdge* he_out)
    :idx{idx}, p{x, y, z}, he_out{he_out}
    {}

    inline real_t& x() { return p[0]; }
    inline const real_t& x() const { return p[0]; }
    inline real_t& y() { return p[1]; }
    inline const real_t& y() const { return p[1]; }
    inline real_t& z() { return p[2]; }
    inline const real_t& z() const { return p[2]; }
};

struct Vector
{
    real_t x;
    real_t y;
    real_t z;

    Vector(const real_t& x, const real_t& y, const real_t& z)
    : x{x}, y{y}, z{z}
    {}

    Vector(const Vector& vec)
    : x{vec.x}, y{vec.y}, z{vec.z}
    {}

    inline real_t length() const
    {
        return std::sqrt(x * x + y * y + z * z);
    }

    inline real_t norm() const
    {
        return x * x + y * y + z * z;
    }

    Vector operator-() const
    {
        return {-x, -y, -z};
    }

    Vector& operator+=(const Vector& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    Vector& operator-=(const Vector& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    Vector operator*(const real_t& r) const
    {
        return {x * r, y * r, z * r};
    }

    inline real_t operator* (const Vector& v) const
    {
        return x * v.x + y * v.y + z * v.z;
    }
};

Vector unit_vec(const Vector& vec)
{
    real_t l = vec.length();
    return {vec.x / l, vec.y / l, vec.z / l};
}

inline Vector operator- (const Vertex& lhs, const Vertex& rhs)
{
    return {lhs.x() - rhs.x(), lhs.y() - rhs.y(), lhs.z() - rhs.z()};
}

inline Vector operator+ (const Vector& lhs, const Vector& rhs)
{
    return Vector{lhs} += rhs;
}

inline Vector operator- (const Vector& lhs, const Vector& rhs)
{
    return Vector{lhs} -= rhs;
}

inline Vector operator^ (const Vector& lhs, const Vector& rhs)
{
    return {lhs.y * rhs.z - lhs.z * rhs.y,
            lhs.z * rhs.x - lhs.x * rhs.z,
            lhs.x * rhs.y - lhs.y * rhs.x};
}

bool operator==(const Vector& lhs, const Vector& rhs)
{
    return (lhs.x == rhs.x || lhs.y == rhs.y || lhs.z == rhs.z);
}


struct Edge
{
    std::size_t idx;
    HalfEdge* he;

    Edge(const std::size_t& idx, HalfEdge * const he)
    : idx{idx}, he{he}
    {}

    real_t norm() const
    {
        const Vector v{*he->to - *he->from};
        return v.x * v.x + v.y * v.y + v.z * v.z;
    }

    real_t length() const
    {
        return std::sqrt(norm());
    }

    std::pair<Vertex*, Vertex*> vertex()
    {
        return {he->from, he->to};
    }
    const std::pair<Vertex*, Vertex*> vertex() const
    {
        return {he->from, he->to};
    }
};

struct Face
{
    std::size_t idx;
    HalfEdge* he;

    Face(const std::size_t& idx, HalfEdge * const he)
    :idx{idx}, he{he}
    {}

    std::tuple<Vertex*, Vertex*, Vertex*> vertex()
    {
        return {he->from, he->to, he->next->to};
    }
    const std::tuple<Vertex* const, Vertex* const, Vertex* const> vertex() const
    {
        return {he->from, he->to, he->next->to};
    }

    std::tuple<Edge*, Edge*, Edge*> edge()
    {
        return {he->edge, he->next->edge, he->prev->edge};
    }
    const std::tuple<Edge* const, Edge* const, Edge* const> edge() const
    {
        return {he->edge, he->next->edge, he->prev->edge};
    }

    real_t area() const
    {
        const Vector v0{*he->to - *he->from};
        const Vector v1{*he->next->to - *he->to};
        const Vector v2{*he->from - *he->prev->from};
        const real_t a{v0.norm()};
        const real_t b{v1.norm()};
        const real_t c{v2.norm()};
        const real_t d{2*a*b + 2*b*c + 2*c*a - a*a - b*b - c*c};

        if (d <= 0) return 0;
        return std::sqrt(d * 0.0625);
    }

    Vector normal() const
    {
        const Vector v0{*he->to - *he->from};
        const Vector v1{*he->prev->from - *he->prev->to};
        const Vector v2{v0^v1};
        return {v2.x / v2.length(), v2.y / v2.length(), v2.z / v2.length()};
    }
};

std::ostream& operator<<(std::ostream& os, const HalfEdge& he)
{
    os << "he" << he.idx + 1 << ": " << he.from->idx + 1 << "->" << he.to->idx + 1;
    return os;
}

std::ostream& operator<<(std::ostream& os, const Vertex& v)
{
    os << "v" << v.idx + 1 << ": {"
       << std::fixed << std::setprecision(6)
       <<  v.x() << ", " << v.y() << ", " << v.z()
       << std::defaultfloat << "}";
    return os;
}

std::ostream& operator<<(std::ostream& os, const Vector& v)
{
    os << "vec " << "{" << v.x << ", " << v.y << ", " << v.z << "}";
    return os;
}

std::ostream& operator<<(std::ostream& os, const Edge& e)
{
    os << "e" << e.idx + 1 << ": " << e.he->from->idx + 1 << "-" << e.he->to->idx + 1;
    return os;
}

std::ostream& operator<<(std::ostream& os, const Face& f)
{
    os << "f" << f.idx + 1 << ": " << f.he->from->idx + 1 << "-" << f.he->to->idx + 1 << "-" << f.he->next->to->idx + 1;
    return os;
}

} // end of namespace cnthd
