#pragma once

#include <vector>
#include <array>
#include <cmath>
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
struct Mesh;


//-----------------------------------------------------------------------------
// definition
//-----------------------------------------------------------------------------

struct HalfEdge {
public:
    std::size_t idx;
    Vertex* from;
    Vertex* to;
    HalfEdge* prev;
    HalfEdge* pair;
    HalfEdge* next;
    Edge* edge;
    Face* face;

private:
    friend Mesh;
    explicit HalfEdge(
        const std::size_t& idx, Edge * const edge,
        HalfEdge * const prev, HalfEdge * const next, HalfEdge * const pair,
        Face * const face)
    : idx{idx}, edge{edge}, prev{prev}, next{next}, pair{pair}, face{face}
    {}
};

struct Vertex {
public:
    std::size_t idx;
    std::array<real_t, 3> p;
    HalfEdge* he_out;

    real_t x() { return p[0]; }
    const real_t x() const { return p[0]; }
    real_t y() { return p[1]; }
    const real_t y() const { return p[1]; }
    real_t z() { return p[2]; }
    const real_t z() const { return p[2]; }

private:
    friend Mesh;
    Vertex(std::size_t idx, real_t x, real_t y, real_t z, HalfEdge* he_out)
    :idx{idx}, p{x, y, z}, he_out{he_out}
    {}
};

struct Vector {
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

    real_t x;
    real_t y;
    real_t z;
};

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

inline real_t operator* (const Vector& lhs, const Vector& rhs)
{
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

inline Vector operator^ (const Vector& lhs, const Vector& rhs)
{
    return {lhs.y * rhs.z - lhs.z * rhs.y,
            lhs.z * rhs.x - lhs.x * rhs.z,
            lhs.x * rhs.y - lhs.y * rhs.x};
}

struct Edge {
public:
    std::size_t idx;
    HalfEdge* he;

    real_t norm() const
    {
        const Vector v{*he->to - *he->from};
        return v.x * v.x + v.y * v.y + v.z * v.z;
    }

    real_t length() const
    {
        return std::sqrt(norm());
    }

    std::array<Vertex*, 2> vertex()
    {
        return {he->from, he->to};
    }
    const std::array<Vertex*, 2> vertex() const
    {
        return {he->from, he->to};
    }

private:
    friend Mesh;
    explicit Edge(const std::size_t& idx, HalfEdge * const he)
    : idx{idx}, he{he}
    {}
};

struct Face {
public:
    std::size_t idx;
    HalfEdge* he;

    std::array<Vertex*, 3> vertex()
    {
        return {he->from, he->to, he->next->to};
    }
    const std::array<Vertex*, 3> vertex() const
    {
        return {he->from, he->to, he->next->to};
    }
    std::array<Edge*, 3> edge()
    {
        return {he->edge, he->next->edge, he->prev->edge};
    }
    const std::array<Edge*, 3> edge() const
    {
        return {he->edge, he->next->edge, he->prev->edge};
    }

    real_t area() const
    {
        const Vector v0{*he->from - *he->to};
        const Vector v1{*he->from - *he->to};
        const Vector v2{*he->from - *he->to};
        const real_t a{v0.norm()};
        const real_t b{v1.norm()};
        const real_t c{v2.norm()};
        const real_t d{2*a*b + 2*b*c + 2*c*a - a*a - b*b - c*c};

        if (d <= 0) return 0;
        return std::sqrt(d) * 0.625;
    }

    Vector normal() const
    {
        const Vector v0{*he->to - *he->from};
        const Vector v1{*he->prev->from - *he->prev->to};
        const Vector v2{v0^v1};
        return {v2.x / v2.length(), v2.y / v2.length(), v2.z / v2.length()};
    }

private:
    friend Mesh;
    explicit Face(const std::size_t& idx, HalfEdge * const he)
    :idx{idx}, he{he}
    {}
};

} // end of namespace cnthd
