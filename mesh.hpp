#pragma once

#include <array>
#include <vector>
#include <string>
#include <utility>
#include <tuple>
#include <iostream>
#include <fstream>
#include <exception>
#include <algorithm>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Core>
#include <boost/format.hpp>
// #include <boost/math/constants/constants.hpp>
#include "structure.hpp"
#include "util.hpp"

namespace cnthd
{

struct Mesh
{
    std::vector<Vertex> vertex;
    std::vector<Edge> edge;
    std::vector<HalfEdge> halfedge;
    std::vector<Face> face;
    std::size_t nV;
    std::size_t nE;
    std::size_t nF;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Eigen::SparseMatrix<bool, Eigen::RowMajor> adj_mat;
    std::vector<std::vector<std::size_t>> adj_list;

    Mesh(const std::vector<std::array<real_t, 3>>& raw_vertices,
         const std::vector<std::array<std::size_t, 3>>& raw_faces);

    int num_genus() const
    {
        return 1 - (nV - nE + nF)/2;
    }

    Mesh& move(const Vector& vec)
    {
        for (Vertex& v: vertex)
        {
            v.p[0] = v.p[0] + vec.x;
            v.p[1] = v.p[1] + vec.y;
            v.p[2] = v.p[2] + vec.z;
        }
        return *this;
    }

    Mesh& rotate(const std::array<real_t, 4>& quaternion)
    {
        const auto& q = quaternion;
        for (Vertex& v: vertex)
        {
            const real_t x = v.x();
            const real_t y = v.y();
            const real_t z = v.z();
            // the formula is calculated using sympy
            v.p[0] = q[0]*q[0]*x + 2*q[0]*q[2]*z - 2*q[0]*q[3]*y + q[1]*q[1]*x + 2*q[1]*q[2]*y + 2*q[1]*q[3]*z - q[2]*q[2]*x - q[3]*q[3]*x;
            v.p[1] = q[0]*q[0]*y - 2*q[0]*q[1]*z + 2*q[0]*q[3]*x - q[1]*q[1]*y + 2*q[1]*q[2]*x + q[2]*q[2]*y + 2*q[2]*q[3]*z - q[3]*q[3]*y;
            v.p[2] = q[0]*q[0]*z + 2*q[0]*q[1]*y - 2*q[0]*q[2]*x - q[1]*q[1]*z + 2*q[1]*q[3]*x - q[2]*q[2]*z + 2*q[2]*q[3]*y + q[3]*q[3]*z;
        }
        return *this;
    }

    Mesh& rotate(const Eigen::Matrix3d& Mat)
    {
        for (Vertex& v: vertex)
        {
            Eigen::Vector3d tmp{v.x(), v.y(), v.z()};
            tmp = Mat * tmp;
            v.p[0] = tmp[0];
            v.p[1] = tmp[1];
            v.p[2] = tmp[2];
        }
        return *this;
    }

    // not implemented
    Mesh& fix_orientation();
};

Mesh::Mesh(const std::vector<std::array<real_t, 3>>& raw_vertices,
           const std::vector<std::array<std::size_t, 3>>& raw_faces)
:nV{raw_vertices.size()}, nF{raw_faces.size()}
{
    //-------------------------------------------------------------------------
    nE = nF * 3 / 2; // triangular mesh
    vertex.reserve(nV);
    edge.reserve(nE);
    halfedge.reserve(nE * 2);
    face.reserve(nF);

    //-------------------------------------------------------------------------
    // build adjacent matrix and adjacent list
    //-------------------------------------------------------------------------
    adj_mat.resize(nV, nV);
    if (nV != 0)
    {
        adj_mat.reserve(Eigen::VectorXi::Constant(nV, 2 * nE / nV + 5));
    }
    adj_list.resize(nV);

    //-------------------------------------------------------------------------
    // build vertex
    //-------------------------------------------------------------------------
    for (std::size_t i = 0; i < nV; ++i)
    {
        const auto& rv = raw_vertices[i];
        vertex.emplace_back(i, rv[0], rv[1], rv[2], /* <HalfEdge*> */ nullptr);
    }

    //-------------------------------------------------------------------------
    // build halfedge structure
    //-------------------------------------------------------------------------
    // Elem = std::size_t
    // unorderedmap = std::pair<Key, Elem>
    // Key is a pair of index for vertex[index]
    // order of vertices is used to check if the other pair of edges is properly stored
    using Key = std::pair<std::size_t, std::size_t>;
    // store already built edge
    unordered_map existedge(nE);

    // loop 1 -----------------------------------------------------------------
    for (std::size_t i = 0; i < nF; ++i)
    {
        const auto& rf = raw_faces[i];
        face.emplace_back(i, /* <HalfEdge*> */ nullptr);

        for (std::size_t j = 0; j < 3; ++j)
        {
            halfedge.emplace_back(
                i*3 + j,
                /* <Vertex*>from= */ &vertex[rf[j%3]], /*<Vertex*>to= */ &vertex[rf[(j+1)%3]],
                /* <HalfEdge*>prev= */ nullptr, /* <HalfEdge*>next= */ nullptr, /* <HalfEdge*>pair= */ nullptr,
                /* <Edge*>edge= */ nullptr, &face[i]);

            Key key{};
            if (rf[j%3] > rf[(j+1)%3])
            {
                key = {rf[j%3], rf[(j+1)%3]};
            }
            else
            {
                key = {rf[(j+1)%3], rf[j%3]};
            }

            const auto iter = existedge.find(key);
            if (iter != existedge.end())
            {
                halfedge[i*3 + j].pair = edge[iter->second].he;
                halfedge[i*3 + j].edge = &edge[iter->second];
            }
            else
            {
                edge.emplace_back(edge.size(), &halfedge[i*3 + j]);
                existedge.emplace(key, edge.size() - 1);
                adj_mat.insert(rf[j%3], rf[(j+1)%3]) = 1;
                adj_mat.insert(rf[(j+1)%3], rf[j%3]) = 1;
                adj_list[rf[j%3]].push_back(rf[(j+1)%3]);
                adj_list[rf[(j+1)%3]].push_back(rf[j%3]);
            }
        }
        face[i].he = &halfedge[i*3];
    }

    // loop 2 -----------------------------------------------------------------
    for (std::size_t i = 0; i < nF; ++i)
    {
        halfedge[i*3].prev = &halfedge[i*3 + 2];
        halfedge[i*3].next = &halfedge[i*3 + 1];
        halfedge[i*3 + 1].prev = &halfedge[i*3];
        halfedge[i*3 + 1].next = &halfedge[i*3 + 2];
        halfedge[i*3 + 2].prev = &halfedge[i*3 + 1];
        halfedge[i*3 + 2].next = &halfedge[i*3];
    }
}

// not implemented
Mesh& Mesh::fix_orientation()
{
    // draft version is below
    // https://github.com/tak0kada/cnthd/commit/e83fdae303cf77ccdc51c4271f7f5d2da1160539#diff-dab052ca81d83cdfb81dc79cf4dcc38d
    // note the algorithm has bugs in it.
    std::cerr << "NOTE: Mesh::fix_orientation is not implmented!" << std::endl;
    return *this;
}

Mesh read_obj(const std::string& path)
{
    std::vector<std::array<real_t, 3>> raw_vertices;
    std::vector<std::array<std::size_t, 3>> raw_faces;

    std::ifstream ifs{path};
    if (ifs.fail())
    {
        throw std::runtime_error("ERROR: cannot open obj file: " + path + ".");
    }

    std::string buf;
    while (ifs >> buf)
    {
        if (buf == "v")
        {
            real_t x, y, z;
            ifs >> x >> y >> z;
            raw_vertices.push_back(std::array<real_t, 3>{x, y, z});
        }
        else if (buf == "f")
        {
            std::size_t v[3] = {};
            for (int i = 0; i < 3; ++i)
            {
                std::string s;
                ifs >> s;
                std::size_t pos{s.find_first_of("/")};
                if (pos == std::string::npos)
                {
                    pos = s.size();
                }
                std::sscanf(s.substr(0, pos).c_str(), "%zu", &v[i]);
            }
            assert (v[0] != 0 && v[1] != 0 && v[2] != 0); // each face has three vetices
            raw_faces.push_back(std::array<std::size_t, 3>{--v[0], --v[1], --v[2]});
        }
        // not implemeted
        else if (buf == "#") {}
        else if (buf == "vt") {}
        else if (buf == "vn") {}
        else if (buf == "o") {}
        else if (buf == "g") {}
        else if (buf == "s") {}
        else if (buf == "mtllib") {}
        else if (buf == "usemtl") {}
        else {}

        // discard the rest of the line
        ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    Mesh mesh{raw_vertices, raw_faces};
    // mesh.fix_orientation();
    return mesh;
}

void write_obj(const Mesh& m, const std::string& path)
{
    std::ofstream ofs;
    try
    {
        ofs.open(path);
    }
    catch (std::ios_base::failure& e)
    {
        std::cerr << e.what() << std::endl;
    }

    for (std::size_t i = 0; i < m.vertex.size(); ++i)
    {
        ofs << "v " << m.vertex[i].x()
            << " "  << m.vertex[i].y()
            << " "  << m.vertex[i].z() << "\n";
    }
    for (const auto& f: m.face)
    {
        const auto v{f.vertex()};
        ofs << "f " << std::get<0>(v)->idx + 1
            << " "  << std::get<1>(v)->idx + 1
            << " "  << std::get<2>(v)->idx + 1 << "\n";
    }
    ofs.close();
}

std::ostream& operator<<(std::ostream& os, const Mesh& mesh)
{
    os << boost::format("nV: %zu, nE: %zu, nF: %zu, g = 1 - 0.5 * (nV - nE + nF) = %d")
          % mesh.nV % mesh.nE % mesh.nF % mesh.num_genus()
       << "\n";

    os << "vertex:\n";
    for (const auto& v: mesh.vertex)
    {
        os << "  " << v << "\n";
    }

    os << "edge:\n";
    for (const auto& e: mesh.edge)
    {
        os << "  " <<  e << "\n";
    }

    os << "halfedge:\n";
    for (const auto&he: mesh.halfedge)
    {
        os << "  "  << he << "\n";
    }

    os << "face:\n";
    for (const auto&f: mesh.face)
    {
        os << "  "  << f << "\n";
    }

    os << "adjacent matrix:\n"
       << static_cast<Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(mesh.adj_mat)
       << "\n"
       << "adjacent list:\n"
       << mesh.adj_list
       << std::endl;

    return os;
}

} // end of namespace cnthd
