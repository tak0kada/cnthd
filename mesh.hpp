#pragma once

#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
// #include <filesystem>
// #include <any>
#include <boost/math/constants/constants.hpp>
#include <eigen3/Eigen/Sparse>
#include "structure.hpp"
#include "util.hpp"

namespace cnthd
{

struct Mesh {
    std::vector<Vertex> vertex;
    std::vector<Edge> edge;
    std::vector<HalfEdge> halfedge;
    std::vector<Face> face;
    std::size_t nV;
    std::size_t nE;
    std::size_t nF;
    Eigen::SparseMatrix<bool, Eigen::RowMajor> adj_mat;
    std::vector<std::vector<std::size_t>> adj_list;

    Mesh(const std::vector<std::array<real_t, 3>>& raw_vertices,
         const std::vector<std::array<std::size_t, 3>>& raw_faces);

    int num_genus() {
        return nV - nE + nF;
    }
};

Mesh::Mesh(const std::vector<std::array<real_t, 3>>& raw_vertices,
           const std::vector<std::array<std::size_t, 3>>& raw_faces)
:nV{raw_vertices.size()}, nF{raw_faces.size()}
{
    //-------------------------------------------------------------------------
    // build adjacent matrix and adjacent list
    //-------------------------------------------------------------------------
    adj_mat.resize(nV, nV);
    adj_list.resize(nV);
    for (const auto& rf: raw_faces) {
        adj_mat.insert(rf[0], rf[1]) = 1; adj_mat.insert(rf[1], rf[0]) = 1;
        adj_mat.insert(rf[0], rf[2]) = 1; adj_mat.insert(rf[2], rf[0]) = 1;
        adj_mat.insert(rf[1], rf[2]) = 1; adj_mat.insert(rf[2], rf[1]) = 1;

        adj_list[rf[0]].push_back(rf[1]); adj_list[rf[1]].push_back(rf[0]);
        adj_list[rf[0]].push_back(rf[2]); adj_list[rf[2]].push_back(rf[0]);
        adj_list[rf[1]].push_back(rf[2]); adj_list[rf[2]].push_back(rf[1]);
    }

    //-------------------------------------------------------------------------
    nE = nF * 3 / 2; // triangular mesh
    vertex.reserve(nV);
    edge.reserve(nE);
    halfedge.reserve(nE * 2);
    face.reserve(nF);

    //-------------------------------------------------------------------------
    // build vertex
    //-------------------------------------------------------------------------
    for (std::size_t i = 0; i < nV; ++i) {
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
    unordered_map hemap0; hemap0.reserve(nE);
    // we use two maps because if the file is corrupted there are inside-out oriented faces
    unordered_map hemap1; hemap1.reserve(nE);

    // loop 1 -----------------------------------------------------------------
    for (std::size_t i = 0; i < nF; ++i) {
        const auto& rf = raw_faces[i];
        face.emplace_back(i, /* <HalfEdge*> */ nullptr);

        for (std::size_t j = 0; j < 3; ++j) {
            halfedge.emplace_back(
                i*3 + j,
                /* <Vertex*>from= */ &vertex[rf[j%3]], /*<Vertex*>to= */ &vertex[(j+1)%3],
                /* <Edge*>prev= */ nullptr, /* <Edge*>next= */ nullptr, /* <Edge*>pair= */ nullptr,
                &face[i]);
            vertex[rf[j%3]].he_out = &halfedge[i*3 + j];

            // keys below provide access to the opposite side (= pair) halfedge
            // key to properly oriented halfedge pair
            Key key0{rf[(j+1)%3], rf[j%3]};
            // key to ill-oriented halfedge pair (or current halfedge itself)
            Key key1{rf[j%3], rf[(j+1)%3]};
            bool exist00{has_key(hemap0, key0)};
            bool exist01{has_key(hemap0, key1)};
            bool exist10{has_key(hemap1, key0)};
            bool exist11{has_key(hemap1, key1)};

            if (exist00) {
                hemap0[key1] = i*3 + j;
                halfedge[i*3 + j].pair = edge[hemap0[key0]].he;
                edge[hemap0[key0]].he->pair = &halfedge[i*3 + j];
            }
            else if (exist01) {
                hemap1[key1] = i*3 + j;
                halfedge[i*3 + j].pair = edge[hemap0[key1]].he;
                edge[hemap0[key1]].he->pair = &halfedge[i*3 + j];
            }
            else if (exist10) {
                hemap1[key1] = i*3 + j;
                halfedge[i*3 + j].pair = edge[hemap1[key0]].he;
                edge[hemap1[key0]].he->pair = &halfedge[i*3 + j];
            }
            else if (exist11) {
                hemap0[key1] = i*3 + j;
                halfedge[i*3 + j].pair = edge[hemap1[key1]].he;
                edge[hemap1[key1]].he->pair = &halfedge[i*3 + j];
            }
            else { // no opposite halfedge is built yet
                edge.emplace_back(edge.size(), &halfedge[i*3 + j]);
            }
        }
        face[i].he = &halfedge[i*3];
    }

    // loop 2 -----------------------------------------------------------------
    for (std::size_t i = 0; i < nF; ++i) {
        const auto& rf = raw_faces[i];
        halfedge[i*3].prev = &halfedge[i*3 + 2];
        halfedge[i*3].next = &halfedge[i*3 + 1];
        halfedge[i*3 + 1].prev = &halfedge[i*3];
        halfedge[i*3 + 1].next = &halfedge[i*3 + 2];
        halfedge[i*3 + 2].prev = &halfedge[i*3 + 1];
        halfedge[i*3 + 2].next = &halfedge[i*3];
    }

    //-------------------------------------------------------------------------
    // check whether mesh is inside-out or properly oriented
    //-------------------------------------------------------------------------
    // an iterator to point to one of the vertex which has the largest x coordinate
    const auto it = std::max_element(vertex.cbegin(), vertex.cend(),
            [](const Vertex& v0, const Vertex& v1) -> bool { return v0.x() > v1.x(); });

    // collect the halfedges that rotate around the face surrounding the vertex (*it)
    std::vector<HalfEdge*> rot_he;
    HalfEdge* he_vert{it->he_out};
    HalfEdge* he_rot{it->he_out->next};
    do {
        rot_he.push_back(he_rot);

        HalfEdge* tmp = he_vert;
        // temporal substitution
        he_vert = tmp->pair->next;
        he_rot = tmp->pair->prev;

        if (he_vert->from->idx != it->idx || he_vert->to->idx != it->idx)
        {
            std::swap(he_vert, he_rot);
        }
    } while (he_vert != it->he_out);

    // generate normalized and aligned vectors correspond to the orientation of hemap0
    std::vector<Vector> rot_vec;
    for (const auto& he: rot_he) {
        if (has_key(hemap0, {he->from->idx, he->to->idx})) {
            Vector v{*he->to - *he->from};
            const real_t l = v.length();
            v.x /= l; v.y /= l; v.z /= l;
            rot_vec.push_back(v);
        }
        else {
            Vector v{*he->from - *he->to};
            const real_t l = v.length();
            v.x /= l; v.y /= l; v.z /= l;
            rot_vec.push_back(v);
        }
    }

    std::size_t N{rot_vec.size()};
    Vector rot{0, 0, 0};
    for (std::size_t i = 0 ; i < N; ++i) {
        rot = rot + (rot_vec[i] ^ rot_vec[(i+1)%N]);
    }

    // inversion --------------------------------------------------------------
    // if rot.x > 0, hemap0 is the correct orientation
    // if rot.x < 0, hemap1 is the correct orientation (inside-out)
    // rot.x == 0 can not happen
    bool insideout;
    if (rot.x > 0) insideout = false;
    else insideout = true;

    if (insideout /* && hemap0.size() > 0 */) {
        for (const std::pair<Key, std::size_t>& p : hemap0) {
            HalfEdge& he{halfedge[p.second]};
            Vertex * const v_tmp{he.to};
            he.to = he.from;
            he.from = v_tmp;
            he.from->he_out = &he;
            HalfEdge * const h_tmp{he.next};
            he.next = he.prev;
            he.prev = h_tmp;
        }
    }
    else if (!insideout && hemap1.size() > 0) {
        for (const auto& p : hemap1) {
            HalfEdge& he{halfedge[p.second]};
            Vertex * const v_tmp{he.to};
            he.to = he.from;
            he.from = v_tmp;
            he.from->he_out = &he;
            HalfEdge * const h_tmp{he.next};
            he.next = he.prev;
            he.prev = h_tmp;
        }
    }
}

Mesh read_obj(const std::string& path)
{
    std::vector<std::array<real_t, 3>> raw_vertices;
    std::vector<std::array<std::size_t, 3>> raw_faces;

    std::ifstream ifs{path};
    std::string buf;

    while (ifs.good()) {
        ifs >> buf;
        if (buf == "v") {
            real_t x, y, z;
            ifs >> x >> y >> z;
            raw_vertices.emplace_back(x, y, z);
        }
        else if (buf == "f") {
            std::size_t v0{}, v1{}, v2{};
            ifs >> v0 >> v1 >> v2;
            assert (v0 != 0 && v1 != 0 && v2 != 0); // each face has three vetices
            raw_faces.emplace_back(--v0, --v1, --v2);
        }
        else if (buf == "#") {
            continue;
        }
        else {
            // not implemeted
        }
        // discard the rest of the line
        ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    return Mesh{raw_vertices, raw_faces};
}

void write_obj(const Mesh& m, const std::string& path)
{
    std::ofstream ofs;
    try {
        ofs.open(path);
    }
    catch (std::ios_base::failure& e){
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
        const std::array<Vertex*, 3> v{f.vertex()};
        ofs << "f " << v[0]->idx + 1
            << " "  << v[1]->idx + 1
            << " "  << v[2]->idx + 1 << "\n";
    }
    ofs.close();
}

} // end of namespace cnthd
