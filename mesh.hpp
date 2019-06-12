#include <array>
#include <vector>
#include <fstream>
#include <string>
#include <unordered_map>
#include <functional>
#include <algorithm>
// #include <filesystem>
// #include <any>
#include <boost/math/constants/constants.hpp>
#include <eigen3/Eigen/Sparse>

namespace std {
    template <>
    struct hash<pair<size_t, size_t>> {
        size_t operator() (const std::pair<size_t, size_t> x) const {
            return hash<size_t>()(x.first) ^ hash<size_t>()(x.second);
        }
    };
}

// taken from scientific name of trigger fish (アミモンガラ); Canthidermis maculata
namespace cnthd {

using real_t = double;
constexpr real_t pi = boost::math::constants::pi<real_t>();

template <typename Key, typename Elem>
bool has_key(const std::unordered_map<Key, Elem>& map, const Key& k) {
    return map.find(k) != map.end();
}


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
struct Vertex {
public:
    std::size_t idx;
    std::array<real_t, 3> p;
    HalfEdge* he;

    real_t x() { return p[0]; }
    const real_t x() const { return p[0]; }
    real_t y() { return p[1]; }
    const real_t y() const { return p[1]; }
    real_t z() { return p[2]; }
    const real_t z() const { return p[2]; }
    std::vector<Vertex&> neighbors();
    std::vector<Edge&> edge();
    std::vector<Face&> face();

private:
    friend Mesh;
    Vertex(std::size_t idx, real_t x, real_t y, real_t z, HalfEdge* he)
    :idx{idx}, p{x, y, z}, he{he}
    {}
};

struct Vector {
    Vector(const real_t& x, const real_t& y, const real_t& z)
    : x{x}, y{y}, z{z}
    {}

    Vector(const Vector& vec)
    : x{vec.x}, y{vec.y}, z{vec.z}
    {}

    inline real_t length() const {
        return std::sqrt(x * x + y * y + z * z);
    }
    inline real_t norm() const {
        return x * x + y * y + z * z;
    }

    inline Vector& operator+=(const Vector& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    inline Vector& operator-=(const Vector& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    real_t x;
    real_t y;
    real_t z;
};

inline Vector operator- (const Vertex& lhs, const Vertex& rhs) {
    return {lhs.x() - rhs.x(), lhs.y() - rhs.y(), lhs.z() - rhs.z()};
}

inline Vector operator+ (const Vector& lhs, const Vector& rhs) {
    return Vector{lhs} += rhs;
}

inline Vector operator- (const Vector& lhs, const Vector& rhs) {
    return Vector{lhs} -= rhs;
}

inline real_t operator* (const Vector& lhs, const Vector& rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

inline Vector operator^ (const Vector& lhs, const Vector& rhs) {
    return {lhs.y * rhs.z - lhs.z * rhs.y,
            lhs.z * rhs.x - lhs.x * rhs.z,
            lhs.x * rhs.y - lhs.y * rhs.x};
}

struct Edge {
    std::size_t idx;
    HalfEdge* he;

    real_t length();
    std::array<Vertex&, 2> vertex();
    void flip();

private:
    friend Mesh;
    explicit Edge(const std::size_t& idx, HalfEdge * const he)
    : idx{idx}, he{he}
    {}
};

struct HalfEdge {
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

struct Face {
    std::size_t idx;
    HalfEdge* he;
    std::array<Vertex&, 3> vertex();
    std::array<Edge&, 3> edge();
    real_t area() const;
    Vector normal() const;
private:
    friend Mesh;
    explicit Face(const std::size_t& idx, HalfEdge * const he)
    :idx{idx}, he{he}
    {}
};

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

    static Mesh read_obj(const std::string& path) {
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

    static void write_obj(const Mesh& m, std::string path) {}

    int num_genus() {
        return nV - nE + nF;
    }

private:
    Mesh() = delete;
    Mesh(std::vector<std::array<real_t, 3>> raw_vertices, std::vector<std::array<std::size_t, 3>> raw_faces)
    :nV{raw_vertices.size()}, nF{raw_faces.size()}
    {
        // build adjacent matrix and adjacent list
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

        // consider only triangular mesh
        nE = nF * 3 / 2;
        vertex.reserve(nV);
        edge.reserve(nE);
        halfedge.reserve(nE * 2);
        face.reserve(nF);

        for (std::size_t i = 0; i < nV; ++i) {
            const auto& rv = raw_vertices[i];
            vertex.emplace_back(i, rv[0], rv[1], rv[2], /* <HalfEdge*> */ nullptr);
        }

        // pair of index for vertex[index] that is a key to edge
        // order of vertices is used to check if the other pair of edges is properly stored
        using Key = std::pair<std::size_t, std::size_t>;
        // make correspondence between edge and index for halfedge[index]
        std::unordered_map<Key, std::size_t> hemap0; hemap0.reserve(nE);
        // we use two maps because if the file is corrupted there are inside-out oriented faces
        std::unordered_map<Key, std::size_t> hemap1; hemap1.reserve(nE);

        for (std::size_t i = 0; i < nF; ++i) {
            const auto& rf = raw_faces[i];
            face.emplace_back(face.size(), /* <HalfEdge*> */ nullptr);

            for (std::size_t j = 0; j < 3; ++j) {
                halfedge.emplace_back(
                    i*3 + j,
                    /* <Vertex*>from= */ &vertex[rf[j%3]], /*<Vertex*>to= */ &vertex[(j+1)%3],
                    /* <Edge*>prev= */ nullptr, /* <Edge*>next= */ nullptr, /* <Edge*>pair= */ nullptr,
                    &face[face.size()]);
                vertex[rf[j%3]].he = &halfedge[halfedge.size() - 1];

                Key key0{rf[(j+1)%3], rf[j%3]};
                Key key1{rf[j%3], rf[(j+1)%3]};
                bool exist00 = has_key(hemap0, key0);
                bool exist01 = has_key(hemap0, key1);
                bool exist10 = has_key(hemap1, key0);
                bool exist11 = has_key(hemap1, key1);

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
            face[face.size() - 1].he = &halfedge[i*3];
        }

        // make orientation of faces consistent
        for (std::size_t i = 0; i < nF; ++i) {
            const auto& rf = raw_faces[i];
            halfedge[i*3].prev = &halfedge[i*3 + 2];
            halfedge[i*3].next = &halfedge[i*3 + 1];
            halfedge[i*3 + 1].prev = &halfedge[i*3];
            halfedge[i*3 + 1].next = &halfedge[i*3 + 2];
            halfedge[i*3 + 2].prev = &halfedge[i*3 + 1];
            halfedge[i*3 + 2].next = &halfedge[i*3];
        }

        // check whether mesh is inside-out or properly oriented --------------
        // it is an iterator to point the vertex which has the largest x coordinate
        const auto it = std::max_element(vertex.cbegin(), vertex.cend(),
                [](const Vertex& v0, const Vertex& v1) -> bool { return v0.x() > v1.x(); });

        // collect the halfedges that rotate around the vertex
        std::vector<HalfEdge*> rot_he;
        HalfEdge* he{it->he};
        do {
            // he is outgoing from the base vertex(*it)
            if (he->prev->to->idx == it->idx) {
                rot_he.push_back(he->next);
                he = he->prev->pair;
            }
            // he is incoming to the base vertex(*it)
            else if (he->next->from->idx == it->idx) {
                rot_he.push_back(he->prev);
                he = he->next->pair;
            }
        } while (he != it->he);

        // collect normalized and aligned vectors correspond to rot_he
        std::vector<Vector> rot_v;
        for (auto& he: rot_he) {
            if (has_key(hemap0, {he->from->idx, he->to->idx})) {
                Vector v{*he->to - *he->from};
                real_t l = v.length();
                v.x /= l; v.y /= l; v.z /= l;
                rot_v.push_back(v);
            }
            else {
                Vector v{*he->from - *he->to};
                real_t l = v.length();
                v.x /= l; v.y /= l; v.z /= l;
                rot_v.push_back(v);
            }
        }

        std::size_t N{rot_v.size()};
        Vector normal{0, 0, 0};
        for (std::size_t i = 0 ; i < N; ++i) {
            normal = normal + (rot_v[i] ^ rot_v[(i+1)%N]);
        }

        bool insideout;
        if (normal.x > 0) insideout = false;
        else insideout = true;

        if (insideout && hemap0.size() > 0) {
            for (const auto& p : hemap0) {
                HalfEdge& he{halfedge[p.second]};
                Vertex * const v_tmp{he.to};
                he.to = he.from;
                he.from = v_tmp;
                he.from->he = &he;
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
                he.from->he = &he;
                HalfEdge * const h_tmp{he.next};
                he.next = he.prev;
                he.prev = h_tmp;
            }
        }
    }
};

}
