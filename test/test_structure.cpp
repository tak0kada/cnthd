#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/detail/tolerance_manip.hpp>
#include <iostream>
#include <algorithm>

#include "structure.hpp"
#include "mesh.hpp"
#include "util.hpp"
#include "unitsp.hpp"

using namespace cnthd;

BOOST_AUTO_TEST_CASE(test_Vector)
{
    Vector x{1, 0, 0};
    Vector y{0, 1, 0};
    Vector z{0, 0, 1};

    BOOST_TEST((x^y) == z);
    BOOST_TEST((y^z) == x);
    BOOST_TEST((z^x) == y);
}

BOOST_AUTO_TEST_CASE(test_Face_area)
{
    Mesh mesh = read_obj("./test/data/cube_small.obj");

    // std::cout << mesh << std::endl;
    for (const auto& f: mesh.face) {
        BOOST_TEST(f.area() == 0.5, boost::test_tools::tolerance(1e-12));
    }
}
