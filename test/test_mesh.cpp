#include <iostream>
#include <string>
#include "mesh.hpp"

int main(int argc, char* argv[])
{
    using namespace cnthd;

    if (argc == 1)
    {
        std::cout << "./a.out input.obj output.obj" << std::endl;
        return 0;
    }

    std::string infile{argv[1]};
    std::string outfile{argv[2]};

    Mesh mesh = read_obj(infile);
    std::cout << mesh << std::endl;
    write_obj(mesh, outfile);

    return 0;
}
