
#include <iostream>
#include <filesystem>
#include "cage_deform_gc_2D.h"
#include<fstream>

using namespace std;
using namespace GreenCoordinates;

int main(int argc, char* argv[])
{
    if (argc != 3) {
        cerr << "usage: " << argv[0] << " model.obj cage.2d\n";
        return __LINE__;
    }
   

    green_coords_2d def;
    def.load_model_points(argv[1]);
    def.load_cage(argv[2]);
    def.saveToFile_cage("./green/cage_origin.vtk");
    def.calculate_green_coords();
    {
        const double pos[2] = { -0.1, 1.1 };
        def.move_cage(3, pos, false);
    }
    {
        const double pos[2] = { 0.9, 1.1 };
        def.move_cage(4, pos, false);
    }
    {
        const double pos[2] = { 0.9, -0.1 };
        def.move_cage(5, pos, false);
    }
    def.saveToFile("./green/origin.vtk");
    def.deform();

    def.saveToFile("./green/mesh.vtk");
    def.saveToFile_cage("./green/cage.vtk");
    def.saveToFile_normal("./green/cage_normal.vtk");

    cout << "done\n";
    return 0;
}

