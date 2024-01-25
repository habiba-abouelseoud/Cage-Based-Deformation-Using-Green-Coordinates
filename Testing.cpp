

#include <iostream>
#include <filesystem>
#include "cage_deform_gc_2D_3D.h"
#include <fstream>

using namespace std;
using namespace GreenCoordinates;

int main(int argc, char* argv[])
{
    if (argc != 3) {
        cerr << "usage: " << argv[0] << " model.obj cage.2d\n";
        return __LINE__;
    }

    // Create an instance of the "green_coords_2d" class
    green_coords_2d def;

    // Load the model points from the provided file
    def.load_model_points(argv[1]);

    // Load the cage from the provided file
    def.load_cage(argv[2]);

    // Save the original cage to a obj file
    def.saveToFile_cage("./green/cage_origin.obj");

    // Calculate the Green Coordinates for the loaded model and cage
    def.calculate_green_coords();

    {
        // Move the third vertex of the cage to a specified position
        const double pos[2] = { -0.1, 1.1 };
        def.move_cage(3, pos, false);
    }
    {
        // Move the fourth vertex of the cage to a specified position
        const double pos[2] = { 0.9, 1.1 };
        def.move_cage(4, pos, false);
    }
    {
        // Move the fifth vertex of the cage to a specified position
        const double pos[2] = { 0.9, -0.1 };
        def.move_cage(5, pos, false);
    }

    // Save the model with the original cage to a obj file
    def.saveToFile("./green/origin.obj");

    // Perform the deformation
    def.deform();

    // Save the deformed mesh to a obj file
    def.saveToFile("./green/mesh.obj");

    // Save the modified cage to a obj file
    def.saveToFile_cage("./green/cage.obj");

    // Save the normals of the modified cage to a obj file
    def.saveToFile_normal("./green/cage_normal.obj");

    cout << "done\n";

    return 0;
}

