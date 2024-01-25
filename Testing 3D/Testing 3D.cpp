// Testing 3D.cpp : This file contains the 'main' function. Program execution begins and ends there.
//


#include <iostream>
#include <filesystem>
#include "cage_deform_gc_2D_3D.h"
#include <fstream>

using namespace std;
using namespace GreenCoordinates;

int main(int argc, char* argv[])
{
    if (argc != 3) {
        cerr << "usage: " << argv[0] << " model.obj cage.obj\n";
        return __LINE__;
    }


    // Create an instance of the "green_deform_3d" class
    green_coords_3d def;

    // Load the sample points from the provided file
    def.load_model_points(argv[1]);

    // Load the cage from the provided file
    def.load_cage(argv[2]);

    // Calculate the Green Coordinates for the loaded sample points and cage
    def.calculate_green_coords();

    // Move specific vertices of the cage to specified positions (bar example)
    const double pos4[3] = { -1.031485, 4.966156, 1.048004 };
    const double pos5[3] = { -1.031485, 4.966156, -1.046975 };
    const double pos6[3] = { 1.063854, 4.966156, -1.046975 };
    const double pos7[3] = { 1.063854, 4.966156, 1.048004 };
    def.move_cage(4, pos5, false);
    def.move_cage(5, pos6, false);
    def.move_cage(6, pos7, false);
    def.move_cage(7, pos4, false);

    // Perform the deformation
    def.deform();

    // Save the deformed mesh to a obj file
    def.saveToFile("./green3d/mesh.obj");

    // Save the modified cage to a VTK file
    def.saveToFile_cage("./green3d/cage.obj");

    // Save the normals of the modified cage to a VTK file
    def.saveToFile_normal("./green3d/cage_normal.obj");

    cout << "done\n";

    return 0;
}
