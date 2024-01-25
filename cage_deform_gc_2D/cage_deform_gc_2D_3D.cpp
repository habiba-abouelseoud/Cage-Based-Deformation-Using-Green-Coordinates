

#include <iostream>
#include "cage_deform_gc_2D_3D.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <omp.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <cmath>
#include <Eigen/Core>
#include <math.h>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matd_t;


using namespace std;
using namespace Eigen;
#define M_PI 3.14159

namespace GreenCoordinates {

	green_coords_2d::green_coords_2d() {}

	int green_coords_2d::load_model_points(const char* file) {
		std::ifstream infile(file);
		if (!infile) {
			return __LINE__; // Unable to open the file, return the line number as an error code
		}

		// Data containers for nodal and cell information
		std::vector<std::vector<double>> nods_data;  // Nodal coordinates (x, z)
		std::vector<std::vector<size_t>> cell_data;  // Cell connectivity (vertex indices)
		std::string line;

		// Read the file line by line
		while (std::getline(infile, line)) {
			std::istringstream iss(line);
			std::string keyword;
			iss >> keyword;

			if (keyword == "v") {
				// Process a line containing nodal coordinates
				double x, y, z;
				iss >> x >> y >> z;
				nods_data.push_back({ x, z }); // Store x and z coordinates of the node
			}
			else if (keyword == "f") {
				// Process a line containing cell connectivity
				size_t v1, v2, v3;
				iss >> v1 >> v2 >> v3;
				cell_data.push_back({ v1, v2, v3 }); // Store vertex indices of the cell
			}
		}

		infile.close(); // Close the input file

		// Resize data containers based on the loaded data
		cell_.resize(cell_data.size(), 3); // Resize cell_ matrix to the number of cells with 3 vertices each
		nods_.resize(2, nods_data.size()); // Resize nods_ matrix to 2 rows and the number of nodes

		// Populate cell data
		for (size_t i = 0; i < cell_data.size(); ++i) {
			// Store vertex indices of each cell, subtracting 1 to account for zero-based indexing
			cell_(i, 0) = cell_data[i][0] - 1;
			cell_(i, 1) = cell_data[i][1] - 1;
			cell_(i, 2) = cell_data[i][2] - 1;
		}

		// Populate nodal data
		for (size_t i = 0; i < nods_data.size(); ++i) {
			// Store x and z coordinates of each node
			nods_(0, i) = nods_data[i][0];
			nods_(1, i) = nods_data[i][1];
		}

		return 0; // Return 0 to indicate successful loading of model points
	}



	//--------------------------------------------------------
	// testing function for load_model_points 
	//---------------------------------------------------------

	//void test_load_model_points() {
	//    green_coords_2d instance;
	//    const char* file = "model.obj";

	//    int result = instance.load_model_points(file);
	//    if (result != 0) {
	//        std::cout << "Failed to load model points. Error code: " << result << std::endl;
	//        return;
	//    }

	//    // Verify the loaded data
	//    const Matrix<size_t>& cell = instance.get_cell();
	//    const Matrix<double>& nods = instance.get_nods();

	//    std::cout << "Loaded cell matrix:" << std::endl;
	//    std::cout << cell << std::endl;

	//    std::cout << "Loaded nods matrix:" << std::endl;
	//    std::cout << nods << std::endl;

	//    // Perform additional assertions or tests based on the loaded data
	//    // ...

	//    std::cout << "Model points loaded successfully!" << std::endl;
	//}

	//------------------------------------------------------------------------------------

	int green_coords_2d::load_cage(const char* file) {
		// Parse the file manually
		ifstream is(file);
		if (is.fail()) {
			cerr << "# can not open " << file << "\n";
			return __LINE__; // Unable to open the file, return the line number as an error code
		}

		// Read the number of cage elements
		size_t ele_num;
		string CELL;
		is >> CELL >> ele_num;
		cage_cell_.resize(2, ele_num); // Resize cage_cell_ matrix to have 2 rows and ele_num columns

		// Read the cage element vertices
		for (size_t i = 0; i < ele_num; ++i)
			is >> cage_cell_(0, i) >> cage_cell_(1, i); // Store vertex coordinates of each cage element

		// Read the number of cage nodes
		size_t nods_num;
		string NODES;
		is >> NODES >> nods_num;
		cage_nods_.resize(2, nods_num); // Resize cage_nods_ matrix to have 2 rows and nods_num columns

		// Read the cage node coordinates
		for (size_t i = 0; i < nods_num; ++i)
			is >> cage_nods_(0, i) >> cage_nods_(1, i); // Store coordinates of each cage node

		is.close(); // Close the input file

		// Get the rest segment element length and initialize current segment lengths
		rest_len_.resize(cage_cell_.cols()); // Resize rest_len_ vector to have the same number of elements as cage_cell_ columns
		curr_len_.resize(cage_cell_.cols()); // Resize curr_len_ vector to have the same number of elements as cage_cell_ columns
		calculate_cage_edge_length(rest_len_); // Calculate the rest segment element lengths

		// Calculate the initial outward normal vectors of the cage elements
		cage_normal_.resize(2, cage_cell_.cols()); // Resize cage_normal_ matrix to have 2 rows and the same number of columns as cage_cell_
		calc_outward_normal(); // Calculate the outward normal vectors

		return 0; // Return 0 to indicate successful loading of the cage
	}


	//--------------------------------------------------------
	// testing function for load_cage 
	//---------------------------------------------------------

	//void test_load_model_points() {
	//    green_coords_2d instance;
	//    const char* file = "model.obj";

	//    int result = instance.load_model_points(file);
	//    if (result != 0) {
	//        std::cout << "Failed to load model points. Error code: " << result << std::endl;
	//        return;
	//    }

	//    // Verify the loaded data
	//    const Matrix<size_t>& cell = instance.get_cell();
	//    const Matrix<double>& nods = instance.get_nods();

	//    std::cout << "Loaded cell matrix:" << std::endl;
	//    std::cout << cell << std::endl;

	//    std::cout << "Loaded nods matrix:" << std::endl;
	//    std::cout << nods << std::endl;

	//    // Perform tests based on the loaded data
	//    // ...

	//    std::cout << "Model points loaded successfully!" << std::endl;
	//}

	//--------------------------------------------------------------


	int green_coords_2d::calculate_cage_edge_length(VectorXd& len) {
#pragma omp parallel for
		for (size_t i = 0; i < cage_cell_.cols(); ++i) {
			// Calculate the length of each cage edge using the Euclidean norm
			len[i] = (cage_nods_.col(cage_cell_(0, i)) - cage_nods_.col(cage_cell_(1, i))).norm();
		}
		return 0; // Return 0 to indicate successful calculation of cage edge lengths
	}



	//--------------------------------------------------------------------
	// testing function for calculate_cage_edge_length
	//---------------------------------------------------------------------
	//void test_calculate_cage_edge_length() {
	//    green_coords_2d instance;

	//    // Set up test data
	//    Matrix<size_t> cage_cell(2, 3);
	//    cage_cell << 0, 1, 2,
	//        1, 2, 0;

	//    Matrix<double> cage_nods(2, 3);
	//    cage_nods << 0.0, 1.0, 2.0,
	//        0.0, 1.0, 0.0;

	//    instance.set_cage_cell(cage_cell);
	//    instance.set_cage_nods(cage_nods);

	//    VectorXd expected_lengths(3);
	//    expected_lengths << 1.0, 1.41421356, 2.0;

	//    VectorXd calculated_lengths(3);
	//    int result = instance.calculate_cage_edge_length(calculated_lengths);

	//    if (result != 0) {
	//        std::cerr << "Failed to calculate cage edge lengths. Error code: " << result << std::endl;
	//        return;
	//    }

	//    // Verify the calculated lengths
	//    bool lengths_equal = calculated_lengths.isApprox(expected_lengths);

	//    std::cout << "Calculated edge lengths: " << calculated_lengths.transpose() << std::endl;

	//    if (lengths_equal) {
	//        std::cout << "Cage edge lengths calculated successfully!" << std::endl;
	//    }
	//    else {
	//        std::cout << "Calculated cage edge lengths do not match the expected values." << std::endl;
	//    }
	//}

	//-----------------------------------------------------------------------------



	int green_coords_2d::calculate_outward_normal() {
#pragma omp parallel for
		for (size_t i = 0; i < cage_cell_.cols(); ++i) {
			// Calculate the direction vector of the cage edge
			Vector2d dir = cage_nods_.col(cage_cell_(1, i)) - cage_nods_.col(cage_cell_(0, i));
			dir /= dir.norm(); // Normalize the direction vector

			// Calculate the outward normal vector of the cage edge
			cage_normal_(0, i) = -dir[1];
			cage_normal_(1, i) = dir[0];
		}
		return 0; // Return 0 to indicate successful calculation of outward normal vectors
	}


	//-----------------------------------------------------------------------
	// testing function for calculate_outward_normal
	//-----------------------------------------------------------------------

	//void test_calculate_outward_normal() {
	//    green_coords_2d instance;

	//    // Set up test data
	//    Matrix<size_t> cage_cell(2, 3);
	//    cage_cell << 0, 1, 2,
	//        1, 2, 0;

	//    Matrix<double> cage_nods(2, 3);
	//    cage_nods << 0.0, 1.0, 2.0,
	//        0.0, 1.0, 0.0;

	//    instance.set_cage_cell(cage_cell);
	//    instance.set_cage_nods(cage_nods);

	//    Matrix<double> expected_normals(2, 3);
	//    expected_normals << 0.0, 0.70710678, -0.70710678,
	//        1.0, 0.70710678, 0.70710678;

	//    Matrix<double> calculated_normals(2, 3);
	//    int result = instance.calculate_outward_normal(calculated_normals);

	//    if (result != 0) {
	//        std::cerr << "Failed to calculate outward normals. Error code: " << result << std::endl;
	//        return;
	//    }

	//    // Verify the calculated normals
	//    bool normals_equal = calculated_normals.isApprox(expected_normals);

	//    std::cout << "Calculated outward normals:" << std::endl;
	//    std::cout << calculated_normals << std::endl;

	//    if (normals_equal) {
	//        std::cout << "Outward normals calculated successfully!" << std::endl;
	//    }
	//    else {
	//        std::cout << "Calculated outward normals do not match the expected values." << std::endl;
	//    }
	//}

	//--------------------------------------------------------------------



	int green_coords_2d::calculate_green_coords() {
		phi_ = MatrixXd::Zero(cage_nods_.cols(), nods_.cols()); // Initialize phi_ matrix with zeros
		psi_ = MatrixXd::Zero(cage_normal_.cols(), nods_.cols()); // Initialize psi_ matrix with zeros

#pragma omp parallel for
		for (size_t pid = 0; pid < nods_.cols(); ++pid) {
			for (size_t i = 0; i < cage_cell_.cols(); ++i) {
				// Calculate necessary variables for Green coordinates computation
				Vector2d a = cage_nods_.col(cage_cell_(1, i)) - cage_nods_.col(cage_cell_(0, i));
				Vector2d b = cage_nods_.col(cage_cell_(0, i)) - nods_.col(pid);
				double Q, S, R, BA, SRT, L0, L1, A0, A1, A10, L10, psi_entry, phi_entry1, phi_entry0;
				Q = a.dot(a);
				S = b.dot(b);
				R = 2 * a.dot(b);
				BA = b.dot(a.norm() * cage_normal_.col(i));
				SRT = sqrt(4 * S * Q - R * R);
				L0 = log(S);
				L1 = log(S + Q + R);
				A0 = atan2(R, SRT) / SRT;
				A1 = atan2(2 * Q + R, SRT) / SRT;
				A10 = A1 - A0;
				L10 = L1 - L0;
				psi_entry = -a.norm() / (4 * M_PI) * ((4 * S - R * R / Q) * A10 + R / (2 * Q) * L10 + L1 - 2);
				phi_entry1 = -BA / (2 * M_PI) * (L10 / (2 * Q) - A10 * R / Q);
				phi_entry0 = +BA / (2 * M_PI) * (L10 / (2 * Q) - A10 * (2 + R / Q));

#pragma omp critical
				{
					// Accumulate psi, phi entries while ensuring thread safety
					psi_(i, pid) += psi_entry;
					phi_(cage_cell_(1, i), pid) += phi_entry1;
					phi_(cage_cell_(0, i), pid) += phi_entry0;
				}
			}
		}

		// Normalize phi_ columns
		for (size_t j = 0; j < phi_.cols(); ++j) {
			double sum = phi_.col(j).sum();
			phi_.col(j) /= sum;
		}

		return 0; // Return 0 to indicate successful calculation of Green coordinates
	}


	//-------------------------------------------------------------------
	// Testing calculate_green_coords
	//-------------------------------------------------------------------

	//void test_calculate_green_coords() {
	//    green_coords_2d instance;

	//    // Set up test data
	//    Matrix<size_t> cage_cell(2, 3);
	//    cage_cell << 0, 1, 2,
	//        1, 2, 0;

	//    Matrix<double> cage_nods(2, 3);
	//    cage_nods << 0.0, 1.0, 2.0,
	//        0.0, 1.0, 0.0;

	//    Matrix<double> nods(2, 1);
	//    nods << 1.0,
	//        0.5;

	//    instance.set_cage_cell(cage_cell);
	//    instance.set_cage_nods(cage_nods);
	//    instance.set_nods(nods);

	//    MatrixXd expected_phi(3, 1);
	//    expected_phi << -0.2577112,
	//        0.2289056,
	//        0.0288056;

	//    MatrixXd expected_psi(3, 1);
	//    expected_psi << -0.1963495,
	//        0.0,
	//        0.1963495;

	//    int result = instance.calculate_green_coords();

	//    if (result != 0) {
	//        std::cerr << "Failed to calculate Green coordinates. Error code: " << result << std::endl;
	//        return;
	//    }

	//    // Verify the calculated Green coordinates
	//    const MatrixXd& calculated_phi = instance.get_phi();
	//    const MatrixXd& calculated_psi = instance.get_psi();

	//    bool phi_equal = calculated_phi.isApprox(expected_phi);
	//    bool psi_equal = calculated_psi.isApprox(expected_psi);

	//    std::cout << "Calculated phi:" << std::endl;
	//    std::cout << calculated_phi << std::endl;

	//    std::cout << "Calculated psi:" << std::endl;
	//    std::cout << calculated_psi << std::endl;

	//    if (phi_equal && psi_equal) {
	//        std::cout << "Green coordinates calculated successfully!" << std::endl;
	//    }
	//    else {
	//        std::cout << "Calculated Green coordinates do not match the expected values." << std::endl;
	//    }
	//}

	//------------------------------------------------------------------------------------------------





	int green_coords_2d::move_cage(const size_t id, const double* dx, bool disp) {
		Vector2d X(dx[0], dx[1]); // Create a 2D vector from the displacement array

		if (disp)
			cage_nods_.col(id) += X; // If 'disp' is true, add the displacement to the current cage node position
		else
			cage_nods_.col(id) = X; // If 'disp' is false, set the cage node position to the displacement

		return 0; // Return 0 to indicate successful movement of the cage node
	}


	//-----------------------------------------------------------------------
	// testing move_cage function 
	 //---------------------------------------------------------------------

	 //void test_move_cage() {
	 //    green_coords_2d instance;

	 //    // Set up test data
	 //    Matrix<double> cage_nods(2, 3);
	 //    cage_nods << 0.0, 1.0, 2.0,
	 //        0.0, 1.0, 0.0;

	 //    instance.set_cage_nods(cage_nods);

	 //    size_t id = 1;
	 //    double dx[2] = { 0.5, 0.5 };
	 //    bool disp = true;

	 //    int result = instance.move_cage(id, dx, disp);

	 //    if (result != 0) {
	 //        std::cerr << "Failed to move cage. Error code: " << result << std::endl;
	 //        return;
	 //    }

	 //    // Verify the moved cage nods
	 //    const Matrix<double>& updated_cage_nods = instance.get_cage_nods();

	 //    std::cout << "Updated cage nods:" << std::endl;
	 //    std::cout << updated_cage_nods << std::endl;

	 //    
	 //    std::cout << "Cage moved successfully!" << std::endl;
	 //}

	 //------------------------------------------------------------------------


	int green_coords_2d::deform() {
		calc_outward_normal(); // Calculate the outward normal vectors of the cage elements
		calculate_cage_edge_length(curr_len_); // Calculate the current segment lengths of the cage edges
		VectorXd ratio = curr_len_.cwiseQuotient(rest_len_); // Calculate the ratio between current and rest segment lengths

		if (cage_normal_.cols() != ratio.rows()) {
			std::cerr << "Error: cage_normal and ratio dimensions do not match." << std::endl;
			return __LINE__; // Return the line number as an error code if dimensions do not match
		}

		nods_ = cage_nods_ * phi_ + (cage_normal_ * ratio.asDiagonal()) * psi_;
		// Deform the nodal positions using Green coordinates, cage node positions, cage normals, and ratio
		// The new nodal positions are calculated as a combination of the weighted cage node positions and the weighted cage normals

		return 0; // Return 0 to indicate successful deformation
	}


	//----------------------------------------------------
	// Test for deform function 
	//-----------------------------------------------------

	//void test_deform() {
	//    green_coords_2d instance;

	//    // Set up test data
	//    Matrix<size_t> cage_cell(2, 3);
	//    cage_cell << 0, 1, 2,
	//        1, 2, 0;

	//    Matrix<double> cage_nods(2, 3);
	//    cage_nods << 0.0, 1.0, 2.0,
	//        0.0, 1.0, 0.0;

	//    Matrix<double> nods(2, 1);
	//    nods << 1.0,
	//        0.5;

	//    VectorXd rest_len(3);
	//    rest_len << 1.0, 1.0, 1.0;

	//    MatrixXd phi(3, 1);
	//    phi << 0.1,
	//        0.2,
	//        0.3;

	//    MatrixXd psi(3, 1);
	//    psi << 0.4,
	//        0.5,
	//        0.6;

	//    instance.set_cage_cell(cage_cell);
	//    instance.set_cage_nods(cage_nods);
	//    instance.set_nods(nods);
	//    instance.set_rest_len(rest_len);
	//    instance.set_phi(phi);
	//    instance.set_psi(psi);

	//    int result = instance.deform();

	//    if (result != 0) {
	//        std::cerr << "Failed to deform. Error code: " << result << std::endl;
	//        return;
	//    }

	//    // Verify the deformed nods
	//    const Matrix<double>& deformed_nods = instance.get_nods();

	//    std::cout << "Deformed nods:" << std::endl;
	//    std::cout << deformed_nods << std::endl;

	//    // Perform additional assertions or tests based on the deformed nods
	//    // ...

	//    std::cout << "Deformation successful!" << std::endl;
	//}

	//-------------------------------------------------------------------------



	int green_coords_2d::saveToFile(const char* file) {
		MatrixXd nods_3d(3, nods_.cols());
#pragma omp parallel for
		for (size_t i = 0; i < nods_.cols(); ++i) {
			nods_3d(0, i) = nods_(0, i);
			nods_3d(1, i) = 0.0;
			nods_3d(2, i) = nods_(1, i);
		}

		std::ofstream os(file);
		if (!os) {
			std::cerr << "Error: Failed to open file for writing: " << file << std::endl;
			return __LINE__;
		}

		// Save as OBJ format
		// Write vertex coordinates
		for (size_t i = 0; i < nods_3d.cols(); ++i) {
			os << "v " << nods_3d(0, i) << " " << nods_3d(1, i) << " " << nods_3d(2, i) << std::endl;
		}

		// Write faces
		for (size_t i = 0; i < cell_.cols(); ++i) {
			os << "f";
			for (size_t j = 0; j < cell_.rows(); ++j) {
				os << " " << cell_(j, i) + 1;
			}
			os << std::endl;
		}

		os.close();

		// Save as TRI format

		/*
		std::ofstream os(file);
		if (!os) {
			std::cerr << "Error: Failed to open file for writing: " << file << std::endl;
			return __LINE__;
		}

		// Write vertex coordinates
		for (size_t i = 0; i < nods_3d.cols(); ++i) {
			os << nods_3d(0, i) << " " << nods_3d(1, i) << " " << nods_3d(2, i) << std::endl;
		}

		// Write faces
		for (size_t i = 0; i < cell_.cols(); ++i) {
			os << cell_(0, i) + 1 << " " << cell_(1, i) + 1 << " " << cell_(2, i) + 1 << std::endl;
		}

		os.close();
		*/

		return 0;
	}

	//---------------------------------------------
	// testing function saveToFile 
	//-------------------------------------------

	//void test_saveToFile() {
	//    green_coords_2d instance;

	//    // Set up test data
	//    MatrixXd nods(2, 4);
	//    nods << 0.0, 1.0, 1.0, 0.0,
	//        0.0, 0.0, 1.0, 1.0;

	//    Matrix<size_t> cell(3, 2);
	//    cell << 0, 1,
	//        1, 2,
	//        2, 3;

	//    instance.set_nods(nods);
	//    instance.set_cell(cell);

	//    const char* file = "mesh.obj";

	//    int result = instance.saveToFile(file);

	//    if (result != 0) {
	//        std::cerr << "Failed to save mesh data to file. Error code: " << result << std::endl;
	//        return;
	//    }

	//    std::cout << "Mesh data saved to file: " << file << std::endl;
	//}
//--------------------------------------------------------------------------------------




	int green_coords_2d::saveToFile_cage(const char* file) {
		MatrixXd cage_nods_3d(3, cage_nods_.cols());
#pragma omp parallel for
		for (size_t i = 0; i < cage_nods_.cols(); ++i) {
			cage_nods_3d(0, i) = cage_nods_(0, i);
			cage_nods_3d(1, i) = 0.0;
			cage_nods_3d(2, i) = cage_nods_(1, i);
		}

		std::ofstream os(file);
		if (!os) {
			std::cerr << "Error: Failed to open file for writing: " << file << std::endl;
			return __LINE__;
		}

		// Save as OBJ format
		// Write vertex coordinates
		for (size_t i = 0; i < cage_nods_3d.cols(); ++i) {
			os << "v " << cage_nods_3d(0, i) << " " << cage_nods_3d(1, i) << " " << cage_nods_3d(2, i) << std::endl;
		}

		// Write lines
		for (size_t i = 0; i < cage_cell_.cols(); ++i) {
			os << "l " << cage_cell_(0, i) + 1 << " " << cage_cell_(1, i) + 1 << std::endl;
		}

		os.close();

		// Save as TRI format
		/*
		std::ofstream os(file);
		if (!os) {
			std::cerr << "Error: Failed to open file for writing: " << file << std::endl;
			return __LINE__;
		}

		// Write vertex coordinates
		for (size_t i = 0; i < cage_nods_3d.cols(); ++i) {
			os << cage_nods_3d(0, i) << " " << cage_nods_3d(1, i) << " " << cage_nods_3d(2, i) << std::endl;
		}

		// Write lines
		for (size_t i = 0; i < cage_cell_.cols(); ++i) {
			os << cage_cell_(0, i) + 1 << " " << cage_cell_(1, i) + 1 << std::endl;
		}

		os.close();
		*/

		return 0;
	}



	//-------------------------------------------------------------
	// testing functiom saveToFile_cage 
	//-------------------------------------------------------------

	//void test_saveToFile_cage() {
	//    green_coords_2d instance;

	//    // Set up test data
	//    MatrixXd cage_nods(2, 4);
	//    cage_nods << 0.0, 1.0, 1.0, 0.0,
	//        0.0, 0.0, 1.0, 1.0;

	//    Matrix<size_t> cage_cell(2, 2);
	//    cage_cell << 0, 1,
	//        1, 2;

	//    instance.set_cage_nods(cage_nods);
	//    instance.set_cage_cell(cage_cell);

	//    const char* file = "cage.obj";

	//    int result = instance.saveToFile_cage(file);

	//    if (result != 0) {
	//        std::cerr << "Failed to save cage data to file. Error code: " << result << std::endl;
	//        return;
	//    }

	//    std::cout << "Cage data saved to file: " << file << std::endl;
	//}
//----------------------------------------------------------------------


	int green_coords_2d::saveToFile_normal(const char* file) {
		Matrix<size_t, Dynamic, Dynamic> normal_cell(2, cage_cell_.cols());
		for (size_t i = 0; i < normal_cell.size(); ++i)
			normal_cell(i) = i;

		Matrix<double, 3, Dynamic> normal_nods(3, 2 * normal_cell.size());
		normal_nods.setZero();


#pragma omp parallel for
		for (size_t i = 0; i < normal_nods.cols() / 2; ++i) {
			Vector2d mid = 0.5 * (cage_nods_.col(cage_cell_(0, i)) + cage_nods_.col(cage_cell_(1, i)));
			Vector2d end = mid + cage_normal_.col(i);
			normal_nods(0, 2 * i + 0) = mid[0];
			normal_nods(2, 2 * i + 0) = mid[1];
			normal_nods(0, 2 * i + 1) = end[0];
			normal_nods(2, 2 * i + 1) = end[1];
		}

		std::ofstream os(file);
		if (!os) {
			std::cerr << "Error: Failed to open file for writing: " << file << std::endl;
			return __LINE__;
		}

		// Save as OBJ format
		// Write vertex coordinates
		for (size_t i = 0; i < normal_nods.cols(); ++i) {
			os << "v " << normal_nods(0, i) << " " << normal_nods(1, i) << " " << normal_nods(2, i) << std::endl;
		}

		// Write lines
		for (size_t i = 0; i < normal_cell.cols(); ++i) {
			os << "l " << 2 * i + 1 << " " << 2 * i + 2 << std::endl;
		}

		os.close();

		// Save as TRI format
		// Uncomment the following code to save as TRI format
		/*
		std::ofstream os(file);
		if (!os) {
			std::cerr << "Error: Failed to open file for writing: " << file << std::endl;
			return __LINE__;
		}

		// Write vertex coordinates
		for (size_t i = 0; i < normal_nods.cols(); ++i) {
			os << normal_nods(0, i) << " " << normal_nods(1, i) << " " << normal_nods(2, i) << std::endl;
		}

		// Write lines
		for (size_t i = 0; i < normal_cell.cols(); ++i) {
			os << 2 * i + 1 << " " << 2 * i + 2 << std::endl;
		}

		os.close();
		*/

		return 0;
	}




	//------------------------------------------------------
   // Testing function for  saveToFile_normal 
   //--------------------------------------------------

   //void test_saveToFile_normal() {
   //    green_coords_2d instance;

   //    // Set up test data
   //    Matrix<size_t> cage_cell(2, 3);
   //    cage_cell << 0, 1, 2,
   //        1, 2, 0;

   //    Matrix<double> cage_nods(2, 3);
   //    cage_nods << 0.0, 1.0, 2.0,
   //        0.0, 1.0, 0.0;

   //    instance.set_cage_cell(cage_cell);
   //    instance.set_cage_nods(cage_nods);

   //    const char* file = "normal.obj";

   //    int result = instance.saveToFile_normal(file);

   //    if (result != 0) {
   //        std::cerr << "Failed to save normal data to file. Error code: " << result << std::endl;
   //        return;
   //    }

   //    std::cout << "Normal data saved to file: " << file << std::endl;
   //}


   //------------------------------------------------------------------
	////////////////////////////////////////////////////////////////////////////////////////////////
	// implementing green coordinates in 3D , however some function are commented
	// due too errorrs

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	green_coords_3d::green_coords_3d() {}

	int green_coords_3d::load_model_points(const char* file) {
		// Open the input file
		std::ifstream inputFile(file);
		if (!inputFile.is_open())
			return __LINE__;

		// Temporary storage for vertices and faces
		std::vector<std::vector<double>> tempNods;
		std::vector<std::vector<size_t>> tempCell;

		// Read the file line by line
		std::string line;
		while (std::getline(inputFile, line)) {
			std::stringstream ss(line);
			std::string identifier;
			ss >> identifier;

			// Parse vertex information
			if (identifier == "v") {
				double x, y, z;
				ss >> x >> y >> z;
				tempNods.push_back({ x, y, z });
			}
			// Parse face information
			else if (identifier == "f") {
				size_t v1, v2, v3;
				ss >> v1 >> v2 >> v3;
				tempCell.push_back({ v1 - 1, v2 - 1, v3 - 1 }); // Adjust indices to 0-based
			}
		}

		// Resize the nods_ Eigen matrix and transfer data from temporary vectors
		nods_.resize(tempNods.size(), 3);
		for (size_t i = 0; i < tempNods.size(); ++i) {
			nods_(i, 0) = tempNods[i][0];
			nods_(i, 1) = tempNods[i][1];
			nods_(i, 2) = tempNods[i][2];
		}

		// Resize the cell_ Eigen matrix and transfer data from temporary vectors
		cell_.resize(tempCell.size(), 3);
		for (size_t i = 0; i < tempCell.size(); ++i) {
			cell_(i, 0) = tempCell[i][0];
			cell_(i, 1) = tempCell[i][1];
			cell_(i, 2) = tempCell[i][2];
		}

		// Close the input file
		inputFile.close();

		// Return success
		return 0;
	}





	int green_coords_3d::load_cage(const char* file) {
		// Parse the file manually
		ifstream is(file);
		if (is.fail()) {
			cerr << "# can not open " << file << "\n";
			return __LINE__; // Unable to open the file, return the line number as an error code
		}

		// Read the number of cage elements
		size_t ele_num;
		string CELL;
		is >> CELL >> ele_num;
		cage_cell_.resize(3, ele_num); // Resize cage_cell_ matrix to have 3 rows and ele_num columns

		// Read the cage element vertices
		for (size_t i = 0; i < ele_num; ++i)
			is >> cage_cell_(0, i) >> cage_cell_(1, i) >> cage_cell_(2, i); // Store vertex coordinates of each cage element

		// Read the number of cage nodes
		size_t nods_num;
		string NODES;
		is >> NODES >> nods_num;
		cage_nods_.resize(3, nods_num); // Resize cage_nods_ matrix to have 3 rows and nods_num columns

		// Read the cage node coordinates
		for (size_t i = 0; i < nods_num; ++i)
			is >> cage_nods_(0, i) >> cage_nods_(1, i) >> cage_nods_(2, i); // Store coordinates of each cage node

		is.close(); // Close the input file

		// Get the rest segment element length and initialize current segment lengths
		cage_rest_uv_.resize(cage_cell_.cols()); // Resize rest_len_ vector to have the same number of elements as cage_cell_ columns
		cage_curr_uv_.resize(cage_cell_.cols()); // Resize curr_len_ vector to have the same number of elements as cage_cell_ columns
		calculate_cage_edge(cage_rest_uv_); // Calculate the rest segment element lengths

		// Calculate the initial outward normal vectors of the cage elements
		cage_normal_.resize(3, cage_cell_.cols()); // Resize cage_normal_ matrix to have 3 rows and the same number of columns as cage_cell_
		calc_outward_normal(); // Calculate the outward normal vectors

		return 0; // Return 0 to indicate successful loading of the cage
	}





	int green_coords_3d::calculate_cage_edge(matd_t& uv) {
		size_t numFaces = cage_cell_.cols();
		size_t numEdges = 2 * numFaces;

		// Resize the uv matrix if necessary
		if (uv.rows() != 3 || uv.cols() != numEdges)
			uv.resize(3, numEdges);

#pragma omp parallel for
		for (size_t i = 0; i < numFaces; ++i) {
			size_t edgeIndex = 2 * i;
			size_t v0Index = cage_cell_(0, i);
			size_t v1Index = cage_cell_(1, i);
			size_t v2Index = cage_cell_(2, i);

			// Calculate the first edge vector
			uv(0, edgeIndex) = cage_nods_(v1Index, 0) - cage_nods_(v0Index, 0);
			uv(1, edgeIndex) = cage_nods_(v1Index, 1) - cage_nods_(v0Index, 1);
			uv(2, edgeIndex) = cage_nods_(v1Index, 2) - cage_nods_(v0Index, 2);

			edgeIndex++;

			// Calculate the second edge vector
			uv(0, edgeIndex) = cage_nods_(v2Index, 0) - cage_nods_(v1Index, 0);
			uv(1, edgeIndex) = cage_nods_(v2Index, 1) - cage_nods_(v1Index, 1);
			uv(2, edgeIndex) = cage_nods_(v2Index, 2) - cage_nods_(v1Index, 2);
		}

		return 0;
	}



	int green_coords_3d::calculate_stretch_ratio(matd_t& s) {
		size_t numFaces = cage_normal_.cols();
		s.resize(numFaces, 1);
		calculate_cage_edge(cage_curr_uv_);
		assert(cage_curr_uv_.cols() == 2 * s.rows());

#pragma omp parallel for
		for (size_t i = 0; i < numFaces; ++i) {
			matd_t u = cage_rest_uv_.col(2 * i);
			matd_t v = cage_rest_uv_.col(2 * i + 1);
			matd_t ucap = cage_curr_uv_.col(2 * i);
			matd_t vcap = cage_curr_uv_.col(2 * i + 1);
			// Calculate dot products
			double dotUcapUcap = ucap.dot(ucap);
			double dotV = v.dot(v);
			double dotUcapVcap = ucap.dot(vcap);
			double dotU = u.dot(u);
			double dotVcapVcap = vcap.dot(vcap);
			double dotUcapV = ucap.dot(v);
			// Calculate stretch ratio
			s(i) = std::sqrt(dotUcapUcap * dotV - 2 * dotUcapVcap * dotUcapV + dotVcapVcap * dotU)
				/ (std::sqrt(8.0) * rest_area_(i));
		}

		return 0;
	}



	// Function to calculate the cross product of two 3D vectors
	matd_t cross(const matd_t& v1, const matd_t& v2) {
		matd_t result(3, 1);
		result(0) = v1(1) * v2(2) - v1(2) * v2(1);
		result(1) = v1(2) * v2(0) - v1(0) * v2(2);
		result(2) = v1(0) * v2(1) - v1(1) * v2(0);
		return result;
	}

	int green_coords_3d::calculate_outward_normal() {
		cage_normal_.resize(3, cage_cell_.cols());

		for (size_t i = 0; i < cage_cell_.cols(); ++i) {
			// Get the vertices of the current cage element
			matd_t vt(3, 3);
			for (size_t j = 0; j < 3; ++j) {
				size_t vertex_index = cage_cell_(j, i);
				vt(0, j) = cage_nods_(0, vertex_index);
				vt(1, j) = cage_nods_(1, vertex_index);
				vt(2, j) = cage_nods_(2, vertex_index);
			}

			// Calculate the face normal using cross product
			matd_t v1 = vt.col(1) - vt.col(0);
			matd_t v2 = vt.col(2) - vt.col(0);
			matd_t face_normal = cross(v1, v2);

			// Normalize the face normal
			double norm = face_normal.norm();
			if (norm > 0)
				face_normal /= norm;

			// Store the outward normal vector
			for (size_t j = 0; j < 3; ++j)
				cage_normal_(j, i) = face_normal(j);
		}

		return 0;
	}




	int green_coords_3d::saveToFile(const char* file) {
		Eigen::MatrixXd nods_3d(3, nods_.cols());
#pragma omp parallel for
		for (size_t i = 0; i < nods_.cols(); ++i) {
			nods_3d(0, i) = nods_(0, i);
			nods_3d(1, i) = 0.0;
			nods_3d(2, i) = nods_(1, i);
		}

		std::ofstream os(file);
		if (!os) {
			std::cerr << "Error: Failed to open file for writing: " << file << std::endl;
			return __LINE__;
		}

		// Save as OBJ format
		// Write vertex coordinates
		for (size_t i = 0; i < nods_3d.cols(); ++i) {
			os << "v " << nods_3d(0, i) << " " << nods_3d(1, i) << " " << nods_3d(2, i) << std::endl;
		}

		// Write faces
		for (size_t i = 0; i < cell_.cols(); ++i) {
			os << "f";
			for (size_t j = 0; j < cell_.rows(); ++j) {
				os << " " << cell_(j, i) + 1;
			}
			os << std::endl;
		}

		os.close();

		//// Save as TRI format
		//std::string triFileName(file);
		//triFileName += ".tri";
		//std::ofstream triOs(triFileName);
		//if (!triOs) {
		//    std::cerr << "Error: Failed to open TRI file for writing: " << triFileName << std::endl;
		//    return __LINE__;
		//}

		//// Write vertex coordinates
		//for (size_t i = 0; i < nods_3d.cols(); ++i) {
		//    triOs << nods_3d(0, i) << " " << nods_3d(1, i) << " " << nods_3d(2, i) << std::endl;
		//}

		//// Write faces
		//for (size_t i = 0; i < cell_.cols(); ++i) {
		//    triOs << cell_(0, i) + 1 << " " << cell_(1, i) + 1 << " " << cell_(2, i) + 1 << std::endl;
		//}

		//triOs.close();

		return 0;
	}




	int green_coords_3d::saveToFile_cage(const char* file) {
		ofstream objOs(file);
		if (!objOs) {
			cerr << "Error: Failed to open OBJ file for writing: " << file << endl;
			return __LINE__;
		}

		// Write vertex coordinates
		for (size_t i = 0; i < cage_nods_.size(); ++i) {
			objOs << "v " << cage_nods_(0, i) << " " << cage_nods_(1, i) << " " << cage_nods_(2, i) << endl;
		}

		// Write faces
		for (size_t i = 0; i < cage_cell_.size(); ++i) {
			objOs << "f";
			for (size_t j = 0; j < cage_cell_.rows(); ++j) {
				objOs << " " << cage_cell_(j, i) + 1;
			}
			objOs << endl;
		}

		objOs.close();

		return 0;
	}


	//
	//	int green_coords_3d::saveToFile_normal(const char* file) {
	//		mati_t normal_cell(2, cage_cell_.size(2));
	//		matd_t normal_nods(3, 2 * cage_cell_.size(2));
	//#pragma omp parallel for
	//		for (size_t i = 0; i < normal_cell.size(2); ++i) {
	//			normal_cell(0, i) = 2 * i + 0;
	//			normal_cell(1, i) = 2 * i + 1;
	//			normal_nods(colon(), normal_cell(0, i)) = 1.0 / 3 * cage_nods_(colon(), cage_cell_(colon(), i)) * ones<double>(3, 1);
	//			normal_nods(colon(), normal_cell(1, i)) = normal_nods(colon(), normal_cell(0, i)) + cage_normal_(colon(), i);
	//		}
	//
	//		std::ofstream os(file);
	//		if (!os) {
	//			std::cerr << "Error: Failed to open OBJ file for writing: " << file << std::endl;
	//			return __LINE__;
	//		}
	//
	//		// Write vertex coordinates
	//		for (size_t i = 0; i < normal_nods.size(2); ++i) {
	//			os << "v " << normal_nods(0, i) << " " << normal_nods(1, i) << " " << normal_nods(2, i) << std::endl;
	//		}
	//
	//		// Write faces
	//		for (size_t i = 0; i < normal_cell.size(2); ++i) {
	//			os << "f " << normal_cell(0, i) + 1 << " " << normal_cell(1, i) + 1 << std::endl;
	//		}
	//
	//		os.close();
	//		return 0;
	//	}


		//double green_coords_3d::GCTriInt(const matd_t& p, const matd_t& v1, const matd_t& v2, const matd_t& eta) {
		//	// Calculate the angles and distances
		//	double alpha = std::acos(std::min(1.0, std::max(-1.0, dot(v2 - v1, p - v1) / norm(v2 - v1) / norm(p - v1))));
		//	double beta = std::acos(std::min(1.0, std::max(-1.0, dot(v1 - p, v2 - p) / norm(v1 - p) / norm(v2 - p))));
		//	double lambda = dot(p - v1, p - v1) * std::sin(alpha) * std::sin(alpha);
		//	double c = dot(p - eta, p - eta);

		//	// Calculate the angles and constants used in the integral
		//	double theta[2] = { M_PI - alpha, M_PI - alpha - beta };
		//	double I[2];
		//	double sqrt_c = std::sqrt(c);
		//	double sqrt_lambda = std::sqrt(lambda);
		//	for (int i = 0; i < 2; ++i) {
		//		double S = std::sin(theta[i]);
		//		double C = std::cos(theta[i]);

		//		// Calculate the integral terms
		//		I[i] = -std::copysign(1.0, S) / 2.0 * (2.0 * sqrt_c * std::atan2(sqrt_c * C, std::sqrt(lambda + S * S * c))
		//			+ sqrt_lambda * std::log(2.0 * sqrt_lambda * S * S / ((1 - C) * (1 - C))
		//				* (1.0 - 2 * c * C / (c + c * C + lambda + std::sqrt(lambda * lambda + lambda * c * S * S)))));
		//	}

		//	// Calculate the Green Coordinates integral
		//	return -1.0 / (4 * M_PI) * std::fabs(I[0] - I[1] - sqrt_c * beta);
		//}



	//
	//	int green_coords_3d::calculate_green_coords() {
	//		matd_t zero3d = zeros<double>(3, 1);
	//		phi_ = zeros<double>(cage_nods_.size(2), nods_.size(2));
	//		psi_ = zeros<double>(cage_normal_.size(2), nods_.size(2));
	//
	//		for (size_t pid = 0; pid < nods_.size(2); ++pid) {
	//			for (size_t i = 0; i < cage_cell_.size(2); ++i) {
	//				// Calculate the relative positions and normals
	//				matd_t VJ = cage_nods_(colon(), cage_cell_(colon(), i)) - nods_(colon(), pid) * ones<double>(1, 3);
	//				matd_t p = dot(VJ(colon(), 0), cage_normal_(colon(), i)) * cage_normal_(colon(), i);
	//				matd_t s(3, 1), I(3, 1), II(3, 1), N(3, 3);
	//
	//				for (size_t j = 0; j < 3; ++j) {
	//					// Calculate the sign, Green Coordinates integrals, and normals
	//					s[j] = sign(dot(cross(VJ(colon(), j) - p, VJ(colon(), (j + 1) % 3) - p), cage_normal_(colon(), i)));
	//					I[j] = GCTriInt(p, VJ(colon(), j), VJ(colon(), (j + 1) % 3), zero3d);
	//					II[j] = GCTriInt(zero3d, VJ(colon(), (j + 1) % 3), VJ(colon(), j), zero3d);
	//					matd_t q = cross(VJ(colon(), (j + 1) % 3), VJ(colon(), j));
	//					N(colon(), j) = q / norm(q);
	//				}
	//
	//				double I_ = -fabs(dot(s, I));
	//				psi_(i, pid) += -I_;
	//				matd_t w = I_ * cage_normal_(colon(), i) + N * II;
	//
	//				if (norm(w) > 1e-8) {
	//					// Calculate the weighted contributions to phi_
	//					for (size_t j = 0; j < 3; ++j)
	//						phi_(cage_cell_(j, i), pid) += dot(N(colon(), (j + 1) % 3), w) / dot(N(colon(), (j + 1) % 3), VJ(colon(), j));
	//				}
	//			}
	//		}
	//
	//		// Perform normalization for translation invariance
	//#pragma omp parallel for
	//		for (size_t j = 0; j < phi_.size(2); ++j) {
	//			double col_sum = sum(phi_(colon(), j));
	//			phi_(colon(), j) /= col_sum;
	//		}
	//
	//		return 0;
	//	}






		//int green_coords_3d::move_cage(const size_t id, const double* dx, bool disp) {
		//	// Create a matrix view of the displacement vector
		//	itr_matrix<const double*> X(3, 1, dx);
		//	// Move the cage node by adding the displacement vector
		//	if (disp)
		//		cage_nods_(colon(), id) += X;
		//	else
		//		cage_nods_(colon(), id) = X;
		//	return 0;
		//}


		//int green_coords_3d::deform() {
		//	// Calculate the outward normal vectors of the cage elements
		//	calc_outward_normal();

		//	// Calculate the stretch ratios
		//	matd_t s;
		//	calc_stretch_ratio(s);
		//	ASSERT(s.size() == cage_normal_.size(2));

		//	// Apply the stretch ratios to psi_
		//	matd_t stretch_psi = psi_;
		//	for (size_t row = 0; row < stretch_psi.size(1); ++row)
		//		stretch_psi(row, colon()) *= s[row];

		//	// Calculate the deformed positions of the mesh nodes
		//	nods_ = cage_nods_ * phi_ + cage_normal_ * stretch_psi;

		//	return 0;
		//}













}