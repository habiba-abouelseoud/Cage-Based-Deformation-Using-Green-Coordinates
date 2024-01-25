//cage_deform_gc_2D_3D.h
//Copyright © Habibatallah Abuelseoud.
//Created by Habibatallah Abouelseoud 19/11/2022.

/**@file
cage_deform_gc_2D_3D.h
The header file for source file of the cage_deform_gc_2D , in which the code structure uses the inheritance
as class can be derived from more than one classes, which means it can inherit data and functions
from multiple base classes.this file conatins the defintions  of the Parent class.
*/

#ifndef CAGE_DEFORM_GC_2D
#define CAGE_DEFORM_GC_2D
#include <Eigen/Sparse>
namespace GreenCoordinates {

	/**
	 * @brief Abstract base class for Green coordinates deformation.
	 */

	class green_coordinates
	{
	public:
		/**
		 * @brief Default constructor for green_coordinates.
		 */
		green_coordinates() {}

		/**
		 * @brief Virtual destructor for green_coordinates.
		 */
		virtual ~green_coordinates() {}

		/**
		 * @brief Load model points from a file.
		 * @param file The file path to load the model points from.
		 * @return The result of the operation (0 for success, non-zero for failure).
		 */
		virtual int load_model_points(const char* file) = 0;

		/**
		 * @brief Load the cage from a file.
		 * @param file The file path to load the cage from.
		 * @return The result of the operation (0 for success, non-zero for failure).
		 */
		virtual int load_cage(const char* file) = 0;

		/**
		 * @brief Calculate the Green coordinates.
		 * @return The result of the operation (0 for success, non-zero for failure).
		 */
		virtual int calculate_green_coords() = 0;

		/**
		 * @brief Move a cage node by a given displacement.
		 * @param id The ID of the cage node to move.
		 * @param dx The displacement vector [dx, dy].
		 * @param disp Whether to apply the displacement (true) or set the position directly (false).
		 * @return The result of the operation (0 for success, non-zero for failure).
		 */
		virtual int move_cage(const size_t id, const double* dx, bool disp) = 0;

		/**
		 * @brief Perform the deformation using Green coordinates.
		 * @return The result of the operation (0 for success, non-zero for failure).
		 */
		virtual int deform() = 0;

		/**
		 * @brief Save the deformed model to a file.
		 * @param file The file path to save the deformed model to.
		 * @return The result of the operation (0 for success, non-zero for failure).
		 */
		virtual int saveToFile(const char* file) = 0;

	protected:
		/**
		 * @brief Calculate the outward normal vectors of the cage elements.
		 * @return The result of the operation (0 for success, non-zero for failure).
		 */
		virtual int calc_outward_normal() = 0;
	};

	/**
	 * @brief Concrete class for 2D Green coordinates deformation.
	 */
	class green_coords_2d : public green_coordinates
	{
	public:
		/**
		 * @brief Default constructor for green_coords_2d.
		 */
		green_coords_2d();

		/**
		* @brief Load model points from a file.
		*
		* This function reads nodal and cell data from a file and populates the corresponding data containers.
		* The file format should follow the "v" (nodal coordinates) and "f" (cell connectivity) keywords.
		*
		* @param file The file path to load the model points from.
		* @return The result of the operation (0 for success, non-zero for failure).
		*         In case of failure, the line number where the error occurred is returned.
		*/
		int load_model_points(const char* file);

		/**
		 * @brief Load the cage from a file.
		 *
		 * This function manually parses the file to extract the cage element and node information.
		 * The file should contain the number of cage elements, followed by the cage element vertices,
		 * and then the number of cage nodes and their coordinates.
		 *
		 * @param file The file path to load the cage from.
		 * @return The result of the operation (0 for success, non-zero for failure).
		 *         In case of failure, the line number where the error occurred is returned.
		 */

		int load_cage(const char* file);


		/**
	   * @brief Calculate the Green coordinates.
	   *
	   * This function computes the Green coordinates for each point in the `nods_` matrix.
	   * It iterates over each point and each cage element to calculate the necessary variables,
	   * such as vector differences and logarithmic terms, to compute the Green coordinates.
	   * The resulting Green coordinates are stored in the `phi_` and `psi_` matrices.
	   *
	   * @return The result of the operation (0 for success, non-zero for failure).
	   */

		int calculate_green_coords();

		/**
		 * @brief Move a cage node by a given displacement.
		 *
		 * This function moves a cage node identified by its index (`id`) by a specified displacement (`dx`).
		 * The displacement is provided as a pointer to a double array containing the x and y components.
		 * If the `disp` parameter is true, the displacement is added to the current cage node position.
		 * If `disp` is false, the cage node position is set to the displacement.
		 *
		 * @param id The index of the cage node to move.
		 * @param dx Pointer to a double array containing the displacement in x and y components.
		 * @param disp Flag indicating whether to add the displacement to the current position or set it as the new position.
		 * @return The result of the operation (0 for success, non-zero for failure).
		 */
		int move_cage(const size_t id, const double* dx, bool disp);

		/**
		 * @brief Deform the model points using Green coordinates and cage information.
		 *
		 * This function deforms the model points based on the calculated Green coordinates, cage node positions,
		 * cage normals, and the ratio between current and rest segment lengths. It performs the following steps:
		 *
		 * 1. Calculate the outward normal vectors of the cage elements.
		 * 2. Calculate the current segment lengths of the cage edges.
		 * 3. Calculate the ratio between the current and rest segment lengths.
		 * 4. Check if the dimensions of the cage normals and the ratio match. If they don't, return an error code.
		 * 5. Deform the nodal positions using the Green coordinates, cage node positions, cage normals, and ratio.
		 *    The new nodal positions are calculated as a combination of the weighted cage node positions and the
		 *    weighted cage normals.
		 *
		 * @return The result of the operation (0 for success, non-zero for failure).
		 *         In case of failure, the line number where the error occurred is returned.
		 */
		int deform();

		/**
		 * @brief Save the normal information to a file.
		 *
		 * This function saves the midpoints and endpoints of the cage edges with the corresponding outward normals to a file.
		 * The file can be saved in either OBJ or TRI format.
		 *
		 * @param file The file path to save the normal information to.
		 * @return The result of the operation (0 for success, non-zero for failure).
		 *         In case of failure, the line number where the error occurred is returned.
		 */

		int saveToFile(const char* file);

		/**
		 * @brief Save the normal information to a file.
		 *
		 * This function saves the midpoints and endpoints of the cage edges with the corresponding outward normals to a file.
		 * The file can be saved in either OBJ or TRI format.
		 *
		 * @param file The file path to save the normal information to.
		 * @return The result of the operation (0 for success, non-zero for failure).
		 *         In case of failure, the line number where the error occurred is returned.
		 */
		int saveToFile_normal(const char* file);

		/**
		* @brief Save the cage information to a file.
		*
		* This function saves the cage node positions and the cage connectivity information to a file.
		* The file can be saved in either OBJ or TRI format.
		*
		* @param file The file path to save the cage information to.
		* @return The result of the operation (0 for success, non-zero for failure).
		*         In case of failure, the line number where the error occurred is returned.
		*/

		int saveToFile_cage(const char* file);


	protected:
		/**
		 * @brief Calculate the outward normal vectors of the cage edges.
		 *
		 * This function calculates the outward normal vector for each cage edge.
		 * It first calculates the direction vector of each cage edge and then normalizes it.
		 * The resulting outward normal vectors are stored in the `cage_normal_` matrix.
		 *
		 * @return The result of the operation (0 for success, non-zero for failure).
		 */
		int calculate_outward_normal();

		/**
		* @brief Calculate the lengths of the cage edges.
		*
		* This function calculates the length of each cage edge using the Euclidean norm.
		* The result is stored in the provided `len` vector.
		*
		* @param len The vector to store the calculated edge lengths.
		* @return The result of the operation (0 for success, non-zero for failure).
		*
		*/

		int calculate_cage_edge_length(Eigen::VectorXd& len);

	private:
		Eigen::MatrixXi cell_;
		Eigen::MatrixXd nods_;

		Eigen::MatrixXi cage_cell_;
		Eigen::MatrixXd cage_nods_;
		Eigen::MatrixXd cage_normal_;

		Eigen::VectorXd rest_len_;
		Eigen::VectorXd curr_len_;
		Eigen::MatrixXd phi_;
		Eigen::MatrixXd psi_;
	};






	/**
	 * @brief Concrete class for 3D Green coordinates deformation.
	 */
	class green_coords_3d : public green_coordinates
	{
	public:
		typedef Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mati_t;
		typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matd_t;


		/**
		 * @brief Default constructor for green_coords_3d.
		 */
		green_coords_3d();

		/**
		 * @brief Load the model points from a file.
		 *
		 * This function loads the model points from a file and populates the nods_ and cell_ matrices.
		 * It reads the file line by line and parses the vertex and face information.
		 * Vertex lines start with "v" and contain the x, y, and z coordinates of a vertex.
		 * Face lines start with "f" and contain the indices of the vertices that form a face.
		 * The indices are adjusted to be 0-based.
		 * The vertex and face information is temporarily stored in vectors and then transferred to the nods_ and cell_ matrices.
		 *
		 * @param[in] file   The file to load the model points from.
		 * @return          0 on success, or the line number if the file cannot be opened.
		 */
		int load_model_points(const char* file);

		/**
		 * @brief Load the cage from a file.
		 *
		 * This function loads the cage from a file and populates the cage_cell_ and cage_nods_ matrices.
		 * It reads the number of cage elements and vertices, followed by the vertex coordinates.
		 * The cage_cell_ matrix is resized to have 3 rows and the number of cage elements columns,
		 * and the cage_nods_ matrix is resized to have 3 rows and the number of cage nodes columns.
		 * The coordinates of the cage elements and nodes are stored in the respective matrices.
		 * The rest segment element lengths are calculated and stored in the cage_rest_uv_ matrix.
		 * The initial outward normal vectors of the cage elements are calculated and stored in the cage_normal_ matrix.
		 *
		 * @param[in] file   The file to load the cage from.
		 * @return          0 on success, or the line number if the file cannot be opened.
		 */

		int load_cage(const char* file);

		/**
		 * @brief Calculate the Green Coordinates for each point in the mesh.
		 *
		 * This function computes the Green Coordinates for each point in the mesh based on the provided cage and normal information.
		 * The Green Coordinates represent the deformation influence of the cage on each point.
		 *
		 * @return 0 on success.
		 */
		int calculate_green_coords();

		/**
		 * @brief Move a cage node by applying a displacement vector.
		 *
		 * This function moves a cage node by applying a displacement vector to its position.
		 *
		 * @param id    The index of the cage node to be moved.
		 * @param dx    Pointer to the displacement vector.
		 * @param disp  Flag indicating whether to apply the displacement incrementally (true) or directly (false).
		 * @return      0 on success.
		 */
		int move_cage(const size_t id, const double* dx, bool disp);

		/**
		* @brief Deform the mesh based on Green Coordinates and stretch ratios.
		*
		* This function deforms the mesh based on the computed Green Coordinates and stretch ratios.
		* It calculates the deformed positions of the mesh nodes by applying the transformations using the cage information.
		*
		* @return  0 on success.
		*/

		int deform();

		/**
		 * @brief Save the model to a file.
		 *
		 * This function saves the model to a file in OBJ format.
		 * It creates a 3D version of the nods_ matrix by setting the y-coordinate to 0.
		 * The 3D vertex coordinates are then written to the file, followed by the face information.
		 * The faces are written as 1-based indices.
		 *
		 * @param[in] file   The file to save the model to.
		 * @return          0 on success, or the line number if the file cannot be opened for writing.
		 */
		int saveToFile(const char* file);

		/**
		 * @brief Save the cage with normals to a file.
		 *
		 * This function saves the cage with normals to a file in OBJ format.
		 * It creates a new set of vertices for the normals by offsetting the original cage vertices.
		 * The vertex coordinates of the new cage with normals are written to the file, followed by the face information.
		 * The faces are written as 1-based indices.
		 *
		 * @param[in] file   The file to save the cage with normals to.
		 * @return          0 on success, or the line number if the file cannot be opened for writing.
		 */
		int saveToFile_normal(const char* file);

		/**
		 * @brief Save the cage to a file.
		 *
		 * This function saves the cage to a file in OBJ format.
		 * It opens the file for writing and checks if it was opened successfully.
		 * The vertex coordinates of the cage are written to the file, followed by the face information.
		 * The faces are written as 1-based indices.
		 *
		 * @param[in] file   The file to save the cage to.
		 * @return          0 on success, or the line number if the file cannot be opened for writing.
		 */

		int saveToFile_cage(const char* file);

	protected:
		/**
		 * @brief Calculate the outward normal vectors for each cage element.
		 *
		 * This function calculates the outward normal vectors for each cage element based on the vertices.
		 * It computes the face normal using the cross product and normalizes it.
		 * The calculated outward normal vectors are stored in the cage_normal_ matrix.
		 *
		 * @return  0 on success.
		 */
		int calculate_outward_normal();


		/**
		 * @brief Calculate the stretch ratios for each cage face.
		 *
		 * This function calculates the stretch ratios for each cage face based on the current and rest positions of the cage edges.
		 * It iterates over each face and calculates the dot products and stretch ratio using the formula:
		 *     sqrt(dotUcapUcap * dotV - 2 * dotUcapVcap * dotUcapV + dotVcapVcap * dotU) / (sqrt(8.0) * rest_area_(i))
		 * The calculated stretch ratios are stored in the s vector.
		 *
		 * @param[out] s   Vector to store the stretch ratios.
		 * @return        0 on success.
		 */
		int calculate_stretch_ratio(matd_t& s);
		//double calculateFaceArea(const matd_t& vt);

		 /**
		  * @brief Calculate the edge vectors for each cage face.
		  *
		  * This function calculates the edge vectors for each cage face based on the vertex positions.
		  * It iterates over each face, calculates the edge vectors, and stores them in the uv matrix.
		  *
          * @param[out] uv  Matrix to store the edge vectors.
		  * @return        0 on success.
		  */

		int calculate_cage_edge(matd_t& uv);


		/**
		 * @brief Calculate the Green Coordinates integral for a triangular element.
		 *
		 * This function computes the Green Coordinates integral for a given point `p` with respect to a triangular element defined by vertices `v1` and `v2`.
		 * The Green Coordinates quantify the influence of the triangular element on the deformation of the point `p`.
		 *
		 * @param p     The point of interest.
		 * @param v1    The first vertex of the triangular element.
		 * @param v2    The second vertex of the triangular element.
		 * @param eta   The reference point used in the integral calculation.
		 * @return      The Green Coordinates integral value.
		 */
		double GCTriInt(const matd_t& p, const matd_t& v1,
			const matd_t& v2, const matd_t& eta);

	private:

		mati_t cell_;
		matd_t nods_;

		mati_t cage_cell_;
		matd_t cage_nods_;
		matd_t cage_normal_;
		matd_t cage_rest_uv_;
		matd_t cage_curr_uv_;

		matd_t rest_area_;
		matd_t phi_;
		matd_t psi_;
	};


}


#endif 
