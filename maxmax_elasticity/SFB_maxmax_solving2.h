#include <iostream>
#include <fstream>

#include "getfem/getfem_assembling.h" /* import assembly methods      */
#include "getfem/getfem_export.h"     /* export functions             */
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_mesh_fem_level_set.h"
#include "getfem/getfem_partial_mesh_fem.h"
#include "getfem/getfem_import.h"
//#include "getfem/getfem_inter_element.h"
#include "gmm/gmm.h"
#include "gmm/gmm_def.h"



using namespace std;
using bgeot::base_small_vector;
using bgeot::base_vector;
using bgeot::base_node;
using bgeot::scalar_type; /* = double */
using bgeot::short_type;  /* = short */
using bgeot::size_type;   /* = unsigned long */
using bgeot::dim_type;
using bgeot::base_matrix; /* small dense matrix. */


typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;
typedef gmm::row_matrix<sparse_vector> sparse_row_matrix;



extern scalar_type Tmax;
extern size_type N_time;
extern double alphax;
extern double betax;
extern double lambdax;
extern double mux;
extern double nux;
extern double thetax;


void compute_gradients(std::vector<plain_vector> &Ju, scalar_type &Jtau,
                       getfem::mesh &mesh, const getfem::mesh_fem mf, getfem::mesh_im &mim, 
                       const std::vector<plain_vector> y0, const std::vector<plain_vector> y1,
                       const std::vector<plain_vector> u, const scalar_type tau, 
                       const std::vector<plain_vector> p0, const std::vector<plain_vector> p1);


void solve_direct(std::vector<plain_vector> &y0, std::vector<plain_vector> &y1, std::vector<plain_vector> &pcond, getfem::mesh &mesh, const getfem::mesh_fem mf, getfem::mesh_im &mim, 
                  const plain_vector y00, const plain_vector y10,
                  const std::vector<plain_vector> u, const scalar_type tau);


void solve_adjoint(std::vector<plain_vector> &p0, std::vector<plain_vector> &p1, getfem::mesh &mesh, const getfem::mesh_fem mf, getfem::mesh_im &mim, 
		   const std::vector<plain_vector> y0, const std::vector<plain_vector> y1, const scalar_type tau);



















