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


scalar_type change_pi(const scalar_type, const scalar_type);

scalar_type change_pi_dot(const scalar_type, const scalar_type);

scalar_type change_pi_dot_tau(const scalar_type, const scalar_type);

scalar_type integrate_trapeze(const plain_vector, const scalar_type);

scalar_type solution_initiale(const base_node);
scalar_type solution_initiale_velocity(const base_node);

scalar_type Hamiltonian(getfem::mesh &mesh, const getfem::mesh_fem mf, getfem::mesh_im &mim, 
                        const plain_vector y0, const plain_vector y1, const plain_vector u, const plain_vector p0, const plain_vector p1);

scalar_type Super_norm(const std::vector<plain_vector> Ju, const scalar_type Jtau);

void matlab_export(const scalar_type tau, const scalar_type tol);

void matlab_export_solution(const std::vector<plain_vector> y);

void matlab_export_solution_velocity(const std::vector<plain_vector> y);

void matlab_export_pressure(const std::vector<plain_vector> pcond);

void matlab_export_control(const std::vector<plain_vector> u);

void matlab_export_time(const plain_vector time_vec);

void asm_navier_stokes_nl_tgm(sparse_matrix &TGM, getfem::mesh_im &mim,
                              const getfem::mesh_fem mf, const plain_vector U,
                              const getfem::mesh_region &rg);
























