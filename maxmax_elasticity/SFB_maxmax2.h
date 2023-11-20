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
//#include "getfem/getfem_Coulomb_friction.h"
#include "getfem/getfem_import.h"
//#include "getfem/getfem_inter_element.h"
#include "gmm/gmm.h"
#include "gmm/gmm_def.h"

	//#include "advanced_inversion.h"
	//#include "advanced_inversion_stabilization.h"

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


#include "SFB_maxmax_routines2.h"
#include "SFB_maxmax_solving2.h"





