#include <iostream>
#include <fstream>

#include "getfem/getfem_assembling.h" /* import assembly methods      */
#include "getfem/getfem_export.h"     /* export functions             */
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_mesh_fem_level_set.h"
#include "getfem/getfem_nonlinear_elasticity.h"
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


#include "SFB_maxmax2.h"

scalar_type Tmax;
size_type N_time;
double alphax;
double betax;
double lambdax;
double mux;
double nux;
double thetax;


int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.


  // Load data
  bgeot::md_param PARAM;
  PARAM.read_command_line(argc, argv);
  betax = PARAM.real_value("BETAX", "Density of the material");
  lambdax = PARAM.real_value("LAMBDAX", "Elastic coefficient");
  mux = PARAM.real_value("MUX", "Elasticity Coefficient");
  alphax = PARAM.real_value("ALPHAX", "Cost parameter");
  nux = PARAM.real_value("NUX", "Viscosity");
  Tmax = PARAM.real_value("TMAX", "Final time");
  N_time = int(PARAM.int_value("N_TIME", "Number of time steps"));
  thetax = PARAM.real_value("THETA", "Parameter for the theta-method");

  int Nmax = int(PARAM.int_value("NMAX", "Maximum number of iterations"));
  double eps = PARAM.real_value("EPS", "Accuracy");
  double tau0 = PARAM.real_value("TAU0", "Initial optimal time");

  // The mesh
  getfem::mesh mesh;
  std::string MESH_FILE = PARAM.string_value("MESH_FILE", "Mesh file");
  getfem::import_mesh(MESH_FILE, mesh);
  	//dim_type N = mesh.dim();
  scalar_type h = mesh.minimal_convex_radius_estimate();
  cout << "h = " << h << endl;

  // The mesh_fem
  getfem::mesh_fem mf(mesh);
  std::string FEM = PARAM.string_value("FEM", "Finite element method for the displacement");
  	//mf.set_qdim(N);
  mf.set_finite_element(mesh.convex_index(), getfem::fem_descriptor(FEM));
  size_type nb_dof = mf.nb_dof();
  cout << "nb_dof = " << nb_dof << endl;

  getfem::mesh_im mim(mesh);
  std::string IM = PARAM.string_value("IM", "Mesh file");
  mim.set_integration_method(mesh.convex_index(), getfem::int_method_descriptor(IM));



  // Initialize
  std::vector< plain_vector > y0(N_time);
  std::vector< plain_vector > y1(N_time);
  std::vector< plain_vector > pcond(N_time);
  std::vector< plain_vector > u(N_time);
  std::vector< plain_vector > p0(N_time);
  std::vector< plain_vector > p1(N_time);
  scalar_type tau = tau0;//Tmax*0.5;

  std::vector< plain_vector > Ju(N_time);
  scalar_type Jtau = 0.0;

  for (size_type i=0; i<N_time; ++i) { 
  	gmm::resize(y0[i], nb_dof);
  	gmm::resize(y1[i], nb_dof);
  	gmm::resize(pcond[i], 1);
  	gmm::resize(u[i], nb_dof);
  	gmm::resize(p0[i], nb_dof);
  	gmm::resize(p1[i], nb_dof);
  	gmm::resize(Ju[i], nb_dof);
  }

  plain_vector y00(nb_dof);
  getfem::interpolation_function(mf, y00, solution_initiale);
  plain_vector y10(nb_dof);
  getfem::interpolation_function(mf, y10, solution_initiale_velocity);
  //gmm::clear(y10);


  // Intitialize control
  for (size_type i=0; i<N_time; ++i) {

        gmm::clear(u[i]);

  }

  solve_direct(y0, y1, pcond, mesh, mf, mim, y00, y10, u, tau);
  solve_adjoint(p0, p1, mesh, mf, mim, y0, y1, tau);
  compute_gradients(Ju, Jtau, mesh, mf, mim, y0, y1, u, tau, p0, p1);

  std::vector< plain_vector > u_mem(N_time);
  std::vector< plain_vector > Ju_mem(N_time);
  for (size_type i=0; i<N_time; ++i) { 
	gmm::resize(u_mem[i], nb_dof);
	gmm::copy(u[i], u_mem[i]);
	gmm::resize(Ju_mem[i], nb_dof);
	gmm::copy(Ju[i], Ju_mem[i]);
  }
  scalar_type tau_mem = tau;
  scalar_type Jtau_mem = Jtau;

  scalar_type tol = Super_norm(Ju, Jtau);


  // First export (initialize)
  std::ofstream evoltau;
  evoltau.open("./MATLAB_SFB/evol_tau2.txt", ios::out|ios::trunc);
  evoltau.close();
  std::ofstream evoltol;
  evoltol.open("./MATLAB_SFB/evol_tol2.txt", ios::out|ios::trunc);
  evoltol.close();

  matlab_export(tau, tol);







  // Gradient-steps: Armijo rule style
  bool test = (tol < eps);
  int compteur = 0;

  std::vector< plain_vector > uu(N_time);
  std::vector< plain_vector > Juu(N_time);
  for (size_type i=0; i<N_time; ++i) {
	gmm::resize(uu[i], nb_dof); 
	gmm::copy(u[i], uu[i]);
	gmm::resize(Juu[i], nb_dof);
	gmm::copy(Ju[i], Juu[i]);
  }
  scalar_type ttau = tau;
  scalar_type Jttau = Jtau;
  double amijo_coeff = 1.0;

  for (size_type k=0; k<1; ++k) {

	while ((!test) && (compteur<1)) {

		compteur++;
		amijo_coeff *= 0.1;

		// The unknowns
		for (size_type i=0; i<N_time; ++i) {
			gmm::add(gmm::scaled(Juu[i], amijo_coeff), uu[i]);
		}
		ttau += amijo_coeff*Jttau;

		// The states
		solve_direct(y0, y1, pcond, mesh, mf, mim, y00, y10, uu, ttau);
		solve_adjoint(p0, p1, mesh, mf, mim, y0, y1, ttau);

		// The gradients
		compute_gradients(Juu, Jttau, mesh, mf, mim, y0, y1, uu, ttau, p0, p1);            
		tol = Super_norm(Juu, Jttau);
		test = (tol < eps);
            
	}

	
	for (size_type i=0; i<N_time; ++i) { 
		gmm::copy(uu[i], u[i]);
		gmm::copy(Juu[i], Ju[i]);

	}
	tau = ttau;
	Jtau = Jttau;
	
	cout << "Compteur amijo = " << compteur << endl;
	cout << "Tol Amijo = " << tol << endl;

  }



matlab_export(tau, tol);




/***************************************************************************************/
/********************************* Barzilai - Borwein **********************************/
/***************************************************************************************/

scalar_type alpha_u1 = 0.0;
scalar_type alpha_u2 = 0.0;
scalar_type alpha_tau1 = 0.0;
scalar_type alpha_tau2 = 0.0;

std::vector< plain_vector > Su(N_time);
std::vector< plain_vector > Gu(N_time);
for (size_type i=0; i<N_time; ++i) {
	gmm::resize(Su[i], nb_dof); 
	gmm::resize(Gu[i], nb_dof);
}
scalar_type Gtau = 0.0;
scalar_type Stau = 0.0;

scalar_type au1_up, au1_down, au2_up, au2_down;

compteur = 0;


while ((tol > eps) && (compteur < Nmax)) {

	compteur++;

	// Compute scalar-products
	for (size_type i=0; i<N_time; ++i) {
		gmm::clear(Su[i]);
		gmm::add(u[i], gmm::scaled(u_mem[i], -1.0), Su[i]);
		gmm::clear(Gu[i]);
		gmm::add(Ju[i], gmm::scaled(Ju_mem[i], -1.0), Gu[i]);
	}
	Stau = tau - tau_mem;
	Gtau = Jtau - Jtau_mem;
    
	au1_up = 0.0; 
	au1_down = 0.0;
	au2_up = 0.0; 
	au2_down = 0.0;
	for (size_type i=0; i<N_time; ++i) {
		au1_up += gmm::vect_sp(Su[i], Su[i]);
		au1_down += gmm::vect_sp(Su[i], Gu[i]);
		au2_up += gmm::vect_sp(Su[i], Gu[i]);
		au2_down += gmm::vect_sp(Gu[i], Gu[i]);
	}

	cout << "au1_down = " << au1_down << endl;
	cout << "au2_down = " << au2_down << endl;
	cout << "Stau = " << Stau << endl;
	cout << "Gtau = " << Gtau << endl;

	if (gmm::abs(au1_down) > 0) alpha_u1 = au1_up/au1_down;
	if (gmm::abs(au2_down) > 0) alpha_u2 = au2_up/au2_down;
    
	if ( (gmm::abs(Stau) > 0) && (gmm::abs(Gtau) > 0) ) {
		alpha_tau1 = (Stau*Stau)/(Stau*Gtau);
		alpha_tau2 = (Stau*Gtau)/(Gtau*Gtau);
	}


	// Save the previous quantities
	for (size_type i=0; i<N_time; ++i) {
		gmm::copy(u[i], u_mem[i]);
		gmm::copy(Ju[i], Ju_mem[i]);
	}
	tau_mem = tau;
	Jtau_mem = Jtau;


	// Update of the unknowns
	//if (compteur & 1) {
		for (size_type i=0; i<N_time; ++i) {
			gmm::add(gmm::scaled(Ju[i], -alpha_u2), u[i]);
		}
		tau = tau - alpha_tau2*Jtau;
	//}
    	/*else {
		for (size_type i=0; i<N_time; ++i) {
			gmm::add(gmm::scaled(Ju[i], -alpha_u2), u[i]);
		}
		tau = tau - alpha_tau2*Jtau;
	}*/


	/* Update the control
	for (size_type i=0; i<N_time; ++i) { 
		gmm::add(u[i], du[i], u[i]);
	}
	tau += dtau;*/

	// Update the solution and the adjoint-state
cout << "OUOUOUOUOUUOOUOUUOUOOUOUOUOUOUOUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU" << endl;
	solve_direct(y0, y1, pcond, mesh, mf, mim, y00, y10, u, tau);
cout << "OUOUOUOUOUUOOUOUUOUOOUOUOUOUOUOUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU 222222222222" << endl;
	solve_adjoint(p0, p1, mesh, mf, mim, y0, y1, tau);

	// Update the gradient and the tolerance
	compute_gradients(Ju, Jtau, mesh, mf, mim, y0, y1, u, tau, p0, p1);
	tol = Super_norm(Ju, Jtau);
	cout << "Tol Barzilai-Borwein = " << tol << endl;
	matlab_export(tau, tol);

}


/****************************************************************************************/
/************************************* Export *******************************************/
/****************************************************************************************/


  // Post-processing

  plain_vector time_vec(N_time);
  double step_ds = 1.0/double(N_time);


  for (size_type i=0; i<N_time; ++i) { 

	time_vec[i] = change_pi(tau, step_ds*double(i+1));

  }


  // Export the result - vtk format

  cout << "Exporting..." << endl;

  for (int i=0; i<int(N_time); ++i) { 

	{
  	std::string str_y = "y_solution";
  		//char sy[128]; sprintf(sy, "VTK_SFB/%s%d.vtk", str_y.c_str(), time_vec[i]);
	char sy[128]; sprintf(sy, "VTK_SFB/%s%d.vtk", str_y.c_str(), i);
  	getfem::vtk_export exp(sy, (2==1));
  	exp.exporting(mf); 
  	exp.write_point_data(mf, y0[i], "solution");
	}

	{
  	std::string str_u = "u_control";
  		//char su[128]; sprintf(su, "VTK_SFB/%s%d.vtk", str_u.c_str(), time_vec[i]);
	char su[128]; sprintf(su, "VTK_SFB/%s%d.vtk", str_u.c_str(), i);
  	getfem::vtk_export exp2(su, (2==1));
  	exp2.exporting(mf); 
  	exp2.write_point_data(mf, u[i], "control");
	}

  }


  cout << "Export of the solution in vtk format done." << endl;


  // Export for matlab vizu
  matlab_export_solution(y0);
  matlab_export_solution_velocity(y1);
  matlab_export_pressure(pcond);
  matlab_export_control(u);
  matlab_export_time(time_vec);





  return 0; 

}




























