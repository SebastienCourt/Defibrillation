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

#include "SFB_maxmax2.h"




extern scalar_type Tmax;
extern size_type N_time;
extern double alphax;
extern double betax;
extern double lambdax;
extern double mux;
extern double nux;
extern double thetax;



scalar_type change_pi(const scalar_type tau, const scalar_type s) {

  scalar_type valux = 0.0;
  if (s < 0.5) {
	valux = 2.0*s*tau;
  }
  else {
	valux = 2.0*(Tmax-tau)*s + 2.0*tau - Tmax;
  }

  return valux;

}


scalar_type change_pi_dot(const scalar_type tau, const scalar_type s) {

  scalar_type valux = 0.0;
  if (s < 0.5) {
	valux = 2.0*tau;
  }
  else {
	valux = 2.0*(Tmax-tau);
  }

  return valux;

}


scalar_type change_pi_dot_tau(const scalar_type tau, const scalar_type s) {

  scalar_type valux = 0.0;
  if (s < 0.5) {
	valux = 2.0;
  }
  else {
	valux = -2.0;
  }

  return valux;

}




scalar_type integrate_trapeze(const plain_vector y, const scalar_type ds) {

  scalar_type integrale = 0.0;

  for (size_type i = 0; i < y.size()-1; ++i) { 

	integrale += 0.5*ds*(y[i]+y[i+1]);

  }

  return integrale;

} 


scalar_type solution_initiale(const base_node p) {

  scalar_type y0 = 0.0;

  y0 = 0.0;//sin(2.0*M_PI*p[0]);

  //y0 = -3.0/(4.0*M_PI)*sin((p[0]/3.0)*M_PI*(p[0]+3.0));//p[0];//10.0*(1.0 - exp(-(p[0])))*(exp(-(p[0])) - exp(-1.0));
  	//y0 = sin(M_PI*p[0]);

  /*if (p[0] < 0.333333333333333333333333333) {
	y0 = sin(M_PI*3.0*p[0]);
  }*/

  return y0;


}


scalar_type solution_initiale_velocity(const base_node p) {

  scalar_type y1 = 0.0;

  y1 = sin(0.5*M_PI*p[0]);

  //y0 = -3.0/(4.0*M_PI)*sin((p[0]/3.0)*M_PI*(p[0]+3.0));//p[0];//10.0*(1.0 - exp(-(p[0])))*(exp(-(p[0])) - exp(-1.0));
  	//y0 = sin(M_PI*p[0]);

  /*if (p[0] < 0.333333333333333333333333333) {
	y0 = sin(M_PI*3.0*p[0]);
  }*/

  return y1;


}


scalar_type Hamiltonian(getfem::mesh &mesh, const getfem::mesh_fem mf, getfem::mesh_im &mim, 
			 const plain_vector y0, const plain_vector y1, const plain_vector u, const plain_vector p0, const plain_vector p1) {


  scalar_type H = -0.5*alphax*gmm::vect_sp(u, u);

  size_type nb_dof = mf.nb_dof();
  getfem::model MS;
  MS.add_fem_variable("w", mf);
  MS.add_initialized_scalar_data("betaxxx", betax);
  MS.add_initialized_scalar_data("nuxxx", nux);
  MS.add_initialized_scalar_data("muxxx", mux);
  MS.add_initialized_scalar_data("lambdaxxx", lambdax);
  MS.add_initialized_fem_data("y0", mf, y0);
  MS.add_initialized_fem_data("y1", mf, y1);

  //getfem::add_generic_elliptic_brick(MS, mim, "w", "nuxxx");
  getfem::add_linear_generic_assembly_brick(MS, mim, "nuxxx*Grad_y1:Grad_Test_w", -1);
  getfem::add_linear_generic_assembly_brick(MS, mim, "betaxxx*((Id(meshdim)+Grad_y0)*( muxxx*(Grad_y0'*Grad_y0+Grad_y0+Grad_y0') )):Grad_Test_w", -1);
  getfem::add_linear_generic_assembly_brick(MS, mim, "betaxxx*((Id(meshdim)+Grad_y0)*( 0.5*lambdaxxx*Trace(Grad_y0'*Grad_y0+Grad_y0+Grad_y0')*Id(qdim(w)) )):Grad_Test_w", -1);

  	// Dirichlet condition with simplification
  		getfem::mesh_region r;
  		getfem::outer_faces_of_mesh(mesh,r);
  		int DIRICHLET_BOUNDARY_NUM = 11;
  		for (getfem::mr_visitor i(r); !i.finished(); ++i) {
  
  	                base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
	                un /= gmm::vect_norm2(un);
	                if ( gmm::abs(un[0]+1.0) < 1e-12 ) {
		                mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(),i.f());
	                }
                        //mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(),i.f());
  		}
  	plain_vector HDC(nb_dof);
  	gmm::clear(HDC);
  	MS.add_initialized_fem_data("Dirichletdata", mf, HDC);	
  	add_Dirichlet_condition_with_simplification(MS, "w", DIRICHLET_BOUNDARY_NUM, "Dirichletdata");

  MS.assembly(getfem::model::BUILD_ALL);

  //sparse_matrix A(nb_dof, nb_dof);
  //gmm::copy(MS.real_tangent_matrix(), A);
  //gmm::scale(A, -1.0);

  sparse_matrix M(nb_dof, nb_dof);
  getfem::asm_mass_matrix(M, mim, mf);

  plain_vector F(nb_dof);
  gmm::copy(MS.real_rhs(), F);

  H += gmm::vect_sp(M, p1, F);

  H -= gmm::vect_sp(M, p0,y1);

  plain_vector uu(nb_dof);
  gmm::copy(u, uu);
  for (size_type k=0; k<nb_dof; ++k) {
	base_node bn = mf.point_of_basic_dof(k);
	if (bn[0] < 0.75) {
		uu[k] = 0.0;
	}
  }

  H -= gmm::vect_sp(M, p1, uu);


  return H;

}




scalar_type Super_norm(const std::vector<plain_vector> Ju, const scalar_type Jtau) {

  scalar_type SN = Jtau*Jtau;

  for (size_type i=0; i< N_time; ++i) {

	SN += gmm::vect_sp(Ju[i], Ju[i]);

  }

  return gmm::sqrt(SN);

}


/*********************************************************************/
/*********************** Export for matlab ***************************/
/*********************************************************************/



void matlab_export_time(const plain_vector time_vec) {


  std::ofstream timevex;
  timevex.open("./MATLAB_SFB/time_vec2.m", ios::out|ios::trunc);

  timevex << "function t = time_vec" << endl;
  timevex << endl;
  timevex << "t = [ " ;

  for (size_type i=0; i<time_vec.size()-1; ++i) {

	timevex << time_vec[i] << " ; " ;
	
  }

  timevex << time_vec[time_vec.size()-1] << " ]; " << endl;
  timevex.close();


}


void matlab_export_pressure(const std::vector<plain_vector> pcond) {

  size_type nb_dof = pcond[0].size();
  plain_vector yt(nb_dof);

  	
  std::ofstream solexp;
  solexp.open("./MATLAB_SFB/Pressure.m", ios::out|ios::trunc);

  solexp << "function p = pressure" << endl;
  solexp << endl;
  solexp << "p = [ " ;

  for (size_type i=0; i<N_time-1; ++i) {

	gmm::clear(yt);
	gmm::copy(pcond[i], yt);
	for (size_type k=0; k<nb_dof-1; ++k) {

		solexp << yt[k] << ", " ;

	}
	solexp << yt[nb_dof-1] << " ; " << endl;
	
  }

  gmm::clear(yt);
  gmm::copy(pcond[N_time-1], yt);
  for (size_type k=0; k<nb_dof-1; ++k) {

		solexp << yt[k] << ", " ;

  }
  solexp << yt[nb_dof-1] << " ]; " << endl;

  solexp.close();


}



void matlab_export_solution(const std::vector<plain_vector> y) {

  size_type nb_dof = y[0].size();
  plain_vector yt(nb_dof);

  	
  std::ofstream solexp;
  solexp.open("./MATLAB_SFB/State.m", ios::out|ios::trunc);

  solexp << "function y = state" << endl;
  solexp << endl;
  solexp << "y = [ " ;

  for (size_type i=0; i<N_time-1; ++i) {

	gmm::clear(yt);
	gmm::copy(y[i], yt);
	for (size_type k=0; k<nb_dof-1; ++k) {

		solexp << yt[k] << ", " ;

	}
	solexp << yt[nb_dof-1] << " ; " << endl;
	
  }

  gmm::clear(yt);
  gmm::copy(y[N_time-1], yt);
  for (size_type k=0; k<nb_dof-1; ++k) {

		solexp << yt[k] << ", " ;

  }
  solexp << yt[nb_dof-1] << " ]; " << endl;

  solexp.close();


}



void matlab_export_solution_velocity(const std::vector<plain_vector> y) {

  size_type nb_dof = y[0].size();
  plain_vector yt(nb_dof);

  	
  std::ofstream solexp;
  solexp.open("./MATLAB_SFB/State_velocity.m", ios::out|ios::trunc);

  solexp << "function y = state_velocity" << endl;
  solexp << endl;
  solexp << "y = [ " ;

  for (size_type i=0; i<N_time-1; ++i) {

	gmm::clear(yt);
	gmm::copy(y[i], yt);
	for (size_type k=0; k<nb_dof-1; ++k) {

		solexp << yt[k] << ", " ;

	}
	solexp << yt[nb_dof-1] << " ; " << endl;
	
  }

  gmm::clear(yt);
  gmm::copy(y[N_time-1], yt);
  for (size_type k=0; k<nb_dof-1; ++k) {

		solexp << yt[k] << ", " ;

  }
  solexp << yt[nb_dof-1] << " ]; " << endl;

  solexp.close();


}



void matlab_export_control(const std::vector<plain_vector> y) {

  size_type nb_dof = y[0].size();
  plain_vector yt(nb_dof);

  	
  std::ofstream solexp;
  solexp.open("./MATLAB_SFB/Control.m", ios::out|ios::trunc);

  solexp << "function u = control" << endl;
  solexp << endl;
  solexp << "u = [ " ;

  for (size_type i=0; i<N_time-1; ++i) {

	gmm::clear(yt);
	gmm::copy(y[i], yt);
	for (size_type k=0; k<nb_dof-1; ++k) {

		solexp << yt[k] << ", " ;

	}
	solexp << yt[nb_dof-1] << " ; " << endl;
	
  }

  gmm::clear(yt);
  gmm::copy(y[N_time-1], yt);
  for (size_type k=0; k<nb_dof-1; ++k) {

		solexp << yt[k] << ", " ;

  }
  solexp << yt[nb_dof-1] << " ]; " << endl;

  solexp.close();


}




void matlab_export(const scalar_type tau, const scalar_type tol) {


  std::ofstream evoltau;
  evoltau.open("./MATLAB_SFB/evol_tau2.txt", ios::out|ios::app);
  evoltau << tau << ",   " ;
  evoltau.close();

  std::ofstream evoltol;
  evoltol.open("./MATLAB_SFB/evol_tol2.txt", ios::out|ios::app);
  evoltol << tol << ",   " ;
  evoltol.close();


}


	//template<typename MAT, typename VEC>
void asm_navier_stokes_nl_tgm(sparse_matrix &TGM, getfem::mesh_im &mim,
                              const getfem::mesh_fem mf, const plain_vector U,
                              const getfem::mesh_region &rg=getfem::mesh_region::all_convexes()) {

    getfem::generic_assembly assem;
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_data(U);
    assem.push_mat(TGM);
    assem.set("u=data(#1);"
              "t=comp(vBase(#1).vGrad(#1).vBase(#1));"
              "M(#1, #1) +=u(j).t(:,k,j,k,l,:,l)+u(i).t(:,k,:,k,l,i,l);");
    assem.assembly(rg);

}









