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

extern scalar_type Tmax;
extern size_type N_time;
extern double alphax;
extern double betax;
extern double lambdax;
extern double mux;
extern double nux;
extern double thetax;





/********************************************************************************************************************************/
/***************************************** Compute the gradients ****************************************************************/
/********************************************************************************************************************************/

void compute_gradients(std::vector<plain_vector> &Ju, scalar_type &Jtau,
                       getfem::mesh &mesh, const getfem::mesh_fem mf, getfem::mesh_im &mim,
                       const std::vector<plain_vector> y0, const std::vector<plain_vector> y1, 
                       const std::vector<plain_vector> u, const scalar_type tau, 
                       const std::vector<plain_vector> p0, const std::vector<plain_vector> p1) {


  size_type nb_dof = mf.nb_dof();
  scalar_type ds = 1.0/double(N_time);

  // Compute Jtau
  plain_vector Hami(N_time);
  for (size_type i=0; i<N_time; ++i) { 
        Hami[i] = change_pi_dot_tau(tau,ds*double(i)) * Hamiltonian(mesh, mf, mim, y0[i], y1[i], u[i], p0[i], p1[i]);
  }

  Jtau = integrate_trapeze(Hami, ds);


  // Compute Ju
  for (size_type i=0; i<N_time; ++i) { 
        gmm::copy(gmm::scaled(u[i], -alphax), Ju[i]);
  }
  std::vector<plain_vector> pp(N_time);
  for (size_type i=0; i<N_time; ++i) { 
        gmm::resize(pp[i], nb_dof);
  }


  plain_vector ppk(nb_dof);
  for (size_type i=0; i<N_time; ++i) { 

        gmm::clear(ppk);
        gmm::copy(p1[i], ppk);
        gmm::scale(ppk, -1.0);
        for (size_type k=0; k<nb_dof; ++k) {
                base_node bn = mf.point_of_basic_dof(k);
                if (bn[0] < 0.75) {
                        ppk[k] = 0.0;
                }
        }
        gmm::copy(ppk, pp[i]);

  }

  for (size_type i=0; i<N_time; ++i) { 
        gmm::add(pp[i], Ju[i]);
  }

  for (size_type i=0; i<N_time; ++i) { 
        gmm::scale(Ju[i], change_pi_dot(tau, ds*double(i)));
  }


}





/********************************************************************************************************************************/
/***************************************** Direct solve *************************************************************************/
/********************************************************************************************************************************/


void solve_direct(std::vector<plain_vector> &y0, std::vector<plain_vector> &y1, std::vector<plain_vector> &pcond,
                  getfem::mesh &mesh, const getfem::mesh_fem mf, getfem::mesh_im &mim, 
                  const plain_vector y00, const plain_vector y10, const std::vector<plain_vector> u, const scalar_type tau) {


  size_type nb_dof = mf.nb_dof();

  double step_dt = 1.0/double(N_time);
  double timex = 0.0;

  // Initial condition
  for (size_type i=0; i<N_time; ++i) { 
          gmm::resize(y0[i], nb_dof);
          gmm::resize(y1[i], nb_dof);
  }

  gmm::copy(y00, y0[0]);
  gmm::copy(y10, y1[0]);

  //  Region for the Dirichlet condition
  getfem::mesh_region r;
  getfem::outer_faces_of_mesh(mesh,r);
  int DIRICHLET_BOUNDARY_NUM = 11;
  int DIRICHLET_BOUNDARY_NUM_VELOCITY = 21;
  for (getfem::mr_visitor i(r); !i.finished(); ++i) {
  
        //mesh.region(DIRICHLET_BOUNDARY_NUM_VELOCITY).add(i.cv(),i.f());
        base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
        un /= gmm::vect_norm2(un);
        if ( gmm::abs(un[0]+1.0) < 1e-12 ) {
                mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(),i.f());
                mesh.region(DIRICHLET_BOUNDARY_NUM_VELOCITY).add(i.cv(),i.f());
        }
        //mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(),i.f());
  }


  getfem::model MS;
  MS.add_fem_variable("y", mf);

  plain_vector Pixdot(nb_dof);
  plain_vector NuxPixdot(nb_dof);
  for (size_type k=0; k<nb_dof; ++k) {
        Pixdot[k] = change_pi_dot(tau,timex);
        NuxPixdot[k] = nux*change_pi_dot(tau,timex);
  }
  MS.add_initialized_fem_data("pixdot", mf, Pixdot);
  MS.add_initialized_fem_data("nuxpixdot", mf, NuxPixdot);



  // The pressure
  MS.add_fixed_size_variable("p_cond", 1);
  plain_vector Aire(nb_dof);
  plain_vector Voly(nb_dof);
  for (size_type i = 0; i < nb_dof; ++i) {
       Voly[i] = 1.0; // That is det(I+nabla u_0) when u_0 = 0
       Aire[i] = 1.0;
  }
  //getfem::asm_source_term(Cv, mim, mf, mf, Aire);
  sparse_matrix M(nb_dof, nb_dof);
  getfem::asm_mass_matrix(M, mim, mf, mf);
  plain_vector C_rhs(1);
  C_rhs[0] = gmm::vect_sp(M, Voly, Aire);;
  /*getfem::add_explicit_matrix(MS, "p_cond", "p", gmm::row_vector(Cv));
  getfem::add_explicit_matrix(MS, "p", "p_cond", gmm::transposed(gmm::row_vector(Cv)));
  getfem::add_explicit_rhs(MS, "p_cond", C_rhs);*/
  //MS.add_initialized_scalar_data("Vol", C_rhs);
  MS.add_initialized_fem_data("Aire", mf, Aire);        
  getfem::add_nonlinear_generic_assembly_brick(MS, mim, "p_cond*((1+Grad_y)*Aire)*Test_p_cond");
  getfem::add_explicit_rhs(MS, "p_cond", C_rhs);


  // Nonlinear brick
  MS.add_initialized_scalar_data("betaxxx", betax);
  MS.add_initialized_scalar_data("lambdaxxx", lambdax);
  MS.add_initialized_scalar_data("muxxx", mux);
          //getfem::add_nonlinear_generic_assembly_brick(MS, mim, "-pixdot*betaxxx*(Grad_y.y)*Test_y", -1); // For Burgers
  //getfem::add_linear_generic_assembly_brick(MS, mim, "pixdot*betaxxx*(Grad_y)*Grad_Test_y", -1);
  MS.add_macro("mytensor", "0.5*(Grad_y'*Grad_y+Grad_y+Grad_y')");
  getfem::add_nonlinear_generic_assembly_brick(MS, mim, "pixdot*betaxxx*((Id(meshdim)+Grad_y)*( 2.0*muxxx*mytensor )):Grad_Test_y", -1);
  getfem::add_nonlinear_generic_assembly_brick(MS, mim, "pixdot*betaxxx*((Id(meshdim)+Grad_y)*( lambdaxxx*Trace(mytensor)*Id(qdim(y)) )):Grad_Test_y", -1);


  // Dirichlet condition with simplification
  plain_vector H(nb_dof);
  gmm::clear(H);
  MS.add_initialized_fem_data("Dirichletdata", mf, H);        
  add_Dirichlet_condition_with_simplification(MS, "y", DIRICHLET_BOUNDARY_NUM, "Dirichletdata");

  // Volumic source term: the control
  plain_vector Pix_u(nb_dof);
  gmm::copy(u[0], Pix_u);
  gmm::scale(Pix_u, change_pi_dot(tau,timex));
  MS.add_initialized_fem_data("Control_rhs", mf, Pix_u);
  getfem::add_source_term_brick(MS, mim, "y", "Control_rhs");


  // Transient part 
          //double thetax = 0.5;
          //add_theta_method_for_first_order(MS, "y", thetax); // For Burgers  
  add_theta_method_for_second_order(MS, "y", thetax);
  add_mass_brick(MS, mim, "Dot2_y");


  // The regularizing term
          //getfem::add_Laplacian_brick(MS, mim, "y");
  getfem::add_generic_elliptic_brick(MS, mim, "Dot_y", "nuxpixdot");
  add_Dirichlet_condition_with_simplification(MS, "Dot_y", DIRICHLET_BOUNDARY_NUM_VELOCITY, "Dirichletdata");


  MS.set_time(timex);
  MS.set_time_step(step_dt);


  plain_vector Y0(nb_dof), Y1(nb_dof), P0(1);
  gmm::copy(y00, Y0);
  gmm::copy(Y0, MS.set_real_variable("Previous_y"));
  gmm::copy(y10, Y1);
  gmm::copy(Y1, MS.set_real_variable("Previous_Dot_y"));
          //MS.perform_init_time_derivative(step_dt*0.05);

  gmm::iteration iter(1e-9, 0, 40000);
        //plain_vector Y(nb_dof);
 

  // Time iteration and solving
  for (size_type i=1; i<N_time; ++i) {

        timex += step_dt;

        // Update the data
        gmm::clear(Pixdot);
        for (size_type k=0; k<nb_dof; ++k) {
                Pixdot[k] = change_pi_dot(tau,timex);
        }
        gmm::copy(Pixdot, MS.set_real_variable("pixdot"));
        gmm::copy(gmm::scaled(Pixdot, nux), MS.set_real_variable("nuxpixdot"));
        gmm::clear(Pix_u);
        gmm::copy(gmm::scaled(u[i], change_pi_dot(tau,timex)), Pix_u);
        gmm::copy(Pix_u, MS.set_real_variable("Control_rhs"));
  //getfem::add_nonlinear_generic_assembly_brick(MS, mim, "-pixdot*betaxxx*(Grad_y.y)", -1);

        // Solving
                //cout << "Solving direct for time s = " << timex << " ..."  << endl;
        cout << "Standard solve -- " << i << endl;
        cout << "Changepidot = " << change_pi_dot(tau,timex) << endl;
        cout << "Timex = " << timex << endl;
        cout << "Tau = " << tau << endl;
        getfem::standard_solve(MS, iter);
                //cout << "Solved." << endl;

        // Extract the solution
        gmm::clear(Y0);
        gmm::clear(Y1);
        gmm::clear(P0);
        gmm::copy(MS.real_variable("y"), Y0);
        gmm::copy(MS.real_variable("Dot_y"), Y1);
        gmm::copy(MS.real_variable("p_cond"), P0);
        gmm::copy(Y0, y0[i]);
        gmm::copy(Y1, y1[i]);
        gmm::copy(P0, pcond[i]);

        MS.shift_variables_for_time_integration();

  }

} 





/********************************************************************************************************************************/
/***************************************** Adjoint solve ************************************************************************/
/********************************************************************************************************************************/



void solve_adjoint(std::vector<plain_vector> &p0, std::vector<plain_vector> &p1, 
                   getfem::mesh &mesh, const getfem::mesh_fem mf, getfem::mesh_im &mim, 
                   const std::vector<plain_vector> y0,  const std::vector<plain_vector> y1, const scalar_type tau) {


  size_type nb_dof = mf.nb_dof();

  double step_dt = 1.0/double(N_time);
  double timex = 1.0;


  // Initial condition
  gmm::resize(p0, N_time);
  gmm::resize(p1, N_time);
  for (size_type i=0; i<N_time; ++i) { 
          gmm::resize(p0[i], nb_dof);
          gmm::resize(p1[i], nb_dof);
  }


  // No terminal cost
  gmm::clear(p0[N_time-1]);
  gmm::clear(p1[N_time-1]);

  //  Region for the Dirichlet condition
  getfem::mesh_region r;
  getfem::outer_faces_of_mesh(mesh,r);
  int DIRICHLET_BOUNDARY_NUM = 11;
  for (getfem::mr_visitor i(r); !i.finished(); ++i) {
  
          base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
        un /= gmm::vect_norm2(un);
        //if ( gmm::abs(un[0]+1.0) < 1e-12 ) {
                mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(),i.f());
        //}
        //mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(),i.f());
  }


  // Assembly of the stiffness matrix with Dirichlet boundary condition
  getfem::model MS;
          //getfem::add_generic_elliptic_brick(MS, mim, "y", "pixdot");
  MS.add_fem_variable("p0", mf);
  MS.add_fem_variable("p1", mf);

  plain_vector Nux(nb_dof);
  for (size_type k=0; k<nb_dof; ++k) {
        Nux[k] = nux;
  }
  MS.add_initialized_fem_data("nux", mf, Nux);

  MS.add_initialized_fem_data("y0", mf, y0[N_time-1]);
  MS.add_initialized_fem_data("y1", mf, y1[N_time-1]);


  getfem::add_generic_elliptic_brick(MS, mim, "p1", "nux");
  plain_vector H(nb_dof);
  gmm::clear(H);
  MS.add_initialized_fem_data("Dirichletdata", mf, H);        
  add_Dirichlet_condition_with_simplification(MS, "p1", DIRICHLET_BOUNDARY_NUM, "Dirichletdata");


  MS.add_initialized_scalar_data("betaxxx", betax);
  MS.add_initialized_scalar_data("lambdaxxx", lambdax);
  MS.add_initialized_scalar_data("muxxx", mux);
  MS.add_macro("Ftens", "Id(meshdim)+Grad_y0");
  MS.add_macro("linear_tensor", "0.5*(Grad_p1'*(Id(meshdim)+Grad_y0) + (Id(meshdim)+Grad_y0)'*Grad_p1)");
  MS.add_macro("mytensor", "0.5*(Grad_y0'*Grad_y0+Grad_y0+Grad_y0')");
  getfem::add_linear_generic_assembly_brick(MS, mim, "betaxxx*((Grad_p1)*( 2.0*muxxx*mytensor )):Grad_Test_p0", -1, true, false);
  getfem::add_linear_generic_assembly_brick(MS, mim, "betaxxx*( lambdaxxx*Trace(mytensor) )*Grad_p1:Grad_Test_p0", -1, true, false);
  getfem::add_linear_generic_assembly_brick(MS, mim, "betaxxx*2.0*muxxx*(linear_tensor:(Ftens'*Grad_Test_p0))", -1, true, false);
  getfem::add_linear_generic_assembly_brick(MS, mim, "betaxxx*lambdaxxx*Trace(linear_tensor)*Trace(Ftens'*Grad_Test_p0)", -1, true, false);
  // Assembly of the matrix for the linearized term - linearized version of the nonlinear brick
          /*getfem::model MSA;
          MSA.add_fem_variable("q", mf);
          MSA.add_initialized_scalar_data("betaxxx", betax);
          MSA.add_initialized_fem_data("ygrec", mf, y[N_time-1]);        
          getfem::add_linear_generic_assembly_brick(MSA, mim, "-3.0*ygrec*ygrec*betaxxx*q*Test_q", -1, true, true);*/


  // Useful matrices
  sparse_matrix M(nb_dof, nb_dof);
  getfem::asm_mass_matrix(M, mim, mf, mf);

  getfem::add_explicit_matrix(MS, "p1", "p0", gmm::scaled(M, -1.0));

        //MS.assembly(getfem::model::BUILD_ALL);
        //sparse_matrix A(2*nb_dof, 2*nb_dof);
        //gmm::copy(MS.real_tangent_matrix(), A);
  //gmm::scale(A, -1.0);

  sparse_matrix A1(2*nb_dof, 2*nb_dof);
  gmm::clear(A1);

  gmm::iteration iter(1e-9, 0, 40000);
  scalar_type cpd = 0.0;

  plain_vector P0(nb_dof);
  gmm::clear(P0);
  gmm::copy(p0[N_time-1], P0);
  plain_vector P1(nb_dof);
  gmm::clear(P1);
  gmm::copy(p1[N_time-1], P1);

  plain_vector F(2*nb_dof); 
  gmm::clear(F);
  sparse_matrix MM(2*nb_dof, 2*nb_dof);
  gmm::clear(MM);
  plain_vector P(2*nb_dof);

  gmm::sub_interval I0;
  gmm::sub_interval I1;
  I0 = MS.interval_of_variable("p0");
  I1 = MS.interval_of_variable("p1");

  // Transient part I

  for (size_type i=N_time-1; i>N_time/2; --i) {

        timex -= step_dt;

        // Assembly of the unsteday matrix
        gmm::clear(A1);
        gmm::copy(y0[i-1], MS.set_real_variable("y0"));
        gmm::copy(y1[i-1], MS.set_real_variable("y1"));
                /*gmm::copy(y[i-1], MSA.set_real_variable("ygrec"));
                MSA.assembly(getfem::model::BUILD_ALL);
                gmm::copy(MSA.real_tangent_matrix(), A1);*/
        //asm_navier_stokes_nl_tgm(A1, mim, mf, y[i-1], getfem::mesh_region::all_convexes());
        //gmm::scale(A1, -betax);
        //gmm::add(A, A1);
        MS.assembly(getfem::model::BUILD_ALL);
        gmm::copy(MS.real_tangent_matrix(), A1);
        cpd = change_pi_dot(tau, timex - step_dt);
        gmm::scale(A1, cpd*step_dt);

        gmm::clear(MM);
        gmm::add(A1, MM);
        gmm::add(M, gmm::sub_matrix(MM, I0, I0));
        gmm::add(M, gmm::sub_matrix(MM, I1, I1));

        gmm::clear(F);
        gmm::mult(M, p0[i], gmm::sub_vector(F, I0));
        gmm::mult(M, p1[i], gmm::sub_vector(F, I1));

        // Solving
                  //cout << "Solving adjoint for time s = " << timex << " ..."  << endl;
        gmm::clear(P);
        gmm::MUMPS_solve(MM, P, F);
                  //cout << "Solved." << endl;

          // Extract the solution    
        gmm::copy(gmm::sub_vector(P, I0), p0[i-1]);
        gmm::copy(gmm::sub_vector(P, I1), p1[i-1]);

  }

  
  // Jump condition with the gradient of the functional to maximize
  plain_vector GradH(nb_dof);
  gmm::copy(y1[N_time/2-1], GradH);
        //getfem::compute_gradient(mf, mf, y1[N_time/2-1], GradH);
  for (size_type k=0; k<nb_dof; ++k) {

        GradH[k] += 0.5;
        base_node bn = mf.point_of_basic_dof(k);
        if ((bn[0] < 0.60) || (bn[0] > 0.90)) {
                GradH[k] = 0.0;
        }

  }

  gmm::scale(GradH, 1.0);

  gmm::clear(p0[N_time/2-2]);
  gmm::add(p0[N_time/2-1], GradH);
        //cout << "GradH = " << endl;
        //cout << GradH << endl;
        //cout << 1.0/0.0 << endl;
  gmm::copy(GradH, p0[N_time/2-2]);

  gmm::copy(p1[N_time/2-1], p1[N_time/2-2]);

  // Transient part II

  for (size_type i=N_time/2-2; i>0; --i) {

        timex -= step_dt;

        // Assembly of the unsteday matrix
        gmm::clear(A1);
        gmm::copy(y0[i-1], MS.set_real_variable("y0"));
        gmm::copy(y1[i-1], MS.set_real_variable("y1"));
                /*gmm::copy(y[i-1], MSA.set_real_variable("ygrec"));
                MSA.assembly(getfem::model::BUILD_ALL);
                gmm::copy(MSA.real_tangent_matrix(), A1);*/
        //asm_navier_stokes_nl_tgm(A1, mim, mf, y[i-1], getfem::mesh_region::all_convexes());
        //gmm::scale(A1, -betax);
        //gmm::add(A, A1);
        MS.assembly(getfem::model::BUILD_ALL);
        gmm::copy(MS.real_tangent_matrix(), A1);
        cpd = change_pi_dot(tau, timex - step_dt);
        gmm::scale(A1, cpd*step_dt);

        gmm::clear(MM);
        gmm::add(A1, MM);
        gmm::add(M, gmm::sub_matrix(MM, I0, I0));
        gmm::add(M, gmm::sub_matrix(MM, I1, I1));

        gmm::clear(F);
        gmm::mult(M, p0[i], gmm::sub_vector(F, I0));
        gmm::mult(M, p1[i], gmm::sub_vector(F, I1));

        // Solving
                  //cout << "Solving adjoint for time s = " << timex << " ..."  << endl;
        gmm::clear(P);
        gmm::MUMPS_solve(MM, P, F);
                  //cout << "Solved." << endl;

          // Extract the solution
                
        gmm::copy(gmm::sub_vector(P, I0), p0[i-1]);
        gmm::copy(gmm::sub_vector(P, I1), p1[i-1]);

  }


} 



















