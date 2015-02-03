/** \brief Tests entropy viscosity to compare against MATLAB results
 */
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <sstream>
#include <cstdlib>
#include "../EntropyViscosity.h"

using namespace dealii;

const unsigned int dim = 1; // number of spatial dimensions

int main(int argc, char ** argv) {
   try {
      dealii::deallog.depth_console(0);

      // create data necessary for constructor of entropy viscosity class
      double dt = 0.01;
      unsigned int degree = 1;
      unsigned int refinement_level = 5; // n_cells = 2^refinement_level
      double mu_x = 0.3;
      std::string cross_section_string = "2";
      std::string source_string = "3";
      unsigned int n_quadrature_points = 3;
      std::string entropy_string = "0.5*u*u";
      double entropy_residual_coefficient = 1.0;
      double jump_coefficient = 0.0;
      double x_min = 0.0;
      double x_max = 10.0;

      // dependent input variables
      Triangulation<dim> triangulation;
      GridGenerator::hyper_cube(triangulation, x_min, x_max);
      double domain_volume = std::pow((x_max-x_min),dim);
      triangulation.refine_global(refinement_level);
      unsigned int n_cells = triangulation.n_active_cells();

      Vector<double> new_solution(n_cells+1);
      Vector<double> old_solution(n_cells+1);
      double x = x_min;
      double dx = (x_max - x_min)/n_cells;
      for (unsigned int i = 0; i < n_cells+1; ++i) {
         new_solution(i) = std::sin(2*x);
         old_solution(i) = std::sin(x);
         x += dx;
      }

      FESystem<dim> fe(FE_Q<dim>(degree),1);

      DoFHandler<dim> dof_handler(triangulation);
      dof_handler.distribute_dofs(fe);

      QGauss<dim>   cell_quadrature(n_quadrature_points);
      QGauss<dim-1> face_quadrature(n_quadrature_points);

      Tensor<1,dim> transport_direction(0.0);
      transport_direction[0] = mu_x;

      FunctionParser<dim> cross_section_function;
      FunctionParser<dim> source_function;
      std::map<std::string,double> constants;
      cross_section_function.initialize("x",  cross_section_string,constants,false);
      source_function       .initialize("x,t",source_string,       constants,true);

      // instantiate entropy viscosity object
      EntropyViscosity<dim> EV(fe,
                               n_cells,
                               dof_handler,
                               cell_quadrature,
                               face_quadrature,
                               transport_direction,
                               cross_section_function,
                               source_function,
                               entropy_string,
                               entropy_residual_coefficient,
                               jump_coefficient,
                               domain_volume);                      
      
      // compute entropy viscosity
      Vector<double> entropy_viscosity
         = EV.compute_entropy_viscosity(new_solution,old_solution,dt);

      // output entropy viscosity
      DataOut<dim> visc_out;
      visc_out.attach_dof_handler(dof_handler);
      visc_out.add_data_vector(entropy_viscosity,"Entropy_Viscosity",DataOut<dim>::type_cell_data);
      std::ofstream viscosity_outstream("output/entropy_viscosity_dealii.txt");
      visc_out.build_patches(degree + 1);
      visc_out.write_gnuplot(viscosity_outstream);

      // output solutions
      DataOut<dim> solution_out;
      solution_out.attach_dof_handler(dof_handler);
      solution_out.add_data_vector(new_solution,"new_solution");
      solution_out.add_data_vector(old_solution,"old_solution");
      std::ofstream solutions_outstream("output/solutions_dealii.txt");
      solution_out.build_patches(degree + 1);
      solution_out.write_gnuplot(solutions_outstream);

   } catch (std::exception &exc) {
      std::cerr << std::endl << std::endl
            << "----------------------------------------------------"
            << std::endl;
      std::cerr << "Exception on processing: " << std::endl << exc.what()
            << std::endl << "Aborting!" << std::endl
            << "----------------------------------------------------"
            << std::endl;
      return 1;
   } catch (...) {
      std::cerr << std::endl << std::endl
            << "----------------------------------------------------"
            << std::endl;
      std::cerr << "Unknown exception!" << std::endl << "Aborting!" << std::endl
            << "----------------------------------------------------"
            << std::endl;
      return 1;
   };

   return 0;
}
