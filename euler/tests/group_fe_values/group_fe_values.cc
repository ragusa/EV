#include <algorithm>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/logstream.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/numerics/vector_tools.h>
#include "../../GroupFEValues.h"

using namespace dealii;

DeclException2(ExcSizeMismatch, int, int,
  << "length1 = " << arg1 << ", length2 = " << arg2);

double my_function(const std::vector<double> & solution)
{
  double value = 1.0;
  for (unsigned int i = 0; i < solution.size(); ++i)
    value *= solution[i];
  return value;
}

template <int dim>
class MyFunctionFEValues : public GroupFEValues<dim>
{
public:
  MyFunctionFEValues(
    const unsigned int & n_components,
    const Vector<double> & solution,
    const Triangulation<dim> & triangulation,
    const QGauss<dim> & cell_quadrature) :
    GroupFEValues<dim>(n_components, solution, triangulation, cell_quadrature)
  {
    this->compute_function_dof_values();
  }

private:
  double function(const std::vector<double> & solution) const override
  {
    double value = 1.0;
    for (unsigned int i = 0; i < solution.size(); ++i)
      value *= solution[i];
    return value;
  }
};

template <int dim>
void test()
{
  // create triangulation
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
  const unsigned int refinement_level = 3;
  triangulation.refine_global(refinement_level);

  // create finite elements
  const unsigned int degree = 1;
  const unsigned int n_components = dim + 1;
  const FESystem<dim> fe(FE_Q<dim>(degree), n_components);
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int dofs_per_cell_per_component = dofs_per_cell / n_components;

  // create DoF handler
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  const unsigned int n_dofs = dof_handler.n_dofs();
  const unsigned int n_dofs_per_component = n_dofs / n_components;

  // create function for solution initialization
  std::map<std::string, double> constants;
  std::vector<std::string> solution_strings(n_components);
  solution_strings[0] = "x";
  for (unsigned int d = 0; d < dim; ++d)
    solution_strings[d+1] = "x^" + std::to_string(d+2);
  FunctionParser<dim> solution_function(n_components);
  solution_function.initialize(
    FunctionParser<dim>::default_variable_names(),
    solution_strings,
    constants,
    false);

  // create and initialize solution vector
  Vector<double> solution(n_dofs);
  VectorTools::interpolate(dof_handler, solution_function, solution);

  // create quadrature
  const unsigned int n_q_points_per_dim = 3;
  const QGauss<dim> quadrature(n_q_points_per_dim);
  const unsigned int n_q_points_cell = quadrature.size();

  // initialize group FE values
  MyFunctionFEValues<dim> my_function_fe_values(
    n_components,
    solution,
    triangulation,
    quadrature);
  std::vector<double> function_dofs =
    my_function_fe_values.get_function_dof_values();

  // get mapping of DoFs to support points
  const MappingQ1<dim> mapping_q1;
  std::map<types::global_dof_index, Point<dim>> support_point_map;
  DoFTools::map_dofs_to_support_points(
    mapping_q1, dof_handler, support_point_map);
  
  // get vector of support points on unit cell, which includes duplicates for
  // each component
  std::vector<Point<dim>> unit_support_points = fe.get_unit_support_points();
  Assert(unit_support_points.size() == dofs_per_cell,
    ExcSizeMismatch(unit_support_points.size(), dofs_per_cell));
  // construct vector with unique support points on unit cell
  std::vector<Point<dim>> unique_unit_support_points;
  for (unsigned int i = 0; i < dofs_per_cell; ++i)
    if (std::find(
      unique_unit_support_points.begin(),
      unique_unit_support_points.end(),
      unit_support_points[i]) == unique_unit_support_points.end())
      unique_unit_support_points.push_back(unit_support_points[i]);
  // assert that the vector of unique support points has expected size
  Assert(unique_unit_support_points.size() == dofs_per_cell_per_component,
    ExcSizeMismatch(
      unique_unit_support_points.size(), dofs_per_cell_per_component));

  // mapping for transforming to real points
  MappingQ<dim> mapping_q(degree);

  // vector of unique real points
  std::vector<Point<dim>> unique_real_support_points;

  // loop over cells
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // loop over unique support points in cell and create vector of unique
    // real support points
    for (unsigned int i = 0; i < dofs_per_cell_per_component; ++i)
    {
      // get real point and add to vector if not there already
      Point<dim> real_support_point = mapping_q.transform_unit_to_real_cell(
        cell, unique_unit_support_points[i]);
      if (std::find(
        unique_real_support_points.begin(),
        unique_real_support_points.end(),
        real_support_point) == unique_real_support_points.end())
        unique_real_support_points.push_back(real_support_point);
    }
  }

  // assert that number of unique real support points is equal to the number
  // of DoFs per component
  Assert(unique_real_support_points.size() == n_dofs_per_component,
    ExcSizeMismatch(
      unique_real_support_points.size(), n_dofs_per_component));

  // print table of function at DoF support points
  //---------------------------------------------------------------------------
  // print header
  std::cout << "DoF support points:" << std::endl;
  std::string point_header;
  if (dim == 1)
    std::cout << "     x";
  else if (dim == 2)
    std::cout << "     x      y";
  else
    std::cout << "     x      y      z";
  for (unsigned int i = 0; i < n_components; ++i)
  {
    std::stringstream solution_header;
    solution_header << "u[" << i << "]";
    printf("%11s", solution_header.str().c_str());
  }
  printf("%11s", "f_correct");
  printf("%11s", "f_actual");
  std::cout << std::endl;
  // loop over support points
  for (unsigned int i = 0; i < n_dofs_per_component; ++i)
  {
    // reference the point
    Point<dim> & point = unique_real_support_points[i];

    // print the point
    std::stringstream point_stringstream;
    for (unsigned int d = 0; d < dim; ++d)
    {
      printf("%6.3f", point[d]);
      if (d != dim-1)
        std::cout << " ";
    }
    
    // compute each component of solution at point
    std::vector<double> solution_at_point_correct(n_components);
    for (unsigned int j = 0; j < n_components; ++j)
    {
      solution_at_point_correct[j] = solution_function.value(point, j);
      printf("%11.3e", solution_at_point_correct[j]);
    }

    // compute and print the correct function value
    double function_value_correct = my_function(solution_at_point_correct);
    printf("%11.3e", function_value_correct);

    // print the actual function value
    printf("%11.3e", function_dofs[i]);

    std::cout << std::endl;
  }

  // print table of function at quadrature points
  //---------------------------------------------------------------------------
  // create FE values
  FE_Q<dim> fe_scalar(degree);
  FEValues<dim> fe_values_scalar(fe_scalar, quadrature, update_values);
  DoFHandler<dim> dof_handler_scalar(triangulation);
  dof_handler_scalar.distribute_dofs(fe_scalar);
  typename DoFHandler<dim>::active_cell_iterator
    cell_scalar = dof_handler_scalar.begin_active(),
    endc_scalar = dof_handler_scalar.end();

  // print header
  std::cout << std::endl << "Quadrature points:" << std::endl;
  printf("Cell");
  for (unsigned int q = 0; q < n_q_points_cell; ++q)
    printf("  f_correct[%i]  f_actual[%i]",q,q);
  printf("\n");

  unsigned int i_cell = 0;
  for (cell = dof_handler.begin_active(); cell != endc; ++cell, ++cell_scalar)
  {
    // update FE values
    fe_values_scalar.reinit(cell_scalar);

    // get DoF indices
    std::vector<unsigned int> local_dof_indices(dofs_per_cell_per_component);
    cell_scalar->get_dof_indices(local_dof_indices);

    // compute the expected function values at quadrature points
    std::vector<double> function_values_correct(n_q_points_cell,0.0);
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
      for (unsigned int i = 0; i < dofs_per_cell_per_component; ++i)
        function_values_correct[q] += function_dofs[local_dof_indices[i]]
          * fe_values_scalar.shape_value(i,q);

    // get the actual function values at quadrature points
    std::vector<double> function_values =
      my_function_fe_values.get_function_values(cell);

    // print values
    printf("%4i",i_cell);
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
      printf("%14.3e  %11.3e", function_values_correct[q], function_values[q]);
    std::cout << std::endl;

    // increment cell index
    i_cell++;
  }
}

int main(int, char **)
{
  try
  {
    dealii::deallog.depth_console(0);

    test<1>();
  }
  catch (std::exception & exc)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }

  return 0;
}
