#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

using namespace dealii;

int main(int, char **)
{
  try
  {
    dealii::deallog.depth_console(0);

    // create FE
    const FE_Q<1> fe(1);

    // create triangulation
    const double domain_length = 10.0;
    Triangulation<1> triangulation;
    GridGenerator::hyper_cube(triangulation, 0.0, domain_length);
    triangulation.refine_global(3);
    const unsigned int n_cells = triangulation.n_active_cells();
    const double dx = domain_length / n_cells;

    // create DoF handler
    DoFHandler<1> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe);

    // create dummy quadrature for evaluation point
    const double eval_point = 7.3;
    const unsigned int eval_cell = 5; // cell number of evaluation point
    const double unit_cell_position = (eval_point - eval_cell * dx) / dx;
    const std::vector<Point<1>> point(1, Point<1>(unit_cell_position));
    const Quadrature<1> quadrature(point);

    // create FE values
    FEValues<1> fe_values(fe, quadrature, update_values);

    // create solution vector
    Vector<double> solution(dof_handler.n_dofs());
    const double left_value = 32.4;
    const double right_value = 12.8;
    solution[eval_cell] = left_value;
    solution[eval_cell + 1] = right_value;

    // loop over cells to evaluate solution at evaluation point
    typename DoFHandler<1>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
    for (unsigned int i_cell = 0; cell != endc; ++cell, ++i_cell)
    {
      if (i_cell == eval_cell)
      {
        // compute actual value
        fe_values.reinit(cell);
        std::vector<double> solution_values(1);
        fe_values.get_function_values(solution, solution_values);

        // compute expected value
        const double left_basis_value = 1.0 - unit_cell_position;
        const double right_basis_value = unit_cell_position;
        const double expected_value =
          left_value * left_basis_value + right_value * right_basis_value;

        std::cout << "Evaluation point: " << eval_point << std::endl;
        std::cout << "Left cell point: " << i_cell * dx << std::endl;
        std::cout << "Right cell point: " << (i_cell + 1) * dx << std::endl;
        std::cout << "Left solution value: " << left_value << std::endl;
        std::cout << "Right solution value: " << right_value << std::endl;
        std::cout << "Expected value: " << expected_value << std::endl;
        std::cout << "Actual value: " << solution_values[0] << std::endl;
      }
    }
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
