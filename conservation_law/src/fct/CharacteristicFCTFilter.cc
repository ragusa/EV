/**
 * \file CharacteristicFCTFilter.cc
 * \brief Provides the function definitions for the CharacteristicFCTFilter
 * class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 * \param[in] limiter_  limiter
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 * \param[in] lumped_mass_matrix_  lumped mass matrix \f$\mathbf{M}^L\f$
 */
template <int dim>
CharacteristicFCTFilter<dim>::CharacteristicFCTFilter(
  const RunParameters & run_parameters_,
  const std::shared_ptr<Limiter<dim>> limiter_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const SparseMatrix<double> & lumped_mass_matrix_)
  : ExplicitEulerFCTFilter<dim>(
      run_parameters_, limiter_, dof_handler_, fe_, lumped_mass_matrix_)
{
  // create sparsity pattern for limiter and antidiffusion matrices
  DynamicSparsityPattern dsp(this->n_dofs);
  DoFTools::make_sparsity_pattern(*this->dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  // initialize matrix with sparsity pattern
  characteristic_antidiffusion_matrix.reinit(sparsity_pattern);

  // resize vectors
  old_solution_characteristic.reinit(this->n_dofs);
  tmp_vector.reinit(this->n_dofs);

  // create lists of DoF indices
  create_dof_indices_lists();
}

/**
 * \brief Filters antidiffusive fluxes.
 *
 * \param[in] old_solution  old solution \f$\mathbf{U}^n\f$
 * \param[in] dt  time step size \f$\Delta t\f$
 * \param[in] inviscid_ss_flux  inviscid steady-state flux vector
 *            (entries are \f$(\mathbf{A}\mathbf{U}^n)_i\f$ for scalar case,
 *            \f$\sum\limits_j\mathbf{c}_{i,j}\cdot\mathrm{F}^n_j\f$ for
 *            systems case)
 * \param[in] ss_reaction  steady-state reaction vector \f$\mathbf{\sigma}\f$,
 *            where \f$\sigma_i = \int\limits_{S_i}
 *            \varphi_i(\mathbf{x})\sigma(\mathbf{x})dV
 *            = \sum\limits_j\int\limits_{S_{i,j}}
 *            \varphi_i(\mathbf{x})\varphi_j(\mathbf{x})\sigma(\mathbf{x})dV \f$
 * \param[in] low_order_diffusion_matrix  low-order diffusion matrix
 *            \f$\mathbf{D}^L\f$
 * \param[in] ss_rhs  steady-state right hand side vector \f$\mathbf{b}^n\f$
 * \param[inout] limiter_matrix  limiter matrix \f$\mathbf{L}\f$
 * \param[inout] antidiffusion_matrix  antidiffusion matrix \f$\mathbf{P}\f$
 */
template <int dim>
void CharacteristicFCTFilter<dim>::filter_antidiffusive_fluxes(
  const Vector<double> & old_solution,
  const double & dt,
  const Vector<double> & inviscid_ss_flux,
  const Vector<double> & ss_reaction,
  const SparseMatrix<double> & low_order_diffusion_matrix,
  const Vector<double> & ss_rhs,
  SparseMatrix<double> & limiter_matrix,
  SparseMatrix<double> & antidiffusion_matrix)
{
  // compute characteristic solution bounds \hat{W}- and \hat{W}+
  compute_solution_bounds(old_solution, dt, ss_reaction, ss_rhs);

  // compute characteristic antidiffusion bounds \hat{Q}- and \hat{Q}+
  compute_antidiffusion_bounds(
    old_solution, dt, inviscid_ss_flux, low_order_diffusion_matrix, ss_rhs);

  // enforce antidiffusion bounds signs if requested
  if (this->do_enforce_antidiffusion_bounds_signs)
    this->enforce_antidiffusion_bounds_signs();

  // transform antidiffusion matrix: P -> \hat{P}
  transform_matrix(
    old_solution, antidiffusion_matrix, characteristic_antidiffusion_matrix);

  // limit antidiffusion fluxes
  this->limiter->compute_limiter_matrix(characteristic_antidiffusion_matrix,
                                        this->antidiffusion_bounds,
                                        limiter_matrix);
  this->limiter->apply_limiter_matrix(limiter_matrix, antidiffusion_matrix);
}

/**
 * \brief Checks to see if the FCT bounds were satisfied.
 *
 * \param[in] new_solution  new solution vector \f$\mathbf{U}^{n+1}\f$
 *
 * \return flag that FCT bounds were satisfied for all filters
 */
template <int dim>
bool CharacteristicFCTFilter<dim>::check_bounds(
  const Vector<double> & new_solution)
{
  // create reference
  Vector<double> & solution_characteristic = old_solution_characteristic;

  // transform solution to characteristic variables
  transform_vector(new_solution, new_solution, solution_characteristic);

  // check bounds
  const bool bounds_satisfied =
    this->solution_bounds.check_bounds(solution_characteristic);

  return bounds_satisfied;
}

/**
 * \brief Computes bounds to be imposed on the characteristic FCT solution,
 *        \f$\hat{W}_i^-\f$ and \f$\hat{W}_i^+\f$.
 *
 * The bounds are computed as
 * \f[
 *   \hat{W}_i^-(\mathrm{U}^n) = \hat{U}_{min,i}^n \,,
 * \f]
 * \f[
 *   \hat{W}_i^+(\mathrm{U}^n) = \hat{U}_{max,i}^n \,,
 * \f]
 * where \f$\hat{U}_{min,i}\equiv\min\limits_{j\in \mathcal{I}(S_i)} \hat{U}_j\f$
 * and \f$\hat{U}_{max,i}\equiv\max\limits_{j\in \mathcal{I}(S_i)} \hat{U}_j\f$.
 *
 * \warning These bounds currently assume
 * - no reaction term
 * - no source term
 * To use these terms, it is necessary to perform appropriate transformations
 * on them.
 *
 * \param[in] old_solution  old solution vector \f$\mathbf{U}^n\f$
 * \param[in] dt  time step size \f$\Delta t\f$
 * \param[in] ss_reaction  steady-state reaction vector \f$\mathbf{\sigma}\f$,
 *            where \f$\sigma_i = \int\limits_{S_i}
 *            \varphi_i(\mathbf{x})\sigma(\mathbf{x})dV
 *            = \sum\limits_j\int\limits_{S_{i,j}}
 *            \varphi_i(\mathbf{x})\varphi_j(\mathbf{x})\sigma(\mathbf{x})dV \f$
 * \param[in] ss_rhs  old steady-state right hand side vector \f$\mathbf{b}^n\f$
 */
template <int dim>
void CharacteristicFCTFilter<dim>::compute_solution_bounds(
  const Vector<double> & old_solution,
  const double &,
  const Vector<double> &,
  const Vector<double> &)
{
  // create references
  Vector<double> & solution_min = this->solution_bounds.lower;
  Vector<double> & solution_max = this->solution_bounds.upper;

  // transform old solution vector to characteristic variables
  transform_vector(old_solution, old_solution, old_solution_characteristic);

  // compute minimum and maximum values of solution
  this->compute_min_and_max_of_dof_vector(
    old_solution_characteristic, solution_min, solution_max);
}

/**
 * \brief Computes the characteristic antidiffusion bounds
 *        \f$\hat{\mathbf{Q}}^\pm\f$.
 *
 * \param[in] old_solution  old solution \f$\mathbf{U}^n\f$
 * \param[in] dt  time step size \f$\Delta t\f$
 * \param[in] inviscid_ss_flux  inviscid steady-state flux vector
 *   (entries are \f$(\mathbf{A}\mathbf{U}^n)_i\f$ for scalar case,
 *   \f$\sum\limits_j\mathbf{c}_{i,j}\cdot\mathrm{F}^n_j\f$ for systems case)
 * \param[in] low_order_diffusion_matrix  low-order diffusion matrix
 *            \f$\mathbf{D}^L\f$
 * \param[in] ss_rhs  steady-state right hand side vector \f$\mathbf{b}^n\f$
 */
template <int dim>
void CharacteristicFCTFilter<dim>::compute_antidiffusion_bounds(
  const Vector<double> & old_solution,
  const double & dt,
  const Vector<double> &,
  const SparseMatrix<double> &,
  const Vector<double> &)
{
  // create references
  Vector<double> & Q_minus = this->antidiffusion_bounds.lower;
  Vector<double> & Q_plus = this->antidiffusion_bounds.upper;
  const Vector<double> & solution_min = this->solution_bounds.lower;
  const Vector<double> & solution_max = this->solution_bounds.upper;
  Vector<double> & tmp = this->tmp_vector;
  const SparseMatrix<double> & lumped_mass_matrix = *this->lumped_mass_matrix;

  // start computing Q+
  Q_plus = 0;
  lumped_mass_matrix.vmult(tmp, old_solution);
  Q_plus.add(-1.0 / dt, tmp);
  // Q_plus.add(1.0, inviscid_ss_flux);
  // low_order_diffusion_matrix.vmult(tmp, old_solution);
  // Q_plus.add(1.0, tmp);
  // Q_plus.add(-1.0, ss_rhs);

  // copy current contents of Q+ as these components are identical
  Q_minus = Q_plus;

  // finish computing Q+ and Q-
  lumped_mass_matrix.vmult(tmp, solution_max);
  Q_plus.add(1.0 / dt, tmp);
  lumped_mass_matrix.vmult(tmp, solution_min);
  Q_minus.add(1.0 / dt, tmp);
}

/**
 * \brief Transforms a vector in terms of conservative variables to a vector
 *        in terms of characteristic variables.
 *
 * This function applies a local transformation evaluated with the solution
 * \f$\mathbf{U}\f$ to characteristic variables on a vector \f$\mathbf{y}\f$:
 * \f[
 *   \hat{\mathbf{y}}_i = \mathbf{T}(\mathbf{U}_i)\mathbf{y}_i \,,
 * \f]
 * where \f$\hat{\mathbf{y}}\f$ is the transformed vector, and \f$i\f$ is a
 * node index.
 *
 * \param[in] solution  solution vector \f$\mathbf{U}\f$ at which to evaluate
 *                      transformations
 * \param[in] vector_original  vector \f$\mathbf{y}\f$ to be transformed
 * \param[in] vector_transformed  transformed vector \f$\hat{\mathbf{y}}\f$
 */
template <int dim>
void CharacteristicFCTFilter<dim>::transform_vector(
  const Vector<double> & solution,
  const Vector<double> & vector_original,
  Vector<double> & vector_transformed) const
{
  // loop over nodes
  for (unsigned int n = 0; n < this->n_nodes; ++n)
  {
    // assemble solution vector on node
    Vector<double> solution_node(this->n_components);
    for (unsigned int m = 0; m < this->n_components; ++m)
      solution_node[m] = solution[node_dof_indices[n][m]];

    // compute local transformation matrix inverse on node
    const FullMatrix<double> transformation_matrix =
      compute_transformation_matrix(solution_node);

    // apply local transformation matrix inverse to vector
    // loop over rows of matrix
    for (unsigned int m = 0; m < this->n_components; ++m)
    {
      // initialize sum for matrix-vector product for row
      const unsigned int i = node_dof_indices[n][m];
      vector_transformed[i] = 0.0;

      // loop over columns of matrix to compute matrix-vector product for row
      for (unsigned int mm = 0; mm < this->n_components; ++mm)
      {
        const unsigned int j = node_dof_indices[n][mm];
        vector_transformed[i] +=
          transformation_matrix[m][mm] * vector_original[j];
      }
    }
  }
}

/**
 * \brief Transforms a matrix in terms of conservative variables to a matrix
 *        in terms of characteristic variables.
 *
 * This function applies a local transformation evaluated with the solution
 * \f$\mathbf{U}\f$ to characteristic variables on a matrix \f$\mathbf{A}\f$:
 * \f[
 *   \hat{\mathbf{A}}_{i,j} = \mathbf{T}^{-1}(\mathbf{U}_i)\mathbf{A}_{i,j} \,,
 * \f]
 * where \f$\hat{\mathbf{A}}\f$ is the transformed matrix, and \f$i\f$ and
 * \f$j\f$ are node indices.
 *
 * \param[in] solution  solution vector \f$\mathbf{U}\f$ at which to evaluate
 *                      transformations
 * \param[in] matrix_original  matrix \f$\mathbf{A}\f$ to be transformed
 * \param[in] matrix_transformed  transformed matrix \f$\hat{\mathbf{A}}\f$
 */
template <int dim>
void CharacteristicFCTFilter<dim>::transform_matrix(
  const Vector<double> & solution,
  const SparseMatrix<double> & matrix_original,
  SparseMatrix<double> & matrix_transformed) const
{
  // loop over nodes
  for (unsigned int n_i = 0; n_i < this->n_nodes; ++n_i)
  {
    // assemble solution vector on node
    Vector<double> solution_node(this->n_components);
    for (unsigned int m = 0; m < this->n_components; ++m)
      solution_node[m] = solution[node_dof_indices[n_i][m]];

    // compute local transformation matrix inverse on node
    const FullMatrix<double> transformation_matrix =
      compute_transformation_matrix(solution_node);

    // apply local transformation matrix inverse to vector
    // loop over components
    for (unsigned int m = 0; m < this->n_components; ++m)
    {
      // row index
      const unsigned int i = node_dof_indices[n_i][m];

      // loop over nonzero columns in row
      SparseMatrix<double>::const_iterator it_in = matrix_original.begin(i);
      SparseMatrix<double>::const_iterator it_end = matrix_original.end(i);
      SparseMatrix<double>::iterator it_out = matrix_transformed.begin(i);
      for (; it_in != it_end; ++it_in, ++it_out)
      {
        // get column index
        const unsigned int j = it_in->column();

        // get original matrix value for column
        const double value = it_in->value();

        // assert that if the value is greater than zero, then i and j correspond
        // to the same solution component
        Assert(std::abs(value) < 1.0e-15 || component_indices[j] == m,
               ExcInvalidState());

        // get node index corresponding to column
        const unsigned int n_j = node_indices[j];

        // compute transformed matrix value
        double transformed_value = 0.0;
        for (unsigned int l = 0; l < this->n_components; ++l)
        {
          const unsigned int k = node_dof_indices[n_i][l];
          const unsigned int p = node_dof_indices[n_j][l];
          const double value_original = matrix_original(k, p);

          transformed_value += transformation_matrix[m][l] * value_original;
        }

        // store transformed matrix value
        it_out->value() = transformed_value;
      }
    }
  }
}

/**
 * \brief Creates lists of degree of freedom indices lists such as node list,
 *        component list, and nodal degrees of freedom.
 *
 * \c deal.II documentation states that \e local degree of freedom indices have a
 * a standard ordering: vertex 0 DoFs, vertex 1 DoFs, etc. Within each vertex,
 * it is assumed here that local ordering is by component index.
 */
template <int dim>
void CharacteristicFCTFilter<dim>::create_dof_indices_lists()
{
  // initialize
  node_indices.resize(this->n_dofs);
  component_indices.resize(this->n_dofs);
  node_dof_indices.clear();

  // create vector for flagging whether a node index has been assigned or not
  std::vector<bool> assigned_node(this->n_dofs, false);

  // global degree of freedom indices in cell
  std::vector<unsigned int> local_dof_indices(this->dofs_per_cell);

  unsigned int node_index = 0; // current node index
  Cell cell = this->dof_handler->begin_active();
  Cell endc = this->dof_handler->end();
  for (; cell != endc; ++cell)
  {
    // get list of degrees of freedom on cell
    cell->get_dof_indices(local_dof_indices);

    // loop over nodes on cell
    for (unsigned int i_node = 0; i_node < this->dofs_per_cell_per_component;
         ++i_node)
    {
      // get DoF index of first DoF on node
      const unsigned int i_local0 = i_node * this->n_components;
      const unsigned int i0 = local_dof_indices[i_local0];

      // check if this DoF already has been assigned a node index
      if (!assigned_node[i0])
      {
        // vector of DoF indices on this node
        std::vector<unsigned int> this_node_dof_indices(this->n_components);

        // loop over DoFs on node
        for (unsigned int m = 0; m < this->n_components; ++m)
        {
          // determine global DoF index for component
          const unsigned int i_local = i_node * this->n_components + m;
          const unsigned int i = local_dof_indices[i_local];

          // assign node index for DoF
          node_indices[i] = node_index;

          // assign component index for DoF
          component_indices[i] = m;

          // assign DoF index to list of DoF indices on node
          this_node_dof_indices[m] = i;
        }

        // push back nodal DoF indices list
        node_dof_indices.push_back(this_node_dof_indices);

        // flag that node index has been assigned for the first DoF on node
        assigned_node[i0] = true;

        // increment node index
        node_index++;
      }
    }
  }

  // keep number of nodes
  n_nodes = node_index;
}
