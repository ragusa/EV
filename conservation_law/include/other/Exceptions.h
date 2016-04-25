/**
 * \file Exceptions.h
 * \brief Provides custom exception definitions.
 */
#ifndef Exceptions_h
#define Exceptions_h

using namespace dealii;

/** \brief Exception for an assumption being violated */
DeclException0(ExcAssumptionViolated);

/** \brief Exception for a negative diagonal entry when not expected */
DeclException2(
  ExcNegativeDiagonal, int, double, << "Row " << arg1 << ": " << arg2);

/** \brief Exception for a NaN being encountered */
DeclException0(ExcNaNEncountered);

/** \brief Exception for trying to run a non-steady-state problem */
DeclException0(ExcNotASteadyStateProblem);

/** \brief Exception for a modulus not returning zero */
DeclException2(ExcModulusNotZero, int, int, << arg1 << " % " << arg2 << " != 0");

/** \brief Exception for critical or supercritical flow being encountered */
DeclException1(
  ExcFlowNotSubcritical, int, << "Flow not subcritical: Fr = " << arg1 << " > 1");

/** \brief Exception for having too many transient output files, which are
 *         numbered with 4 digits; therefore having 10001 output files is too
 *         many. */
DeclException0(ExcTooManyTransientOutputFiles);

/** \brief Exception for a file not existing */
DeclException1(
  ExcFileDoesNotExist, std::string, << "The file: " << arg1 << " does not exist");

/** \brief Exception for a directory not being able to be opened */
DeclException1(ExcDirectoryCannotBeOpened,
               std::string,
               << "The directory " << arg1 << " cannot be opened");

/** \brief Exception for sizes being inconsistent */
DeclException2(ExcSizesInconsistent,
               unsigned int,
               unsigned int,
               << "Size " << arg1 << " does not match size " << arg2);

/** \brief Exception for diffusion type being invalid */
DeclException0(ExcInvalidDiffusionType);

/** \brief Exception for negativity encountered */
DeclException2(
  ExcNegativity, std::string, double, << arg1 << " is negative: " << arg2);

/** \brief exception for reaching the maximum iteration */
DeclException1(
  ExcMaxIterationReached, unsigned int, << "Max iteration reached: " << arg1);

/** \brief Exception for antidiffusion upper bound being negative */
DeclException2(ExcAntidiffusionUpperBoundNegative,
               int,
               double,
               << "Q+[" << arg1 << "] = " << arg2 << " < 0");

/** \brief Exception for antidiffusion lower bound being positive */
DeclException2(ExcAntidiffusionLowerBoundPositive,
               int,
               double,
               << "Q-[" << arg1 << "] = " << arg2 << " > 0");

#endif
