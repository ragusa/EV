/**
 * \file Exceptions.h
 * \brief Provides custom exception definitions.
 */
#ifndef Exceptions_h
#define Exceptions_h

using namespace dealii;

/** \brief Exception for a NaN being encountered */
DeclException0(ExcNaNEncountered);

/** \brief Exception for a modulus not returning zero */
DeclException2(ExcModulusNotZero, int, int, << arg1 << " % " << arg2 << " != 0");

/** \brief Exception for critical or supercritical flow being encountered */
DeclException1(
  ExcFlowNotSubcritical, int, << "Flow not subcritical: Fr = " << arg1 << " > 1");

#endif
