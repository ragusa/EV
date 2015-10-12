/**
 * \file Exceptions.h
 * \brief Provides custom exception definitions.
 */
#ifndef Exceptions_h
#define Exceptions_h

using namespace dealii;

/** \brief Exception for a modulus not returning zero */
DeclException2(ExcModulusNotZero, int, int, << arg1 << " % " << arg2 << " != 0");

#endif
