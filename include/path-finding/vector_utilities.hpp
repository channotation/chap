#ifndef VECTOR_UTILITIES_HPP
#define VECTOR_UTILITIES_HPP

#include <cmath>

#include <gromacs/trajectoryanalysis.h>

/*
 * Turns 3-vector into unit vector pointing the same direction.
 */
void normaliseVector(gmx::RVec &vec);

#endif

