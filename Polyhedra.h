#pragma once
////
//  polyhedra - regular polyhedra, symmetry groups, lattices.
//  D. P. Mitchell
//
#include "Affine.h"
//
//  Vertices of regular polyhedra, unit distance from the origin
//
extern const Vector3 ML_rgvTetrahedron[4];
extern const Vector3 ML_rgvCube[8];
extern const Vector3 ML_rgvOctahedron[6];
extern const Vector3 ML_rgvDodecahedron[20];
extern const Vector3 ML_rgvIcosahedron[12];
//
//  The polyhedral groups
//
extern const Matrix3 ML_rgmTetrahedralGroup[12];
extern const Matrix3 ML_rgmOctahedralGroup[24];
extern const Matrix3 ML_rgmIcosahedralGroup[60];
//
//  Basis lattice vectors, normalized to 1 sample per unit volume
//
extern const Vector3   ML_rgvSimpleCubic[3];
extern const Vector3   ML_rgvFaceCenteredCubic[3];
