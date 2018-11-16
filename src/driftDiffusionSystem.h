#pragma once

#include "diffusionSystem.h"

//! System representing diffusion of a single solution component.
//! Represents an ODE in one of the following forms:
//!     \f[ \dot{y} = B \frac{d}{dx}\left(D \frac{dy}{dx}\right) + C \f]
//!     \f[ \dot{y} = \frac{B}{r} \frac{d}{dr}\left(r D \frac{dy}{dr}\right) + C \f]
//! The ODE in this form may be written as a linear system:
//!     \f[ \dot{y} = Ay + k \f]
//! where the entries of the matrix *A* are determined by the prefactor *B*,
//! the diffusion coefficients *D*, the boundary conditions and the finite
//! difference formulas used.
class DriftDiffusionSystem : public DiffusionSystem
{
public:
    DriftDiffusionSystem();

    //! Build the matrix *A* describing the linear %ODE. `a[j]`, `b[j]`, and
    //! `c[j]` are respectivley the subdiagonal, diagonal, and superdiagonal
    //! elements of row `j`.
    virtual void get_A(dvec& a, dvec& b, dvec& c);

    virtual void resize(size_t N);

private:
    // constants used to compute the matrix coefficients
    dvec c3;
};
