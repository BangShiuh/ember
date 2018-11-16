#include "driftDiffusionSystem.h"
#include <assert.h>

DriftDiffusionSystem::DriftDiffusionSystem()
{
}

void DriftDiffusionSystem::get_A(dvec& a, dvec& b, dvec& c)
{
    DiffusionSystem::get_A(a, b, c);
    assert(mathUtils::notnan(mu));
    assert(mathUtils::notnan(E));
    assert(mathUtils::notnan(B));

    for (size_t j=0; j<=N-2; j++) {
        c3[j] = rphalf[j]*(mu[j]+mu[j+1])/2.0*E[j];
    }
    assert(mathUtils::notnan(c1));
    assert(mathUtils::notnan(c3));

    // Left boundary value
    size_t jStart;
    if (grid.leftBC == BoundaryCondition::FixedValue) {
        jStart = 1;
    } else if (grid.leftBC == BoundaryCondition::ControlVolume) {
        jStart =  1;
        double c0 = B[0] * (grid.alpha + 1) * (mu[0]+mu[1]) / (4 * hh[0]) * E[0];
        b[0] += c0;
        c[0] += c0;
    } else if (grid.leftBC == BoundaryCondition::WallFlux) {
        jStart = 1;
        double c0 = B[0] * (grid.alpha + 1) * (mu[0]+mu[1]) / (4 * hh[0]) * E[0];
        b[0] += c0;
        c[0] += c0;
    } else  { // (leftBC == BoundaryCondition::ZeroGradient)
        // In the case of a zero gradient boundary condition, the boundary value
        // is not computed, and the value one point in is computed by substituting
        // y[0] = y[1] in the finite difference formula.
        jStart = 2;
        b[1] += c1[1]*(c3[1] - 2.0*c3[0]);
        c[1] += c1[1]*c3[1];
    }

    // Right boundary value
    size_t jStop;
    if (grid.rightBC == BoundaryCondition::FixedValue) {
        jStop = N-1;
    } else { // (rightBC == BoundaryCondition::ZeroGradient)
        // In the case of a zero gradient boundary condition, the boundary value
        // is not computed, and the value one point in is computed by substituting
        // y[N-1] = y[N-2] in the finite difference formula.
        jStop = N-2;
        a[N-2] += -c1[N-2]*c3[N-3];
        b[N-2] += c1[N-2]*(2.0*c3[N-2] - c3[N-3]);
    }

    // Intermediate points
    for (size_t j=jStart; j<jStop; j++) {
        a[j] += -c1[j]*c3[j-1];
        b[j] += c1[j]*(c3[j] - c3[j-1]);
        c[j] += c1[j]*c3[j];
    }
    double aa;
    std::cin>>aa;
}

void DriftDiffusionSystem::resize(size_t N_)
{
    DiffusionSystem::resize(N_);
    N = N_;
    mu.setConstant(N, NaN);
    E.setConstant(N, NaN);
    c3.setConstant(N, 0);
}
