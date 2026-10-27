import numpy as np

from common.poisson_boltzmann import PoissonBoltzmannConfig, solve_poisson_boltzmann

def test_zero_surface_potential_is_neutral_solution():
    solution = solve_poisson_boltzmann(0.0)
    assert solution.converged
    assert solution.iterations == 1
    assert solution.residual_inf < 1e-12
    assert np.max(np.abs(solution.potential)) < 1e-14
    assert abs(solution.surface_charge) < 1e-20

def test_polarity_symmetry_and_finite_residual():
    config = PoissonBoltzmannConfig(points=101)
    positive = solve_poisson_boltzmann(0.12, config)
    negative = solve_poisson_boltzmann(-0.12, config)
    assert positive.converged and negative.converged
    assert positive.residual_inf < 1e-9 and negative.residual_inf < 1e-9
    np.testing.assert_allclose(negative.potential, -positive.potential, rtol=2e-8, atol=2e-10)
    assert np.isfinite(positive.surface_charge)
    np.testing.assert_allclose(negative.surface_charge, -positive.surface_charge, rtol=2e-7, atol=1e-12)

def test_grid_refinement_reduces_surface_charge_change():
    coarse = solve_poisson_boltzmann(0.08, PoissonBoltzmannConfig(points=101))
    fine = solve_poisson_boltzmann(0.08, PoissonBoltzmannConfig(points=201))
    finest = solve_poisson_boltzmann(0.08, PoissonBoltzmannConfig(points=401))
    assert coarse.converged and fine.converged and finest.converged
    coarse_change = abs(fine.surface_charge - coarse.surface_charge)
    fine_change = abs(finest.surface_charge - fine.surface_charge)
    assert fine_change < coarse_change
