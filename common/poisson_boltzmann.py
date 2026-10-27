"""One-dimensional Poisson--Boltzmann solver for the S1 baseline."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

Q_E = 1.602176634e-19
EPSILON_0 = 8.8541878128e-12
K_B = 1.380649e-23

@dataclass(frozen=True)
class PoissonBoltzmannConfig:
    temperature: float = 300.0
    intrinsic_density: float = 1.0e16
    relative_permittivity: float = 11.7
    domain_debye_lengths: float = 20.0
    points: int = 201

    def __post_init__(self) -> None:
        if self.temperature <= 0 or self.intrinsic_density <= 0:
            raise ValueError("temperature and intrinsic density must be positive")
        if self.relative_permittivity <= 0 or self.domain_debye_lengths <= 0:
            raise ValueError("material and domain parameters must be positive")
        if self.points < 3 or self.points % 2 == 0:
            raise ValueError("points must be an odd integer >= 3")

    @property
    def thermal_voltage(self) -> float:
        return K_B * self.temperature / Q_E

    @property
    def permittivity(self) -> float:
        return self.relative_permittivity * EPSILON_0

    @property
    def debye_length(self) -> float:
        return np.sqrt(self.permittivity * self.thermal_voltage / (Q_E * self.intrinsic_density))

    @property
    def length(self) -> float:
        return self.domain_debye_lengths * self.debye_length

@dataclass(frozen=True)
class PoissonBoltzmannSolution:
    x: np.ndarray
    potential: np.ndarray
    surface_potential: float
    surface_charge: float
    residual_inf: float
    iterations: int
    converged: bool

def _residual(u: np.ndarray, dx: float, debye_length: float) -> np.ndarray:
    return (u[:-2] - 2.0 * u[1:-1] + u[2:]) / dx**2 - np.sinh(u[1:-1]) / debye_length**2

def solve_poisson_boltzmann(surface_potential: float, config: PoissonBoltzmannConfig | None = None,
                            *, max_iterations: int = 80, tolerance: float = 1.0e-10) -> PoissonBoltzmannSolution:
    """Solve u''=sinh(u)/L_D^2 with u(0)=V_s/V_T and u(L)=0."""
    cfg = config or PoissonBoltzmannConfig()
    if not np.isfinite(surface_potential):
        raise ValueError("surface potential must be finite")
    if max_iterations < 1 or tolerance <= 0:
        raise ValueError("invalid Newton controls")
    x = np.linspace(0.0, cfg.length, cfg.points)
    dx = x[1] - x[0]
    u = np.linspace(surface_potential / cfg.thermal_voltage, 0.0, cfg.points)
    converged = False
    iterations = 0
    for iterations in range(1, max_iterations + 1):
        residual = _residual(u, dx, cfg.debye_length)
        norm = float(np.max(np.abs(residual))) * cfg.debye_length**2
        if norm < tolerance:
            converged = True
            break
        interior = u[1:-1]
        n = interior.size
        jacobian = np.diag(-2.0 / dx**2 - np.cosh(interior) / cfg.debye_length**2)
        off = np.full(n - 1, 1.0 / dx**2)
        jacobian += np.diag(off, 1) + np.diag(off, -1)
        delta = np.linalg.solve(jacobian, -residual)
        step = 1.0
        while step >= 2.0**-20:
            trial = u.copy()
            trial[1:-1] += step * delta
            trial_norm = float(np.max(np.abs(_residual(trial, dx, cfg.debye_length)))) * cfg.debye_length**2
            if np.isfinite(trial_norm) and trial_norm < norm:
                u = trial
                break
            step *= 0.5
        else:
            break
    final_residual = _residual(u, dx, cfg.debye_length)
    residual_inf = float(np.max(np.abs(final_residual))) * cfg.debye_length**2
    potential = cfg.thermal_voltage * u
    # Second-order one-sided derivative avoids polluting the charge estimate
    # with the first-order boundary stencil used by a naive implementation.
    derivative = (-3.0 * potential[0] + 4.0 * potential[1] - potential[2]) / (2.0 * dx)
    surface_charge = -cfg.permittivity * derivative
    return PoissonBoltzmannSolution(x, potential, surface_potential, float(surface_charge),
                                    residual_inf, iterations, converged and residual_inf < tolerance)
