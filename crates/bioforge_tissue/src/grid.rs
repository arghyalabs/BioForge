//! Spatial grids, boundary conditions, and finite-difference Laplacian operators in 1D and 2D.

use serde::{Deserialize, Serialize};

use crate::error::TissueError;

/// Spatial boundary condition for reaction-diffusion PDEs.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum BoundaryCondition {
    /// Zero-flux / impermeable boundary ($\frac{\partial C}{\partial n} = 0$).
    ZeroFluxNeumann,
    /// Fixed boundary concentration value ($C = C_{\text{fixed}}$).
    FixedDirichlet(f64),
    /// Periodic / toroidal boundary condition ($C_{-1} = C_{N-1}, C_{N} = C_0$).
    Periodic,
}

/// 1-Dimensional spatial discretization grid.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Grid1D {
    /// Number of discrete spatial points $N$.
    pub num_points: usize,
    /// Spatial grid spacing $\Delta x$ in micrometers ($\mu\text{m}$).
    pub dx_um: f64,
    /// Scalar field values across grid points.
    pub values: Vec<f64>,
    /// Left boundary condition ($x = 0$).
    pub boundary_left: BoundaryCondition,
    /// Right boundary condition ($x = L$).
    pub boundary_right: BoundaryCondition,
}

impl Grid1D {
    /// Create a new 1D spatial grid with homogeneous initial value.
    pub fn new(num_points: usize, dx_um: f64, initial_val: f64) -> Result<Self, TissueError> {
        if num_points < 2 {
            return Err(TissueError::InvalidGridDimensions {
                width: num_points,
                height: 1,
            });
        }
        Ok(Self {
            num_points,
            dx_um: dx_um.max(1e-6),
            values: vec![initial_val; num_points],
            boundary_left: BoundaryCondition::ZeroFluxNeumann,
            boundary_right: BoundaryCondition::ZeroFluxNeumann,
        })
    }

    /// Total physical length of the 1D spatial domain in $\mu\text{m}$.
    #[must_use]
    pub fn total_length_um(&self) -> f64 {
        (self.num_points - 1) as f64 * self.dx_um
    }

    /// Compute the 1D finite-difference Laplacian $\nabla^2 C \approx \frac{C_{i+1} - 2C_i + C_{i-1}}{\Delta x^2}$ in $\mu\text{m}^{-2}$.
    #[must_use]
    pub fn compute_laplacian(&self) -> Vec<f64> {
        let n = self.num_points;
        let mut lap = vec![0.0; n];
        let dx2 = self.dx_um * self.dx_um;

        for i in 0..n {
            let left = if i == 0 {
                match self.boundary_left {
                    BoundaryCondition::ZeroFluxNeumann => self.values[0],
                    BoundaryCondition::FixedDirichlet(val) => val,
                    BoundaryCondition::Periodic => self.values[n - 1],
                }
            } else {
                self.values[i - 1]
            };

            let right = if i == n - 1 {
                match self.boundary_right {
                    BoundaryCondition::ZeroFluxNeumann => self.values[n - 1],
                    BoundaryCondition::FixedDirichlet(val) => val,
                    BoundaryCondition::Periodic => self.values[0],
                }
            } else {
                self.values[i + 1]
            };

            lap[i] = (left - 2.0 * self.values[i] + right) / dx2;
        }

        lap
    }
}

/// 2-Dimensional spatial discretization grid.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Grid2D {
    /// Number of grid columns (X axis).
    pub width: usize,
    /// Number of grid rows (Y axis).
    pub height: usize,
    /// Spatial grid spacing $\Delta x = \Delta y$ in micrometers ($\mu\text{m}$).
    pub dx_um: f64,
    /// Row-major flattened scalar field values: `values[y * width + x]`.
    pub values: Vec<f64>,
    /// Boundary condition applied to all edges.
    pub boundary: BoundaryCondition,
}

impl Grid2D {
    /// Create a new 2D spatial grid with homogeneous initial value.
    pub fn new(
        width: usize,
        height: usize,
        dx_um: f64,
        initial_val: f64,
    ) -> Result<Self, TissueError> {
        if width < 2 || height < 2 {
            return Err(TissueError::InvalidGridDimensions { width, height });
        }
        Ok(Self {
            width,
            height,
            dx_um: dx_um.max(1e-6),
            values: vec![initial_val; width * height],
            boundary: BoundaryCondition::ZeroFluxNeumann,
        })
    }

    /// Get value at coordinate `(x, y)`.
    #[must_use]
    pub fn get(&self, x: usize, y: usize) -> f64 {
        if x < self.width && y < self.height {
            self.values[y * self.width + x]
        } else {
            0.0
        }
    }

    /// Set value at coordinate `(x, y)`.
    pub fn set(&mut self, x: usize, y: usize, val: f64) {
        if x < self.width && y < self.height {
            self.values[y * self.width + x] = val;
        }
    }

    /// Compute 2D finite-difference Laplacian $\nabla^2 C = \frac{\partial^2 C}{\partial x^2} + \frac{\partial^2 C}{\partial y^2}$ using 5-point stencil.
    #[must_use]
    pub fn compute_laplacian(&self) -> Vec<f64> {
        let w = self.width;
        let h = self.height;
        let mut lap = vec![0.0; w * h];
        let dx2 = self.dx_um * self.dx_um;

        for y in 0..h {
            for x in 0..w {
                let center = self.get(x, y);

                let left = if x == 0 {
                    match self.boundary {
                        BoundaryCondition::ZeroFluxNeumann => center,
                        BoundaryCondition::FixedDirichlet(v) => v,
                        BoundaryCondition::Periodic => self.get(w - 1, y),
                    }
                } else {
                    self.get(x - 1, y)
                };

                let right = if x == w - 1 {
                    match self.boundary {
                        BoundaryCondition::ZeroFluxNeumann => center,
                        BoundaryCondition::FixedDirichlet(v) => v,
                        BoundaryCondition::Periodic => self.get(0, y),
                    }
                } else {
                    self.get(x + 1, y)
                };

                let top = if y == 0 {
                    match self.boundary {
                        BoundaryCondition::ZeroFluxNeumann => center,
                        BoundaryCondition::FixedDirichlet(v) => v,
                        BoundaryCondition::Periodic => self.get(x, h - 1),
                    }
                } else {
                    self.get(x, y - 1)
                };

                let bottom = if y == h - 1 {
                    match self.boundary {
                        BoundaryCondition::ZeroFluxNeumann => center,
                        BoundaryCondition::FixedDirichlet(v) => v,
                        BoundaryCondition::Periodic => self.get(x, 0),
                    }
                } else {
                    self.get(x, y + 1)
                };

                lap[y * w + x] = (left + right + top + bottom - 4.0 * center) / dx2;
            }
        }

        lap
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_grid1d_constant_field_has_zero_laplacian() {
        let grid = Grid1D::new(10, 1.0, 5.0).unwrap();
        let lap = grid.compute_laplacian();
        for &l in &lap {
            assert!(l.abs() < 1e-12);
        }
    }

    #[test]
    fn test_grid1d_parabolic_field_has_constant_laplacian() {
        // u(x) = x^2 => d^2 u / dx^2 = 2.0
        let mut grid = Grid1D::new(5, 1.0, 0.0).unwrap();
        for i in 0..5 {
            grid.values[i] = (i as f64).powi(2); // [0, 1, 4, 9, 16]
        }
        let lap = grid.compute_laplacian();
        // Interior points (i = 1, 2, 3) must equal 2.0 exactly
        for i in 1..4 {
            assert!((lap[i] - 2.0).abs() < 1e-12, "at i={}, lap={}", i, lap[i]);
        }
    }

    #[test]
    fn test_grid2d_laplacian_stencil() {
        let mut grid = Grid2D::new(5, 5, 1.0, 0.0).unwrap();
        grid.set(2, 2, 10.0); // Point source in center

        let lap = grid.compute_laplacian();
        // Center point: (0 + 0 + 0 + 0 - 4*10) / 1^2 = -40.0
        assert_eq!(lap[2 * 5 + 2], -40.0);
        // Neighbors at (1, 2), (3, 2), (2, 1), (2, 3): (10 - 0) / 1^2 = +10.0
        assert_eq!(lap[2 * 5 + 1], 10.0);
        assert_eq!(lap[2 * 5 + 3], 10.0);
    }
}
