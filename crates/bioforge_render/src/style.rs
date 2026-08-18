//! Molecular graphical representation styles.

/// Visual representation styles for biological structures.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RenderStyle {
    /// Space-filling van der Waals representation.
    SpaceFilling {
        /// Sphere mesh subdivision level (0-3).
        subdivisions: u32,
    },
    /// Ball-and-stick representation with atoms and covalent bond cylinders.
    BallAndStick {
        /// Atom sphere radius in Ångströms (default $0.35\text{ \AA}$).
        atom_radius: f64,
        /// Bond cylinder radius in Ångströms (default $0.15\text{ \AA}$).
        bond_radius: f64,
        /// Sphere subdivisions.
        subdivisions: u32,
        /// Cylinder radial segments.
        bond_segments: u32,
    },
    /// Continuous $C_\alpha$ peptide backbone ribbon/tube.
    BackboneTrace {
        /// Backbone tube radius in Ångströms (default $0.30\text{ \AA}$).
        tube_radius: f64,
        /// Cylinder segments.
        segments: u32,
    },
}

impl Default for RenderStyle {
    fn default() -> Self {
        Self::ball_and_stick()
    }
}

impl RenderStyle {
    /// Create standard Ball-and-Stick style.
    #[must_use]
    pub fn ball_and_stick() -> Self {
        Self::BallAndStick {
            atom_radius: 0.35,
            bond_radius: 0.15,
            subdivisions: 1,
            bond_segments: 8,
        }
    }

    /// Create standard Space-Filling (CPK van der Waals) style.
    #[must_use]
    pub fn space_filling() -> Self {
        Self::SpaceFilling { subdivisions: 1 }
    }

    /// Create standard Backbone trace style.
    #[must_use]
    pub fn backbone_trace() -> Self {
        Self::BackboneTrace {
            tube_radius: 0.30,
            segments: 8,
        }
    }
}
