//! 3D Visual Scene extraction from SimulationState and Wavefront OBJ exporter.

use bioforge_state::SimulationState;

use crate::color::{color_by_chain, cpk_color_for_element};
use crate::mesh::{generate_cylinder, generate_sphere, generate_split_cylinder, Mesh};
use crate::style::RenderStyle;

/// A complete 3D visual scene generated from a [`SimulationState`] snapshot.
#[derive(Debug, Clone, PartialEq)]
pub struct Scene {
    /// Combined 3D indexed triangle mesh.
    pub mesh: Mesh,
    /// Number of represented atoms.
    pub atom_count: usize,
    /// Number of represented bonds.
    pub bond_count: usize,
    /// The style used to generate this scene.
    pub style: RenderStyle,
}

impl Scene {
    /// Extract a 3D visual scene from an immutable snapshot of [`SimulationState`].
    ///
    /// Per Scientific Principle 2, this function does not modify physical state.
    #[must_use]
    pub fn from_state(state: &SimulationState, style: RenderStyle) -> Self {
        let mut combined_mesh = Mesh::empty();
        let num_atoms = state.num_atoms;
        let mut bond_count = 0;

        match style {
            RenderStyle::SpaceFilling { subdivisions } => {
                for i in 0..num_atoms {
                    let pos = state.positions[i];
                    let elem = &state.elements[i];
                    let radius = elem.vdw_radius;
                    let color = cpk_color_for_element(elem);

                    let sphere = generate_sphere(pos, radius, subdivisions, color);
                    combined_mesh.append(&sphere);
                }
            }
            RenderStyle::BallAndStick {
                atom_radius,
                bond_radius,
                subdivisions,
                bond_segments,
            } => {
                // Generate atom spheres
                for i in 0..num_atoms {
                    let pos = state.positions[i];
                    let elem = &state.elements[i];
                    let color = cpk_color_for_element(elem);

                    let sphere = generate_sphere(pos, atom_radius, subdivisions, color);
                    combined_mesh.append(&sphere);
                }

                // Generate bond cylinders
                for bond in &state.bonds {
                    let i = bond.atom1;
                    let j = bond.atom2;
                    if i < num_atoms && j < num_atoms {
                        let p1 = state.positions[i];
                        let p2 = state.positions[j];
                        let col1 = cpk_color_for_element(&state.elements[i]);
                        let col2 = cpk_color_for_element(&state.elements[j]);

                        let cyl = generate_split_cylinder(
                            p1,
                            p2,
                            bond_radius,
                            bond_segments,
                            col1,
                            col2,
                        );
                        combined_mesh.append(&cyl);
                        bond_count += 1;
                    }
                }
            }
            RenderStyle::BackboneTrace {
                tube_radius,
                segments,
            } => {
                // Find all C-alpha atoms in order
                let mut ca_indices = Vec::new();
                for i in 0..num_atoms {
                    if state.atom_names[i] == "CA" {
                        ca_indices.push(i);
                    }
                }

                // Connect consecutive CA atoms in same chain
                for w in ca_indices.windows(2) {
                    let i = w[0];
                    let j = w[1];
                    let chain_i = state.chain_ids[i];
                    let chain_j = state.chain_ids[j];

                    if chain_i == chain_j {
                        let p1 = state.positions[i];
                        let p2 = state.positions[j];
                        let color = color_by_chain(chain_i);

                        let tube = generate_cylinder(p1, p2, tube_radius, segments, color);
                        combined_mesh.append(&tube);
                        bond_count += 1;
                    }
                }
            }
        }

        Self {
            mesh: combined_mesh,
            atom_count: num_atoms,
            bond_count,
            style,
        }
    }

    /// Export the 3D scene into Wavefront `.obj` format.
    ///
    /// Compatible with 3D modeling tools (Blender, MeshLab) and WebGL viewers (Three.js).
    #[must_use]
    pub fn export_obj(&self) -> String {
        let mut out = String::new();
        out.push_str("# BioForge 3D Scene Export\n");
        out.push_str(&format!(
            "# Atoms: {}, Bonds: {}, Triangles: {}\n",
            self.atom_count,
            self.bond_count,
            self.mesh.triangle_count()
        ));

        // Vertex positions (v x y z r g b)
        for v in &self.mesh.vertices {
            out.push_str(&format!(
                "v {:.4} {:.4} {:.4} {:.3} {:.3} {:.3}\n",
                v.position[0],
                v.position[1],
                v.position[2],
                v.color[0],
                v.color[1],
                v.color[2]
            ));
        }

        // Vertex normals (vn x y z)
        for v in &self.mesh.vertices {
            out.push_str(&format!(
                "vn {:.4} {:.4} {:.4}\n",
                v.normal[0], v.normal[1], v.normal[2]
            ));
        }

        // Triangle faces (f v1//vn1 v2//vn2 v3//vn3) - 1-indexed in OBJ
        for chunk in self.mesh.indices.chunks_exact(3) {
            let i1 = chunk[0] + 1;
            let i2 = chunk[1] + 1;
            let i3 = chunk[2] + 1;
            out.push_str(&format!("f {}//{} {}//{} {}//{}\n", i1, i1, i2, i2, i3, i3));
        }

        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bioforge_biology::pdb::parse_pdb;

    const TINY_PDB: &str = "\
ATOM      1  N   ALA A   1       1.458   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       1.562   2.163   1.252  1.00  0.00           C
ATOM      4  O   ALA A   1       0.735   1.685   2.056  1.00  0.00           O
ATOM      5  CB  ALA A   1       3.529   1.388   0.000  1.00  0.00           C
END
";

    #[test]
    fn test_scene_generation_space_filling() {
        let mol = parse_pdb(TINY_PDB, "ala").unwrap();
        let state = SimulationState::from_molecule(&mol, None);

        let scene = Scene::from_state(&state, RenderStyle::space_filling());
        assert_eq!(scene.atom_count, 5);
        assert!(scene.mesh.vertex_count() > 50);
        assert!(scene.mesh.triangle_count() > 50);
    }

    #[test]
    fn test_scene_generation_ball_and_stick_and_obj_export() {
        let mol = parse_pdb(TINY_PDB, "ala").unwrap();
        let state = SimulationState::from_molecule(&mol, None);

        let scene = Scene::from_state(&state, RenderStyle::ball_and_stick());
        assert_eq!(scene.atom_count, 5);

        let obj = scene.export_obj();
        assert!(obj.starts_with("# BioForge 3D Scene Export"));
        assert!(obj.contains("v "));
        assert!(obj.contains("vn "));
        assert!(obj.contains("f "));
    }
}
