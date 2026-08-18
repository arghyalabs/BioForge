//! Procedural 3D mesh generation (Icospheres, Cylinders, and Triangle Buffers).

use crate::color::Color;

/// A single vertex with 3D position, normal vector, and RGBA color.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vertex {
    pub position: [f32; 3],
    pub normal: [f32; 3],
    pub color: [f32; 4],
}

/// A 3D indexed triangle mesh.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Mesh {
    pub vertices: Vec<Vertex>,
    pub indices: Vec<u32>,
}

impl Mesh {
    /// Create a new empty mesh.
    #[must_use]
    pub fn empty() -> Self {
        Self {
            vertices: Vec::new(),
            indices: Vec::new(),
        }
    }

    /// Number of vertices in the mesh.
    #[must_use]
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Number of triangles in the mesh.
    #[must_use]
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }

    /// Append another mesh into this mesh, offsetting indices automatically.
    pub fn append(&mut self, other: &Mesh) {
        let offset = self.vertices.len() as u32;
        self.vertices.extend_from_slice(&other.vertices);
        for &idx in &other.indices {
            self.indices.push(offset + idx);
        }
    }

    /// Translate all vertices in this mesh by `[dx, dy, dz]`.
    pub fn translate(&mut self, offset: [f32; 3]) {
        for v in &mut self.vertices {
            v.position[0] += offset[0];
            v.position[1] += offset[1];
            v.position[2] += offset[2];
        }
    }
}

/// Generate a 3D sphere mesh at `center` with `radius` using icosphere subdivision.
#[must_use]
pub fn generate_sphere(
    center: [f64; 3],
    radius: f64,
    subdivisions: u32,
    color: Color,
) -> Mesh {
    let r = radius as f32;
    let cx = center[0] as f32;
    let cy = center[1] as f32;
    let cz = center[2] as f32;
    let col = color.to_array();

    // 12 base vertices of an icosahedron
    let t = (1.0 + 5.0_f32.sqrt()) / 2.0;

    let base_verts = [
        [-1.0,  t, 0.0], [ 1.0,  t, 0.0], [-1.0, -t, 0.0], [ 1.0, -t, 0.0],
        [0.0, -1.0,  t], [0.0,  1.0,  t], [0.0, -1.0, -t], [0.0,  1.0, -t],
        [ t, 0.0, -1.0], [ t, 0.0,  1.0], [-t, 0.0, -1.0], [-t, 0.0,  1.0],
    ];

    let mut vertices: Vec<[f32; 3]> = base_verts
        .iter()
        .map(|v| {
            let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            [v[0] / len, v[1] / len, v[2] / len]
        })
        .collect();

    let mut faces: Vec<[u32; 3]> = vec![
        [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
        [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
        [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
        [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
    ];

    // Subdivide faces
    for _ in 0..subdivisions.min(3) {
        let mut new_faces = Vec::with_capacity(faces.len() * 4);
        let mut midpoint_cache = std::collections::HashMap::new();

        for face in faces {
            let a = face[0];
            let b = face[1];
            let c = face[2];

            let ab = get_midpoint(a, b, &mut vertices, &mut midpoint_cache);
            let bc = get_midpoint(b, c, &mut vertices, &mut midpoint_cache);
            let ca = get_midpoint(c, a, &mut vertices, &mut midpoint_cache);

            new_faces.push([a, ab, ca]);
            new_faces.push([b, bc, ab]);
            new_faces.push([c, ca, bc]);
            new_faces.push([ab, bc, ca]);
        }
        faces = new_faces;
    }

    let final_vertices = vertices
        .into_iter()
        .map(|norm| Vertex {
            position: [cx + norm[0] * r, cy + norm[1] * r, cz + norm[2] * r],
            normal: norm,
            color: col,
        })
        .collect();

    let mut indices = Vec::with_capacity(faces.len() * 3);
    for f in faces {
        indices.push(f[0]);
        indices.push(f[1]);
        indices.push(f[2]);
    }

    Mesh {
        vertices: final_vertices,
        indices,
    }
}

fn get_midpoint(
    p1: u32,
    p2: u32,
    vertices: &mut Vec<[f32; 3]>,
    cache: &mut std::collections::HashMap<(u32, u32), u32>,
) -> u32 {
    let key = if p1 < p2 { (p1, p2) } else { (p2, p1) };
    if let Some(&idx) = cache.get(&key) {
        return idx;
    }

    let v1 = vertices[p1 as usize];
    let v2 = vertices[p2 as usize];
    let mid = [
        (v1[0] + v2[0]) * 0.5,
        (v1[1] + v2[1]) * 0.5,
        (v1[2] + v2[2]) * 0.5,
    ];
    let len = (mid[0] * mid[0] + mid[1] * mid[1] + mid[2] * mid[2]).sqrt();
    let norm_mid = [mid[0] / len, mid[1] / len, mid[2] / len];

    let new_idx = vertices.len() as u32;
    vertices.push(norm_mid);
    cache.insert(key, new_idx);
    new_idx
}

/// Generate a 3D cylinder connecting `start` to `end` with specified `radius`.
#[must_use]
pub fn generate_cylinder(
    start: [f64; 3],
    end: [f64; 3],
    radius: f64,
    segments: u32,
    color: Color,
) -> Mesh {
    let p1 = [start[0] as f32, start[1] as f32, start[2] as f32];
    let p2 = [end[0] as f32, end[1] as f32, end[2] as f32];
    let r = radius as f32;
    let n_seg = segments.max(6);
    let col = color.to_array();

    let d = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
    let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    if len <= 1e-6 {
        return Mesh::empty();
    }
    let dir = [d[0] / len, d[1] / len, d[2] / len];

    // Find two orthogonal unit vectors (u, v) perpendicular to dir
    let arbitrary = if dir[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let u = normalize(cross(dir, arbitrary));
    let v = normalize(cross(dir, u));

    let mut vertices = Vec::new();
    let mut indices = Vec::new();

    for i in 0..n_seg {
        let angle = (i as f32) * 2.0 * std::f32::consts::PI / (n_seg as f32);
        let cos_a = angle.cos();
        let sin_a = angle.sin();

        let norm = [
            cos_a * u[0] + sin_a * v[0],
            cos_a * u[1] + sin_a * v[1],
            cos_a * u[2] + sin_a * v[2],
        ];

        // Bottom ring vertex
        vertices.push(Vertex {
            position: [
                p1[0] + norm[0] * r,
                p1[1] + norm[1] * r,
                p1[2] + norm[2] * r,
            ],
            normal: norm,
            color: col,
        });

        // Top ring vertex
        vertices.push(Vertex {
            position: [
                p2[0] + norm[0] * r,
                p2[1] + norm[1] * r,
                p2[2] + norm[2] * r,
            ],
            normal: norm,
            color: col,
        });
    }

    for i in 0..n_seg {
        let next = (i + 1) % n_seg;
        let b1 = i * 2;
        let t1 = i * 2 + 1;
        let b2 = next * 2;
        let t2 = next * 2 + 1;

        indices.push(b1);
        indices.push(b2);
        indices.push(t1);

        indices.push(t1);
        indices.push(b2);
        indices.push(t2);
    }

    Mesh { vertices, indices }
}

/// Generate a split-colored cylinder (half `color1`, half `color2`).
#[must_use]
pub fn generate_split_cylinder(
    start: [f64; 3],
    end: [f64; 3],
    radius: f64,
    segments: u32,
    color1: Color,
    color2: Color,
) -> Mesh {
    let mid = [
        (start[0] + end[0]) * 0.5,
        (start[1] + end[1]) * 0.5,
        (start[2] + end[2]) * 0.5,
    ];

    let mut mesh1 = generate_cylinder(start, mid, radius, segments, color1);
    let mesh2 = generate_cylinder(mid, end, radius, segments, color2);
    mesh1.append(&mesh2);
    mesh1
}

fn cross(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalize(v: [f32; 3]) -> [f32; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-12);
    [v[0] / len, v[1] / len, v[2] / len]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sphere_mesh_generation() {
        let sphere = generate_sphere([0.0, 0.0, 0.0], 1.5, 1, Color::RED);
        assert!(sphere.vertex_count() > 12);
        assert!(sphere.triangle_count() > 20);

        // Check radius
        let p = sphere.vertices[0].position;
        let dist = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        assert!((dist - 1.5).abs() < 1e-4);
    }

    #[test]
    fn test_cylinder_mesh_generation() {
        let cyl = generate_cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 5.0], 0.2, 8, Color::BLUE);
        assert_eq!(cyl.vertex_count(), 16);
        assert_eq!(cyl.triangle_count(), 16);
    }
}
