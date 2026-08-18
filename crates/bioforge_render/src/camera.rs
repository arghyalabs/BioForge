//! 3D Orbit camera and transformation matrix mathematics.

/// 3D Orbit camera supporting interactive rotation, panning, and zooming.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Camera {
    pub eye: [f32; 3],
    pub target: [f32; 3],
    pub up: [f32; 3],
    pub fov_deg: f32,
    pub aspect_ratio: f32,
    pub z_near: f32,
    pub z_far: f32,
}

impl Default for Camera {
    fn default() -> Self {
        Self {
            eye: [0.0, 0.0, 30.0],
            target: [0.0, 0.0, 0.0],
            up: [0.0, 1.0, 0.0],
            fov_deg: 45.0,
            aspect_ratio: 16.0 / 9.0,
            z_near: 0.1,
            z_far: 1000.0,
        }
    }
}

impl Camera {
    /// Create a camera with custom view parameters.
    #[must_use]
    pub fn new(eye: [f32; 3], target: [f32; 3], aspect_ratio: f32) -> Self {
        Self {
            eye,
            target,
            aspect_ratio,
            ..Default::default()
        }
    }

    /// Compute 4x4 Look-At View Matrix.
    #[must_use]
    pub fn view_matrix(&self) -> [[f32; 4]; 4] {
        let f = normalize([
            self.target[0] - self.eye[0],
            self.target[1] - self.eye[1],
            self.target[2] - self.eye[2],
        ]);
        let s = normalize(cross(f, self.up));
        let u = cross(s, f);

        [
            [s[0], u[0], -f[0], 0.0],
            [s[1], u[1], -f[1], 0.0],
            [s[2], u[2], -f[2], 0.0],
            [
                -dot(s, self.eye),
                -dot(u, self.eye),
                dot(f, self.eye),
                1.0,
            ],
        ]
    }

    /// Compute 4x4 Perspective Projection Matrix.
    #[must_use]
    pub fn projection_matrix(&self) -> [[f32; 4]; 4] {
        let fov_rad = self.fov_deg.to_radians();
        let f = 1.0 / (fov_rad * 0.5).tan();
        let ar = self.aspect_ratio;
        let zn = self.z_near;
        let zf = self.z_far;

        [
            [f / ar, 0.0, 0.0, 0.0],
            [0.0, f, 0.0, 0.0],
            [0.0, 0.0, (zf + zn) / (zn - zf), -1.0],
            [0.0, 0.0, (2.0 * zf * zn) / (zn - zf), 0.0],
        ]
    }

    /// Zoom camera relative to target by multiplying distance by `factor`.
    pub fn zoom(&mut self, factor: f32) {
        let d = [
            self.eye[0] - self.target[0],
            self.eye[1] - self.target[1],
            self.eye[2] - self.target[2],
        ];
        let clamped_factor = factor.clamp(0.1, 10.0);
        self.eye = [
            self.target[0] + d[0] * clamped_factor,
            self.target[1] + d[1] * clamped_factor,
            self.target[2] + d[2] * clamped_factor,
        ];
    }
}

fn cross(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn normalize(v: [f32; 3]) -> [f32; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-12);
    [v[0] / len, v[1] / len, v[2] / len]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_camera_view_matrix_orthonormality() {
        let camera = Camera::default();
        let v = camera.view_matrix();

        // 3x3 rotational part must be orthonormal
        let r0 = [v[0][0], v[1][0], v[2][0]];
        let r1 = [v[0][1], v[1][1], v[2][1]];

        let dot_product = dot(r0, r1);
        assert!(dot_product.abs() < 1e-6);
    }
}
