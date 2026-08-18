//! Color representations and standard biophysical CPK color palette.

use bioforge_biology::Element;

/// RGBA color representation with normalized floating-point components $[0.0, 1.0]$.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Color {
    pub r: f32,
    pub g: f32,
    pub b: f32,
    pub a: f32,
}

impl Color {
    pub const WHITE: Color = Color::rgba(0.95, 0.95, 0.95, 1.0);
    pub const BLACK: Color = Color::rgba(0.05, 0.05, 0.05, 1.0);
    pub const GREY: Color = Color::rgba(0.35, 0.35, 0.35, 1.0);
    pub const RED: Color = Color::rgba(0.90, 0.15, 0.15, 1.0);
    pub const GREEN: Color = Color::rgba(0.15, 0.85, 0.25, 1.0);
    pub const BLUE: Color = Color::rgba(0.15, 0.35, 0.90, 1.0);
    pub const YELLOW: Color = Color::rgba(0.95, 0.85, 0.10, 1.0);
    pub const ORANGE: Color = Color::rgba(1.00, 0.55, 0.05, 1.0);
    pub const PURPLE: Color = Color::rgba(0.65, 0.25, 0.85, 1.0);
    pub const CYAN: Color = Color::rgba(0.10, 0.85, 0.90, 1.0);

    /// Create an RGBA color.
    #[must_use]
    pub const fn rgba(r: f32, g: f32, b: f32, a: f32) -> Self {
        Self { r, g, b, a }
    }

    /// Create an opaque RGB color with alpha = 1.0.
    #[must_use]
    pub const fn rgb(r: f32, g: f32, b: f32) -> Self {
        Self { r, g, b, a: 1.0 }
    }

    /// Convert to `[f32; 4]` array.
    #[must_use]
    pub fn to_array(&self) -> [f32; 4] {
        [self.r, self.g, self.b, self.a]
    }

    /// Convert to `[f32; 3]` RGB array.
    #[must_use]
    pub fn to_rgb_array(&self) -> [f32; 3] {
        [self.r, self.g, self.b]
    }

    /// Linear interpolation between two colors.
    #[must_use]
    pub fn lerp(&self, other: &Color, t: f32) -> Color {
        let t_clamped = t.clamp(0.0, 1.0);
        Color {
            r: self.r + (other.r - self.r) * t_clamped,
            g: self.g + (other.g - self.g) * t_clamped,
            b: self.b + (other.b - self.b) * t_clamped,
            a: self.a + (other.a - self.a) * t_clamped,
        }
    }
}

/// Standard Corey-Pauling-Koltun (CPK) elemental coloring for biological macromolecules.
#[must_use]
pub fn cpk_color_for_element(elem: &Element) -> Color {
    match elem.symbol {
        "H" => Color::WHITE,
        "C" => Color::GREY,
        "N" => Color::BLUE,
        "O" => Color::RED,
        "S" => Color::YELLOW,
        "P" => Color::ORANGE,
        "F" | "Cl" => Color::GREEN,
        "Br" => Color::rgba(0.60, 0.15, 0.15, 1.0),
        "I" => Color::PURPLE,
        "Na" | "K" => Color::rgba(0.45, 0.20, 0.80, 1.0),
        "Ca" | "Mg" => Color::rgba(0.20, 0.65, 0.20, 1.0),
        "Fe" => Color::rgba(0.85, 0.45, 0.10, 1.0),
        "Zn" | "Cu" => Color::rgba(0.60, 0.60, 0.75, 1.0),
        _ => Color::CYAN,
    }
}

/// Distinct color generator for protein chains ('A', 'B', etc.).
#[must_use]
pub fn color_by_chain(chain_id: Option<char>) -> Color {
    match chain_id {
        Some('A') => Color::rgba(0.20, 0.60, 0.85, 1.0), // Soft Blue
        Some('B') => Color::rgba(0.85, 0.35, 0.25, 1.0), // Coral Red
        Some('C') => Color::rgba(0.30, 0.75, 0.40, 1.0), // Emerald
        Some('D') => Color::rgba(0.90, 0.70, 0.20, 1.0), // Gold
        Some('E') => Color::rgba(0.70, 0.35, 0.80, 1.0), // Violet
        _ => Color::rgba(0.50, 0.50, 0.50, 1.0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cpk_colors() {
        let c = Element::from_symbol("C").unwrap();
        let o = Element::from_symbol("O").unwrap();
        let n = Element::from_symbol("N").unwrap();
        let h = Element::from_symbol("H").unwrap();

        assert_eq!(cpk_color_for_element(&c), Color::GREY);
        assert_eq!(cpk_color_for_element(&o), Color::RED);
        assert_eq!(cpk_color_for_element(&n), Color::BLUE);
        assert_eq!(cpk_color_for_element(&h), Color::WHITE);
    }

    #[test]
    fn test_color_lerp() {
        let white = Color::WHITE;
        let black = Color::BLACK;
        let mid = white.lerp(&black, 0.5);

        assert!((mid.r - 0.5).abs() < 0.05);
        assert!((mid.g - 0.5).abs() < 0.05);
        assert!((mid.b - 0.5).abs() < 0.05);
    }
}
