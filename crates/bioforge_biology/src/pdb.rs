//! Basic PDB file parser.
//!
//! Parses ATOM and HETATM records from PDB-format files into a [`Molecule`].
//!
//! ## PDB Format Reference
//!
//! The PDB fixed-column format for ATOM records:
//! ```text
//! Columns  1– 6:  Record type ("ATOM  " or "HETATM")
//! Columns  7–11:  Atom serial number
//! Columns 13–16:  Atom name
//! Column     17:  Alternate location indicator
//! Columns 18–20:  Residue name
//! Column     22:  Chain ID
//! Columns 23–26:  Residue sequence number
//! Columns 31–38:  X coordinate (Ångströms)
//! Columns 39–46:  Y coordinate (Ångströms)
//! Columns 47–54:  Z coordinate (Ångströms)
//! Columns 77–78:  Element symbol (right-justified)
//! ```
//!
//! ## Limitations (Scientific Principle 5 — approximations visible)
//!
//! - No CONECT records (bond parsing)
//! - No multi-model support (only first model)
//! - No crystallographic symmetry (CRYST1, SCALE, etc.)
//! - No alternate conformations (uses first conformation only)
//! - No anisotropic temperature factors

use crate::atom::Atom;
use crate::element::Element;
use crate::error::BiologyError;
use crate::molecule::Molecule;

/// Parse PDB-format text content into a [`Molecule`].
///
/// # Arguments
///
/// * `content` — The full text content of a PDB file.
/// * `name` — A name for the resulting molecule.
///
/// # Errors
///
/// Returns [`BiologyError::PdbParseError`] for malformed lines.
pub fn parse_pdb(content: &str, name: &str) -> Result<Molecule, BiologyError> {
    let mut mol = Molecule::new(name);

    for (line_idx, line) in content.lines().enumerate() {
        let line_num = line_idx + 1;

        // Only process ATOM and HETATM records
        if !line.starts_with("ATOM") && !line.starts_with("HETATM") {
            // Stop at END or ENDMDL (only parse first model)
            if line.starts_with("END") {
                break;
            }
            continue;
        }

        // PDB lines must be at least 54 characters for coordinates
        if line.len() < 54 {
            return Err(BiologyError::PdbParseError {
                line: line_num,
                message: format!(
                    "line too short ({} chars, need at least 54)",
                    line.len()
                ),
            });
        }

        let atom = parse_atom_line(line, line_num)?;
        mol.atoms.push(atom);
    }

    Ok(mol)
}

/// Parse a single ATOM/HETATM line into an [`Atom`].
fn parse_atom_line(line: &str, line_num: usize) -> Result<Atom, BiologyError> {
    // Parse atom serial number (columns 7-11, 0-indexed 6-11)
    let serial_str = safe_substr(line, 6, 11).trim();
    let serial: u32 = serial_str.parse().map_err(|_| BiologyError::PdbParseError {
        line: line_num,
        message: format!("invalid atom serial number: '{}'", serial_str),
    })?;

    // Parse atom name (columns 13-16, 0-indexed 12-16)
    let atom_name = safe_substr(line, 12, 16).trim().to_string();

    // Parse residue name (columns 18-20, 0-indexed 17-20)
    let residue_name = safe_substr(line, 17, 20).trim().to_string();

    // Parse chain ID (column 22, 0-indexed 21)
    let chain_id = line.as_bytes().get(21).map(|&b| b as char).filter(|c| !c.is_whitespace());

    // Parse residue sequence number (columns 23-26, 0-indexed 22-26)
    let resid_str = safe_substr(line, 22, 26).trim();
    let residue_id: Option<i32> = if resid_str.is_empty() {
        None
    } else {
        Some(resid_str.parse().map_err(|_| BiologyError::PdbParseError {
            line: line_num,
            message: format!("invalid residue number: '{}'", resid_str),
        })?)
    };

    // Parse coordinates (columns 31-38, 39-46, 47-54, 0-indexed 30-38, 38-46, 46-54)
    let x = parse_coord(line, 30, 38, "X", line_num)?;
    let y = parse_coord(line, 38, 46, "Y", line_num)?;
    let z = parse_coord(line, 46, 54, "Z", line_num)?;

    // Parse element symbol (columns 77-78, 0-indexed 76-78)
    // Fall back to guessing from atom name if element columns are absent
    let element = resolve_element(line, &atom_name, line_num)?;

    let mut atom = Atom::new(serial, element, [x, y, z], atom_name);
    atom.residue_name = if residue_name.is_empty() {
        None
    } else {
        Some(residue_name)
    };
    atom.residue_id = residue_id;
    atom.chain_id = chain_id;

    Ok(atom)
}

/// Parse a floating-point coordinate from fixed columns.
fn parse_coord(
    line: &str,
    start: usize,
    end: usize,
    axis: &str,
    line_num: usize,
) -> Result<f64, BiologyError> {
    let s = safe_substr(line, start, end).trim();
    s.parse().map_err(|_| BiologyError::PdbParseError {
        line: line_num,
        message: format!("invalid {} coordinate: '{}'", axis, s),
    })
}

/// Resolve the element from PDB columns 77-78 or guess from atom name.
fn resolve_element(line: &str, atom_name: &str, line_num: usize) -> Result<Element, BiologyError> {
    // Try columns 77-78 first (if line is long enough)
    if line.len() >= 78 {
        let elem_str = safe_substr(line, 76, 78).trim();
        if !elem_str.is_empty() {
            if let Some(elem) = Element::from_symbol(elem_str) {
                return Ok(elem);
            }
            // Try single-char variant
            if elem_str.len() > 1 {
                if let Some(elem) = Element::from_symbol(&elem_str[..1]) {
                    return Ok(elem);
                }
            }
        }
    }

    // Guess from atom name: strip digits and spaces, take first letter(s)
    // PDB atom names like "CA" = Carbon-alpha, "FE" = Iron
    // Strategy: try 1-char first (C, N, O, H are much more common in biology),
    // then 2-char (Fe, Ca, Zn) only if 1-char fails.
    let guess = atom_name
        .chars()
        .filter(|c| c.is_alphabetic())
        .collect::<String>();

    if !guess.is_empty() {
        // Try 1-char symbol first (C, N, O, H, S, P — most common)
        let one_char: String = guess.chars().next().unwrap().to_uppercase().collect();
        if let Some(elem) = Element::from_symbol(&one_char) {
            return Ok(elem);
        }

        // Try 2-char symbol (Fe, Ca, Zn, etc.)
        if guess.len() >= 2 {
            let two_char = format!(
                "{}{}",
                guess.chars().next().unwrap().to_uppercase(),
                guess.chars().nth(1).unwrap().to_lowercase()
            );
            if let Some(elem) = Element::from_symbol(&two_char) {
                return Ok(elem);
            }
        }
    }

    Err(BiologyError::PdbParseError {
        line: line_num,
        message: format!("could not determine element for atom '{}'", atom_name),
    })
}

/// Safely extract a substring by byte indices, padding with spaces if needed.
fn safe_substr(s: &str, start: usize, end: usize) -> &str {
    let end = end.min(s.len());
    if start >= s.len() {
        return "";
    }
    &s[start..end]
}

/// Parse a PDB file from disk.
///
/// # Errors
///
/// Returns [`BiologyError::IoError`] if the file cannot be read,
/// or [`BiologyError::PdbParseError`] for malformed content.
pub fn parse_pdb_file(path: &str) -> Result<Molecule, BiologyError> {
    let content = std::fs::read_to_string(path)?;
    // Use the filename (without extension) as the molecule name
    let name = std::path::Path::new(path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("molecule");
    parse_pdb(&content, name)
}

// ─── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimal PDB content: 3 atoms of an alanine residue.
    const TINY_PDB: &str = "\
ATOM      1  N   ALA A   1       1.458   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       1.562   2.163   1.252  1.00  0.00           C
END
";

    #[test]
    fn test_parse_tiny_pdb() {
        let mol = parse_pdb(TINY_PDB, "alanine").unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.name, "alanine");
    }

    #[test]
    fn test_atom_properties() {
        let mol = parse_pdb(TINY_PDB, "ala").unwrap();

        let n = &mol.atoms[0];
        assert_eq!(n.id, 1);
        assert_eq!(n.name, "N");
        assert_eq!(n.element.symbol, "N");
        assert_eq!(n.residue_name.as_deref(), Some("ALA"));
        assert_eq!(n.residue_id, Some(1));
        assert_eq!(n.chain_id, Some('A'));
        assert!((n.position[0] - 1.458).abs() < 0.001);
    }

    #[test]
    fn test_coordinates() {
        let mol = parse_pdb(TINY_PDB, "ala").unwrap();
        let ca = &mol.atoms[1];
        assert!((ca.position[0] - 2.009).abs() < 0.001);
        assert!((ca.position[1] - 1.420).abs() < 0.001);
        assert!((ca.position[2] - 0.000).abs() < 0.001);
    }

    #[test]
    fn test_element_from_columns_77_78() {
        let mol = parse_pdb(TINY_PDB, "test").unwrap();
        assert_eq!(mol.atoms[0].element.symbol, "N");
        assert_eq!(mol.atoms[1].element.symbol, "C");
        assert_eq!(mol.atoms[2].element.symbol, "C");
    }

    #[test]
    fn test_hetatm_record() {
        let pdb = "\
HETATM    1  O   HOH A 100       5.000   6.000   7.000  1.00  0.00           O
END
";
        let mol = parse_pdb(pdb, "water").unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(mol.atoms[0].element.symbol, "O");
        assert_eq!(mol.atoms[0].residue_name.as_deref(), Some("HOH"));
    }

    #[test]
    fn test_empty_pdb() {
        let mol = parse_pdb("END\n", "empty").unwrap();
        assert_eq!(mol.atom_count(), 0);
    }

    #[test]
    fn test_line_too_short() {
        let bad = "ATOM      1  N\n";
        let result = parse_pdb(bad, "bad");
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(matches!(err, BiologyError::PdbParseError { .. }));
    }

    #[test]
    fn test_stops_at_end() {
        let pdb = "\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
END
ATOM      2  CA  ALA A   1       4.000   5.000   6.000  1.00  0.00           C
";
        let mol = parse_pdb(pdb, "test").unwrap();
        // Should only get 1 atom (stops at END)
        assert_eq!(mol.atom_count(), 1);
    }

    #[test]
    fn test_total_mass_from_pdb() {
        let mol = parse_pdb(TINY_PDB, "ala").unwrap();
        let mass = mol.total_mass();
        // N(14.007) + C(12.011) + C(12.011) = 38.029
        assert!((mass - 38.029).abs() < 0.01, "got {}", mass);
    }

    #[test]
    fn test_element_guessing_from_name() {
        // PDB without element columns (line < 78 chars)
        let pdb = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00\nEND\n";
        let mol = parse_pdb(pdb, "test").unwrap();
        // "CA" atom name → should guess "C" element (Carbon alpha)
        assert_eq!(mol.atoms[0].element.symbol, "C");
    }
}
