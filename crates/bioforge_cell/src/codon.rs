//! Universal 64-codon genetic code translation table, DNA transcription, and reverse complementation.

use crate::error::CellError;

/// Translate a 3-letter RNA triplet codon into a single-letter amino acid code or Stop (`*`).
#[must_use]
pub fn translate_rna_codon(codon: &str) -> Option<char> {
    if codon.len() != 3 {
        return None;
    }
    let c = codon.to_ascii_uppercase();
    match c.as_str() {
        // Phenylalanine & Leucine
        "UUU" | "UUC" => Some('F'),
        "UUA" | "UUG" | "CUU" | "CUC" | "CUA" | "CUG" => Some('L'),
        // Isoleucine & Methionine (Start)
        "AUU" | "AUC" | "AUA" => Some('I'),
        "AUG" => Some('M'),
        // Valine
        "GUU" | "GUC" | "GUA" | "GUG" => Some('V'),
        // Serine
        "UCU" | "UCC" | "UCA" | "UCG" | "AGU" | "AGC" => Some('S'),
        // Proline
        "CCU" | "CCC" | "CCA" | "CCG" => Some('P'),
        // Threonine
        "ACU" | "ACC" | "ACA" | "ACG" => Some('T'),
        // Alanine
        "GCU" | "GCC" | "GCA" | "GCG" => Some('A'),
        // Tyrosine & Stop
        "UAU" | "UAC" => Some('Y'),
        "UAA" | "UAG" | "UGA" => Some('*'),
        // Histidine & Glutamine
        "CAU" | "CAC" => Some('H'),
        "CAA" | "CAG" => Some('Q'),
        // Asparagine & Lysine
        "AAU" | "AAC" => Some('N'),
        "AAA" | "AAG" => Some('K'),
        // Aspartate & Glutamate
        "GAU" | "GAC" => Some('D'),
        "GAA" | "GAG" => Some('E'),
        // Cysteine & Tryptophan
        "UGU" | "UGC" => Some('C'),
        "UGG" => Some('W'),
        // Arginine
        "CGU" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => Some('R'),
        // Glycine
        "GGU" | "GGC" | "GGA" | "GGG" => Some('G'),
        _ => None,
    }
}

/// Transcribe DNA coding strand to mRNA transcript (replaces Thymine 'T' with Uracil 'U').
pub fn transcribe_dna_to_rna(dna: &str) -> Result<String, CellError> {
    let mut rna = String::with_capacity(dna.len());
    for (i, ch) in dna.chars().enumerate() {
        match ch.to_ascii_uppercase() {
            'A' => rna.push('A'),
            'C' => rna.push('C'),
            'G' => rna.push('G'),
            'T' | 'U' => rna.push('U'),
            invalid => {
                return Err(CellError::InvalidNucleotide {
                    nucleotide: invalid,
                    position: i,
                })
            }
        }
    }
    Ok(rna)
}

/// Reverse complement of a DNA sequence ($5' \to 3'$ antisense strand).
pub fn reverse_complement_dna(dna: &str) -> Result<String, CellError> {
    let mut rev_comp = String::with_capacity(dna.len());
    for (i, ch) in dna.chars().rev().enumerate() {
        let comp = match ch.to_ascii_uppercase() {
            'A' => 'T',
            'T' | 'U' => 'A',
            'C' => 'G',
            'G' => 'C',
            invalid => {
                return Err(CellError::InvalidNucleotide {
                    nucleotide: invalid,
                    position: dna.len() - 1 - i,
                })
            }
        };
        rev_comp.push(comp);
    }
    Ok(rev_comp)
}

/// Translate an mRNA sequence into a peptide protein string using the standard universal genetic code.
pub fn translate_rna_to_protein(rna: &str) -> Result<String, CellError> {
    let clean_rna: String = rna.chars().filter(|c| !c.is_whitespace()).collect();
    if clean_rna.len() % 3 != 0 {
        return Err(CellError::IncompleteCodonSequence {
            length: clean_rna.len(),
        });
    }

    let mut protein = String::with_capacity(clean_rna.len() / 3);
    for i in (0..clean_rna.len()).step_by(3) {
        let codon = &clean_rna[i..i + 3];
        match translate_rna_codon(codon) {
            Some('*') => {
                protein.push('*');
                break; // Stop codon terminates translation
            }
            Some(aa) => protein.push(aa),
            None => {
                let first_invalid = codon.chars().next().unwrap_or('?');
                return Err(CellError::InvalidNucleotide {
                    nucleotide: first_invalid,
                    position: i,
                });
            }
        }
    }

    Ok(protein)
}

/// Translate a coding DNA sequence directly into a peptide protein string.
pub fn translate_dna_to_protein(dna: &str) -> Result<String, CellError> {
    let rna = transcribe_dna_to_rna(dna)?;
    translate_rna_to_protein(&rna)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_transcription_and_reverse_complement() {
        let dna = "ATGCGATCGTAA";
        let rna = transcribe_dna_to_rna(dna).unwrap();
        assert_eq!(rna, "AUGCGAUCGUAA");

        let rev = reverse_complement_dna(dna).unwrap();
        assert_eq!(rev, "TTACGATCGCAT");
    }

    #[test]
    fn test_complete_translation_central_dogma() {
        // DNA: AUG (Met) - GAA (Glu) - GAC (Asp) - UGG (Trp) - UAA (Stop)
        let dna = "ATGGAAGACTGGTAA";
        let protein = translate_dna_to_protein(dna).unwrap();
        assert_eq!(protein, "MEDW*");
    }

    #[test]
    fn test_all_64_codons_resolve() {
        let bases = ['U', 'C', 'A', 'G'];
        let mut count = 0;
        for &b1 in &bases {
            for &b2 in &bases {
                for &b3 in &bases {
                    let codon = format!("{}{}{}", b1, b2, b3);
                    let aa = translate_rna_codon(&codon);
                    assert!(aa.is_some(), "codon {} failed to translate", codon);
                    count += 1;
                }
            }
        }
        assert_eq!(count, 64);
    }
}
