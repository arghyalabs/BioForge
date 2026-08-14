//! # BioForge Lexer
//!
//! Tokenizer for the BioForge language, built on top of [`logos`].
//!
//! ## Design Decision: Unit Tokens
//!
//! Units like `K`, `fs`, `ps`, `nm` are lexed as [`Token::Identifier`].
//! The **parser** determines if an identifier following a number is a unit.
//! This avoids ambiguity between `K` as Kelvin and `K` as a variable name.

use bioforge_diagnostics::{Diagnostic, Span};
use logos::Logos;

/// All tokens in BioForge v0.1.
#[derive(Logos, Debug, Clone, PartialEq)]
#[logos(skip r"[ \t\r\n\f]+")]
#[logos(skip r"//[^\n]*")]
#[logos(skip r"/\*([^*]|\*[^/])*\*/")]
pub enum Token {
    // ─── Keywords ──────────────────────────────────────────────────────
    #[token("experiment")]
    Experiment,
    #[token("protein")]
    Protein,
    #[token("ligand")]
    Ligand,
    #[token("ion")]
    Ion,
    #[token("molecule")]
    Molecule,
    #[token("atom")]
    Atom,
    #[token("environment")]
    Environment,
    #[token("simulate")]
    Simulate,
    #[token("measure")]
    Measure,
    #[token("visualize")]
    Visualize,
    #[token("true")]
    True,
    #[token("false")]
    False,

    // ─── Delimiters & Symbols ──────────────────────────────────────────
    #[token("{")]
    LBrace,
    #[token("}")]
    RBrace,
    #[token("(")]
    LParen,
    #[token(")")]
    RParen,
    #[token("[")]
    LBracket,
    #[token("]")]
    RBracket,
    #[token("=")]
    Equals,
    #[token(",")]
    Comma,
    #[token(".")]
    Dot,
    #[token("+")]
    Plus,
    #[token("-")]
    Minus,
    #[token("*")]
    Star,
    #[token("/")]
    Slash,

    // ─── Literals ──────────────────────────────────────────────────────
    /// Floating-point number: `7.4`, `1.5e-3`
    #[regex(r"[0-9]+\.[0-9]+([eE][+-]?[0-9]+)?", |lex| lex.slice().parse::<f64>().ok())]
    Float(f64),

    /// Integer number: `310`, `42`
    #[regex(r"[0-9]+", |lex| lex.slice().parse::<i64>().ok())]
    Integer(i64),

    /// String literal: `"protein.pdb"`
    #[regex(r#""[^"]*""#, |lex| {
        let s = lex.slice();
        Some(s[1..s.len()-1].to_string())
    })]
    StringLiteral(String),

    // ─── Identifiers ───────────────────────────────────────────────────
    /// Identifier or unit name: `receptor`, `K`, `fs`, `nm`
    #[regex(r"[a-zA-Z_][a-zA-Z0-9_]*", |lex| lex.slice().to_string())]
    Identifier(String),
}

impl Token {
    /// Returns a human-readable name for this token kind.
    pub fn display_name(&self) -> &str {
        match self {
            Token::Experiment => "keyword 'experiment'",
            Token::Protein => "keyword 'protein'",
            Token::Ligand => "keyword 'ligand'",
            Token::Ion => "keyword 'ion'",
            Token::Molecule => "keyword 'molecule'",
            Token::Atom => "keyword 'atom'",
            Token::Environment => "keyword 'environment'",
            Token::Simulate => "keyword 'simulate'",
            Token::Measure => "keyword 'measure'",
            Token::Visualize => "keyword 'visualize'",
            Token::True => "keyword 'true'",
            Token::False => "keyword 'false'",
            Token::LBrace => "'{'",
            Token::RBrace => "'}'",
            Token::LParen => "'('",
            Token::RParen => "')'",
            Token::LBracket => "'['",
            Token::RBracket => "']'",
            Token::Equals => "'='",
            Token::Comma => "','",
            Token::Dot => "'.'",
            Token::Plus => "'+'",
            Token::Minus => "'-'",
            Token::Star => "'*'",
            Token::Slash => "'/'",
            Token::Float(_) => "float literal",
            Token::Integer(_) => "integer literal",
            Token::StringLiteral(_) => "string literal",
            Token::Identifier(_) => "identifier",
        }
    }

    /// Returns true if this token is a keyword.
    pub fn is_keyword(&self) -> bool {
        matches!(
            self,
            Token::Experiment
                | Token::Protein
                | Token::Ligand
                | Token::Ion
                | Token::Molecule
                | Token::Atom
                | Token::Environment
                | Token::Simulate
                | Token::Measure
                | Token::Visualize
                | Token::True
                | Token::False
        )
    }
}

/// Tokenizer for BioForge source code.
pub struct Lexer<'a> {
    source: &'a str,
}

impl<'a> Lexer<'a> {
    /// Create a new lexer for the given source code.
    pub fn new(source: &'a str) -> Self {
        Lexer { source }
    }

    /// Tokenize the entire source, returning tokens with spans and any diagnostics.
    pub fn tokenize(&self) -> (Vec<(Token, Span)>, Vec<Diagnostic>) {
        let mut tokens = Vec::new();
        let mut diagnostics = Vec::new();

        let lex = Token::lexer(self.source);
        for (result, range) in lex.spanned() {
            let span = Span::new(range.start, range.end);
            match result {
                Ok(token) => tokens.push((token, span)),
                Err(()) => {
                    let slice = &self.source[range.start..range.end];
                    diagnostics.push(
                        Diagnostic::error(format!("Unrecognized token: '{}'", slice))
                            .with_label(span, "unexpected character")
                            .with_help("Check for typos or unsupported characters"),
                    );
                }
            }
        }

        (tokens, diagnostics)
    }
}

// ─── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn lex(input: &str) -> Vec<Token> {
        let lexer = Lexer::new(input);
        let (tokens, _) = lexer.tokenize();
        tokens.into_iter().map(|(t, _)| t).collect()
    }

    #[test]
    fn test_keywords() {
        assert_eq!(lex("experiment"), vec![Token::Experiment]);
        assert_eq!(lex("protein"), vec![Token::Protein]);
        assert_eq!(lex("ligand"), vec![Token::Ligand]);
        assert_eq!(lex("ion"), vec![Token::Ion]);
        assert_eq!(lex("molecule"), vec![Token::Molecule]);
        assert_eq!(lex("atom"), vec![Token::Atom]);
        assert_eq!(lex("environment"), vec![Token::Environment]);
        assert_eq!(lex("simulate"), vec![Token::Simulate]);
        assert_eq!(lex("measure"), vec![Token::Measure]);
        assert_eq!(lex("visualize"), vec![Token::Visualize]);
        assert_eq!(lex("true"), vec![Token::True]);
        assert_eq!(lex("false"), vec![Token::False]);
    }

    #[test]
    fn test_numbers() {
        assert_eq!(lex("310"), vec![Token::Integer(310)]);
        assert_eq!(lex("7.4"), vec![Token::Float(7.4)]);
        assert_eq!(lex("1.5"), vec![Token::Float(1.5)]);
    }

    #[test]
    fn test_string_literal() {
        assert_eq!(
            lex(r#""protein.pdb""#),
            vec![Token::StringLiteral("protein.pdb".to_string())]
        );
    }

    #[test]
    fn test_identifiers() {
        assert_eq!(
            lex("receptor"),
            vec![Token::Identifier("receptor".to_string())]
        );
        assert_eq!(lex("K"), vec![Token::Identifier("K".to_string())]);
        assert_eq!(lex("fs"), vec![Token::Identifier("fs".to_string())]);
        assert_eq!(lex("nm"), vec![Token::Identifier("nm".to_string())]);
    }

    #[test]
    fn test_symbols() {
        let tokens = lex("{ } ( ) = , . + - * /");
        assert_eq!(
            tokens,
            vec![
                Token::LBrace,
                Token::RBrace,
                Token::LParen,
                Token::RParen,
                Token::Equals,
                Token::Comma,
                Token::Dot,
                Token::Plus,
                Token::Minus,
                Token::Star,
                Token::Slash,
            ]
        );
    }

    #[test]
    fn test_line_comments() {
        assert_eq!(
            lex("experiment // this is a comment\nprotein"),
            vec![Token::Experiment, Token::Protein]
        );
    }

    #[test]
    fn test_block_comments() {
        assert_eq!(
            lex("experiment /* block comment */ protein"),
            vec![Token::Experiment, Token::Protein]
        );
    }

    #[test]
    fn test_quantity_tokens() {
        // The lexer produces separate tokens; the parser combines them
        assert_eq!(
            lex("310 K"),
            vec![Token::Integer(310), Token::Identifier("K".to_string())]
        );
        assert_eq!(
            lex("1 fs"),
            vec![Token::Integer(1), Token::Identifier("fs".to_string())]
        );
        assert_eq!(
            lex("10.5 nm"),
            vec![Token::Float(10.5), Token::Identifier("nm".to_string())]
        );
    }

    #[test]
    fn test_full_experiment() {
        let source = r#"experiment Test {
    protein p = load_structure("test.pdb")
    environment {
        temperature = 310 K
    }
}"#;
        let lexer = Lexer::new(source);
        let (tokens, diagnostics) = lexer.tokenize();
        assert!(diagnostics.is_empty(), "Expected no diagnostics");
        assert!(!tokens.is_empty());
        // First token should be 'experiment'
        assert_eq!(tokens[0].0, Token::Experiment);
    }

    #[test]
    fn test_unrecognized_token() {
        let lexer = Lexer::new("experiment @");
        let (tokens, diagnostics) = lexer.tokenize();
        assert_eq!(tokens.len(), 1); // 'experiment' only
        assert_eq!(diagnostics.len(), 1);
        assert!(diagnostics[0].message.contains("Unrecognized"));
    }

    #[test]
    fn test_snapshot_hello_bio() {
        let source = r#"experiment HelloBiology {
    protein receptor = load_structure("protein.pdb")
    ligand drug = load("drug.sdf")
    environment {
        temperature = 310 K
        pH = 7.4
    }
    simulate {
        timestep = 1 fs
        duration = 10 ps
    }
}"#;
        let lexer = Lexer::new(source);
        let (tokens, _) = lexer.tokenize();
        let token_names: Vec<&str> = tokens.iter().map(|(t, _)| t.display_name()).collect();
        insta::assert_debug_snapshot!(token_names);
    }
}
