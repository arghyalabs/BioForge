//! # BioForge Parser
//!
//! Hand-written recursive descent parser for the BioLang language.
//!
//! Converts a stream of tokens from [`bioforge_lexer`] into an AST
//! defined in [`bioforge_ast`].

mod parser;

#[cfg(test)]
mod tests;

use bioforge_ast::Program;
use bioforge_diagnostics::Diagnostic;
use bioforge_lexer::Lexer;

/// Parse BioLang source code into a [`Program`] AST.
///
/// Returns the (possibly partial) AST along with any diagnostics.
/// On syntax errors, the parser attempts error recovery to report
/// as many issues as possible in a single pass.
pub fn parse(source: &str) -> (Program, Vec<Diagnostic>) {
    let lexer = Lexer::new(source);
    let (tokens, mut diagnostics) = lexer.tokenize();

    let mut parser = parser::Parser::new(tokens);
    let program = parser.parse_program();
    diagnostics.extend(parser.into_diagnostics());

    (program, diagnostics)
}
