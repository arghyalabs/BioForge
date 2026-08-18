//! Recursive descent parser implementation for BioForge v0.1.

use bioforge_ast::*;
use bioforge_diagnostics::{Diagnostic, Span, Spanned};
use bioforge_lexer::Token;
use bioforge_types::UnitRegistry;

pub(crate) struct Parser {
    tokens: Vec<(Token, Span)>,
    pos: usize,
    diagnostics: Vec<Diagnostic>,
    /// Unit registry — single source of truth for unit recognition.
    /// Used by `try_parse_unit` to determine whether an identifier
    /// following a number is a unit (forming a Quantity) or a regular identifier.
    unit_registry: UnitRegistry,
}

impl Parser {
    pub fn new(tokens: Vec<(Token, Span)>) -> Self {
        Parser {
            tokens,
            pos: 0,
            diagnostics: Vec::new(),
            unit_registry: UnitRegistry::new(),
        }
    }

    pub fn into_diagnostics(self) -> Vec<Diagnostic> {
        self.diagnostics
    }

    // ─── Top-level parsing ─────────────────────────────────────────────

    pub fn parse_program(&mut self) -> Program {
        let mut experiments = Vec::new();
        while !self.at_end() {
            if self.check(&Token::Experiment) {
                if let Some(exp) = self.parse_experiment() {
                    experiments.push(exp);
                }
            } else {
                let span = self.current_span();
                self.diagnostics.push(
                    Diagnostic::error("Expected 'experiment' declaration")
                        .with_label(span, "unexpected token here")
                        .with_help("A BioForge program consists of one or more 'experiment' blocks"),
                );
                self.advance();
            }
        }
        Program { experiments }
    }

    fn parse_experiment(&mut self) -> Option<Spanned<ExperimentDecl>> {
        let start_span = self.expect_token(&Token::Experiment)?;
        let name = self.expect_identifier()?;
        self.expect_token(&Token::LBrace)?;

        let mut body = Vec::new();
        while !self.check(&Token::RBrace) && !self.at_end() {
            if let Some(stmt) = self.parse_statement() {
                body.push(stmt);
            } else {
                // Error recovery: skip to next statement boundary
                self.synchronize();
            }
        }

        let end_span = self.expect_token(&Token::RBrace).unwrap_or(self.current_span());
        let span = start_span.merge(&end_span);

        Some(Spanned::new(ExperimentDecl { name, body }, span))
    }

    // ─── Statement parsing ─────────────────────────────────────────────

    fn parse_statement(&mut self) -> Option<Spanned<Statement>> {
        match self.peek()? {
            Token::Protein => self.parse_entity_decl(EntityKind::Protein),
            Token::Ligand => self.parse_entity_decl(EntityKind::Ligand),
            Token::Ion => self.parse_entity_decl(EntityKind::Ion),
            Token::Molecule => self.parse_entity_decl(EntityKind::Molecule),
            Token::Atom => self.parse_entity_decl(EntityKind::Atom),
            Token::Environment => self.parse_environment_block(),
            Token::Simulate => self.parse_simulate_block(),
            Token::Measure => self.parse_measure_block(),
            Token::Visualize => self.parse_visualize_block(),
            Token::Identifier(_) => self.parse_assignment(),
            _ => {
                let span = self.current_span();
                self.diagnostics.push(
                    Diagnostic::error("Unexpected token in experiment body")
                        .with_label(span, "expected a statement")
                        .with_help("Valid statements: entity declarations (protein, ligand, atom, ...), environment, simulate, measure, visualize"),
                );
                None
            }
        }
    }

    fn parse_entity_decl(&mut self, kind: EntityKind) -> Option<Spanned<Statement>> {
        let (_, kind_span) = self.advance()?;
        let kind = Spanned::new(kind, kind_span);
        let name = self.expect_identifier()?;
        self.expect_token(&Token::Equals)?;
        let initializer = self.parse_expr()?;
        let span = kind.span.merge(&initializer.span);

        Some(Spanned::new(
            Statement::EntityDecl {
                kind,
                name,
                initializer,
            },
            span,
        ))
    }

    fn parse_environment_block(&mut self) -> Option<Spanned<Statement>> {
        let start_span = self.expect_token(&Token::Environment)?;

        // Optional name
        let name = if !self.check(&Token::LBrace) {
            if let Some(Token::Identifier(_)) = self.peek() {
                Some(self.expect_identifier()?)
            } else {
                None
            }
        } else {
            None
        };

        self.expect_token(&Token::LBrace)?;
        let properties = self.parse_properties();
        let end_span = self.expect_token(&Token::RBrace).unwrap_or(self.current_span());
        let span = start_span.merge(&end_span);

        Some(Spanned::new(
            Statement::EnvironmentBlock { name, properties },
            span,
        ))
    }

    fn parse_simulate_block(&mut self) -> Option<Spanned<Statement>> {
        let start_span = self.expect_token(&Token::Simulate)?;
        self.expect_token(&Token::LBrace)?;
        let properties = self.parse_properties();
        let end_span = self.expect_token(&Token::RBrace).unwrap_or(self.current_span());
        let span = start_span.merge(&end_span);

        Some(Spanned::new(Statement::SimulateBlock { properties }, span))
    }

    fn parse_measure_block(&mut self) -> Option<Spanned<Statement>> {
        let start_span = self.expect_token(&Token::Measure)?;
        self.expect_token(&Token::LBrace)?;

        let mut measurements = Vec::new();
        while !self.check(&Token::RBrace) && !self.at_end() {
            if let Some(expr) = self.parse_expr() {
                measurements.push(expr);
            } else {
                self.synchronize();
            }
        }

        let end_span = self.expect_token(&Token::RBrace).unwrap_or(self.current_span());
        let span = start_span.merge(&end_span);

        Some(Spanned::new(Statement::MeasureBlock { measurements }, span))
    }

    fn parse_visualize_block(&mut self) -> Option<Spanned<Statement>> {
        let start_span = self.expect_token(&Token::Visualize)?;
        self.expect_token(&Token::LBrace)?;

        let mut targets = Vec::new();
        while !self.check(&Token::RBrace) && !self.at_end() {
            if let Some(expr) = self.parse_expr() {
                targets.push(expr);
            } else {
                self.synchronize();
            }
        }

        let end_span = self.expect_token(&Token::RBrace).unwrap_or(self.current_span());
        let span = start_span.merge(&end_span);

        Some(Spanned::new(Statement::VisualizeBlock { targets }, span))
    }

    fn parse_assignment(&mut self) -> Option<Spanned<Statement>> {
        let name = self.expect_identifier()?;

        if !self.check(&Token::Equals) {
            // Not an assignment — could be just an expression statement.
            // For now, treat unexpected identifiers as errors.
            let span = name.span;
            self.diagnostics.push(
                Diagnostic::error(format!("Expected '=' after '{}'", name.node))
                    .with_label(span, "expected '=' here"),
            );
            return None;
        }

        self.expect_token(&Token::Equals)?;
        let value = self.parse_expr()?;
        let span = name.span.merge(&value.span);

        Some(Spanned::new(Statement::Assignment { name, value }, span))
    }

    // ─── Property parsing ──────────────────────────────────────────────

    fn parse_properties(&mut self) -> Vec<Spanned<Property>> {
        let mut props = Vec::new();
        while !self.check(&Token::RBrace) && !self.at_end() {
            if let Some(prop) = self.parse_property() {
                props.push(prop);
            } else {
                self.synchronize();
            }
        }
        props
    }

    fn parse_property(&mut self) -> Option<Spanned<Property>> {
        let key = self.expect_identifier()?;
        self.expect_token(&Token::Equals)?;
        let value = self.parse_expr()?;
        let span = key.span.merge(&value.span);

        Some(Spanned::new(Property { key, value }, span))
    }

    // ─── Expression parsing ────────────────────────────────────────────

    fn parse_expr(&mut self) -> Option<Spanned<Expr>> {
        self.parse_additive()
    }

    /// Parse additive expressions: `a + b`, `a - b`
    fn parse_additive(&mut self) -> Option<Spanned<Expr>> {
        let mut left = self.parse_multiplicative()?;

        while let Some(token) = self.peek() {
            let op = match token {
                Token::Plus => BinOp::Add,
                Token::Minus => BinOp::Sub,
                _ => break,
            };
            let (_, op_span) = self.advance().unwrap();
            let right = self.parse_multiplicative()?;
            let span = left.span.merge(&right.span);
            left = Spanned::new(
                Expr::BinaryOp {
                    left: Box::new(left),
                    op: Spanned::new(op, op_span),
                    right: Box::new(right),
                },
                span,
            );
        }

        Some(left)
    }

    /// Parse multiplicative expressions: `a * b`, `a / b`
    fn parse_multiplicative(&mut self) -> Option<Spanned<Expr>> {
        let mut left = self.parse_primary()?;

        while let Some(token) = self.peek() {
            let op = match token {
                Token::Star => BinOp::Mul,
                Token::Slash => BinOp::Div,
                _ => break,
            };
            let (_, op_span) = self.advance().unwrap();
            let right = self.parse_primary()?;
            let span = left.span.merge(&right.span);
            left = Spanned::new(
                Expr::BinaryOp {
                    left: Box::new(left),
                    op: Spanned::new(op, op_span),
                    right: Box::new(right),
                },
                span,
            );
        }

        Some(left)
    }

    /// Parse primary expressions: literals, identifiers, function calls, quantities.
    fn parse_primary(&mut self) -> Option<Spanned<Expr>> {
        let (token, span) = self.advance()?;

        match token {
            // ── Numbers: check if followed by a unit ──
            Token::Integer(n) => {
                let value = n as f64;
                if let Some(unit) = self.try_parse_unit() {
                    let full_span = span.merge(&unit.span);
                    Some(Spanned::new(Expr::Quantity { value, unit }, full_span))
                } else {
                    Some(Spanned::new(Expr::NumberLiteral(value), span))
                }
            }
            Token::Float(n) => {
                if let Some(unit) = self.try_parse_unit() {
                    let full_span = span.merge(&unit.span);
                    Some(Spanned::new(Expr::Quantity { value: n, unit }, full_span))
                } else {
                    Some(Spanned::new(Expr::NumberLiteral(n), span))
                }
            }

            // ── String literal ──
            Token::StringLiteral(s) => Some(Spanned::new(Expr::StringLiteral(s), span)),

            // ── Boolean literals ──
            Token::True => Some(Spanned::new(Expr::BoolLiteral(true), span)),
            Token::False => Some(Spanned::new(Expr::BoolLiteral(false), span)),

            // ── Identifiers and function calls ──
            Token::Identifier(name) => {
                if self.check(&Token::LParen) {
                    // Function call
                    self.parse_function_call(name, span)
                } else {
                    Some(Spanned::new(Expr::Identifier(name), span))
                }
            }

            // ── Keyword-as-identifier: allow entity keywords to be used as identifiers
            // in expression position (e.g., `atom("H")`)
            Token::Atom => {
                let name = "atom".to_string();
                if self.check(&Token::LParen) {
                    self.parse_function_call(name, span)
                } else {
                    Some(Spanned::new(Expr::Identifier(name), span))
                }
            }
            Token::Ion => {
                let name = "ion".to_string();
                if self.check(&Token::LParen) {
                    self.parse_function_call(name, span)
                } else {
                    Some(Spanned::new(Expr::Identifier(name), span))
                }
            }

            // ── Parenthesized expression ──
            Token::LParen => {
                let expr = self.parse_expr()?;
                self.expect_token(&Token::RParen)?;
                Some(expr)
            }

            _ => {
                self.diagnostics.push(
                    Diagnostic::error("Expected expression")
                        .with_label(span, "unexpected token")
                        .with_help(
                            "Expected a number, string, identifier, or function call",
                        ),
                );
                None
            }
        }
    }

    /// Try to consume a unit identifier following a number.
    fn try_parse_unit(&mut self) -> Option<Spanned<String>> {
        if let Some(Token::Identifier(name)) = self.peek() {
            if self.unit_registry.is_known(name) {
                let name = name.clone();
                let (_, span) = self.advance().unwrap();
                return Some(Spanned::new(name, span));
            }
        }
        None
    }

    /// Parse function call arguments: `name(arg1, arg2, ...)`
    fn parse_function_call(&mut self, name: String, name_span: Span) -> Option<Spanned<Expr>> {
        self.expect_token(&Token::LParen)?;

        let mut args = Vec::new();
        if !self.check(&Token::RParen) {
            if let Some(arg) = self.parse_expr() {
                args.push(arg);
            }
            while self.check(&Token::Comma) {
                self.advance(); // consume comma
                if let Some(arg) = self.parse_expr() {
                    args.push(arg);
                }
            }
        }

        let end_span = self.expect_token(&Token::RParen).unwrap_or(self.current_span());
        let span = name_span.merge(&end_span);

        Some(Spanned::new(
            Expr::FunctionCall {
                name: Spanned::new(name, name_span),
                args,
            },
            span,
        ))
    }

    // ─── Helpers ───────────────────────────────────────────────────────

    fn peek(&self) -> Option<&Token> {
        self.tokens.get(self.pos).map(|(t, _)| t)
    }

    fn check(&self, expected: &Token) -> bool {
        match self.peek() {
            Some(t) => std::mem::discriminant(t) == std::mem::discriminant(expected),
            None => false,
        }
    }

    fn advance(&mut self) -> Option<(Token, Span)> {
        if self.pos < self.tokens.len() {
            let result = self.tokens[self.pos].clone();
            self.pos += 1;
            Some(result)
        } else {
            None
        }
    }

    fn expect_token(&mut self, expected: &Token) -> Option<Span> {
        if self.check(expected) {
            let (_, span) = self.advance().unwrap();
            Some(span)
        } else {
            let span = self.current_span();
            let found = self
                .peek()
                .map(|t| t.display_name().to_string())
                .unwrap_or_else(|| "end of file".to_string());
            self.diagnostics.push(
                Diagnostic::error(format!(
                    "Expected {}, found {}",
                    expected.display_name(),
                    found
                ))
                .with_label(span, format!("expected {}", expected.display_name())),
            );
            None
        }
    }

    fn expect_identifier(&mut self) -> Option<Spanned<String>> {
        match self.peek() {
            Some(Token::Identifier(_)) => {
                if let Some((Token::Identifier(name), span)) = self.advance() {
                    Some(Spanned::new(name, span))
                } else {
                    None
                }
            }
            _ => {
                let span = self.current_span();
                let found = self
                    .peek()
                    .map(|t| t.display_name().to_string())
                    .unwrap_or_else(|| "end of file".to_string());
                self.diagnostics.push(
                    Diagnostic::error(format!("Expected identifier, found {}", found))
                        .with_label(span, "expected identifier"),
                );
                None
            }
        }
    }

    fn current_span(&self) -> Span {
        if self.pos < self.tokens.len() {
            self.tokens[self.pos].1
        } else if !self.tokens.is_empty() {
            let last = &self.tokens[self.tokens.len() - 1];
            Span::new(last.1.end, last.1.end)
        } else {
            Span::new(0, 0)
        }
    }

    fn at_end(&self) -> bool {
        self.pos >= self.tokens.len()
    }

    /// Error recovery: skip tokens until we find a `}` or the start of a new statement.
    fn synchronize(&mut self) {
        while !self.at_end() {
            match self.peek() {
                Some(Token::RBrace) => return,
                Some(Token::Protein)
                | Some(Token::Ligand)
                | Some(Token::Ion)
                | Some(Token::Molecule)
                | Some(Token::Atom)
                | Some(Token::Environment)
                | Some(Token::Simulate)
                | Some(Token::Measure)
                | Some(Token::Visualize)
                | Some(Token::Experiment) => return,
                _ => {
                    self.advance();
                }
            }
        }
    }
}
