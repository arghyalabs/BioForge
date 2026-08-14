use ariadne::{Color, Label as AriadneLabel, Report, ReportKind, Source};
use std::fmt;

/// Source location represented as byte offsets
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Span {
    pub start: usize,
    pub end: usize,
}

impl Span {
    pub fn new(start: usize, end: usize) -> Self {
        Span { start, end }
    }

    pub fn len(&self) -> usize {
        self.end.saturating_sub(self.start)
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn contains(&self, other: &Span) -> bool {
        self.start <= other.start && self.end >= other.end
    }

    pub fn merge(&self, other: &Span) -> Span {
        Span {
            start: self.start.min(other.start),
            end: self.end.max(other.end),
        }
    }
}

impl fmt::Display for Span {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}..{}", self.start, self.end)
    }
}

/// A value with its source location
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Spanned<T> {
    pub node: T,
    pub span: Span,
}

impl<T> Spanned<T> {
    pub fn new(node: T, span: Span) -> Self {
        Spanned { node, span }
    }
}

impl<T: fmt::Display> fmt::Display for Spanned<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.node)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiagnosticSeverity {
    Error,
    Warning,
    Info,
    Hint,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Label {
    pub message: String,
    pub span: Span,
    pub color: Option<Color>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Diagnostic {
    pub severity: DiagnosticSeverity,
    pub message: String,
    pub labels: Vec<Label>,
    pub help: Option<String>,
    pub note: Option<String>,
}

impl Diagnostic {
    pub fn error(message: impl Into<String>) -> Self {
        Diagnostic {
            severity: DiagnosticSeverity::Error,
            message: message.into(),
            labels: Vec::new(),
            help: None,
            note: None,
        }
    }

    pub fn warning(message: impl Into<String>) -> Self {
        Diagnostic {
            severity: DiagnosticSeverity::Warning,
            message: message.into(),
            labels: Vec::new(),
            help: None,
            note: None,
        }
    }

    pub fn with_label(mut self, span: Span, message: impl Into<String>) -> Self {
        self.labels.push(Label {
            message: message.into(),
            span,
            color: None,
        });
        self
    }

    pub fn with_color_label(mut self, span: Span, message: impl Into<String>, color: Color) -> Self {
        self.labels.push(Label {
            message: message.into(),
            span,
            color: Some(color),
        });
        self
    }

    pub fn with_help(mut self, help: impl Into<String>) -> Self {
        self.help = Some(help.into());
        self
    }

    pub fn with_note(mut self, note: impl Into<String>) -> Self {
        self.note = Some(note.into());
        self
    }
}

pub fn render_diagnostics(filename: &str, source: &str, diagnostics: &[Diagnostic]) {
    for diag in diagnostics {
        let kind = match diag.severity {
            DiagnosticSeverity::Error => ReportKind::Error,
            DiagnosticSeverity::Warning => ReportKind::Warning,
            DiagnosticSeverity::Info => ReportKind::Custom("Info", Color::Blue),
            DiagnosticSeverity::Hint => ReportKind::Custom("Hint", Color::Cyan),
        };

        let mut builder = Report::build(kind, filename, 0)
            .with_message(&diag.message);

        for label in &diag.labels {
            let mut a_label = AriadneLabel::new((filename, label.span.start..label.span.end))
                .with_message(&label.message);
            
            if let Some(color) = label.color {
                a_label = a_label.with_color(color);
            } else {
                a_label = a_label.with_color(match diag.severity {
                    DiagnosticSeverity::Error => Color::Red,
                    DiagnosticSeverity::Warning => Color::Yellow,
                    _ => Color::Blue,
                });
            }
            
            builder = builder.with_label(a_label);
        }

        if let Some(help) = &diag.help {
            builder = builder.with_help(help);
        }
        
        if let Some(note) = &diag.note {
            builder = builder.with_note(note);
        }

        builder.finish().print((filename, Source::from(source))).unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_span_len() {
        let s = Span::new(10, 20);
        assert_eq!(s.len(), 10);
    }

    #[test]
    fn test_span_merge() {
        let s1 = Span::new(10, 20);
        let s2 = Span::new(15, 25);
        assert_eq!(s1.merge(&s2), Span::new(10, 25));
    }

    #[test]
    fn test_diagnostic_builder() {
        let diag = Diagnostic::error("test error")
            .with_label(Span::new(0, 5), "here")
            .with_help("try this");
            
        assert_eq!(diag.severity, DiagnosticSeverity::Error);
        assert_eq!(diag.message, "test error");
        assert_eq!(diag.labels.len(), 1);
        assert_eq!(diag.help, Some("try this".to_string()));
    }
}
