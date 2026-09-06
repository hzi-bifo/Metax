use std::fmt;

use clap::ValueEnum;

/// Supported input sequence types.
///
/// Values are parsed once by the CLI and carried as an enum throughout the
/// mapper and profiler, so misspellings cannot silently select the wrong
/// parameter preset.
#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub enum Sequencer {
    #[value(name = "Illumina")]
    Illumina,
    #[value(name = "Nanopore", alias = "ONT")]
    Nanopore,
    #[value(name = "PacBio")]
    PacBio,
    #[value(name = "assembly")]
    Assembly,
}

impl Sequencer {
    #[inline]
    pub const fn is_illumina(self) -> bool {
        matches!(self, Self::Illumina)
    }

    #[inline]
    pub const fn is_long_read(self) -> bool {
        matches!(self, Self::Nanopore | Self::PacBio)
    }

    #[inline]
    pub const fn is_read_platform(self) -> bool {
        matches!(self, Self::Illumina | Self::Nanopore | Self::PacBio)
    }

    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Illumina => "Illumina",
            Self::Nanopore => "Nanopore",
            Self::PacBio => "PacBio",
            Self::Assembly => "assembly",
        }
    }
}

impl fmt::Display for Sequencer {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

#[cfg(test)]
mod tests {
    use clap::ValueEnum;

    use super::Sequencer;

    #[test]
    fn supported_values_parse_case_insensitively_when_requested() {
        assert_eq!(
            Sequencer::from_str("illumina", true),
            Ok(Sequencer::Illumina)
        );
        assert_eq!(
            Sequencer::from_str("NANOPORE", true),
            Ok(Sequencer::Nanopore)
        );
        assert_eq!(Sequencer::from_str("ont", true), Ok(Sequencer::Nanopore));
        assert_eq!(Sequencer::from_str("pacbio", true), Ok(Sequencer::PacBio));
        assert_eq!(
            Sequencer::from_str("ASSEMBLY", true),
            Ok(Sequencer::Assembly)
        );
    }

    #[test]
    fn unknown_values_are_rejected() {
        assert!(Sequencer::from_str("illuminaa", true).is_err());
        assert!(Sequencer::from_str("unknown", true).is_err());
    }
}
