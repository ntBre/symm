use std::{fmt::Display, str::FromStr};

use serde::Deserialize;

#[derive(Debug, Deserialize, PartialEq, Eq, PartialOrd, Ord, Copy, Clone)]
pub enum Irrep {
    // C1
    A,
    // C2
    B,
    // Cs - p = prime
    Ap,
    App,
    // C2v
    A1,
    B2,
    B1,
    A2,
    // D2h
    Ag,
    B1g,
    B2g,
    B3g,
    Au,
    B1u,
    B2u,
    B3u,
}

impl Display for Irrep {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.pad(&format!(
            "{}",
            match self {
                Irrep::A => "A",
                Irrep::B => "B",
                Irrep::Ap => "A'",
                Irrep::App => "A''",
                Irrep::A1 => "A1",
                Irrep::B1 => "B1",
                Irrep::B2 => "B2",
                Irrep::A2 => "A2",
                Irrep::Ag => "Ag",
                Irrep::B1g => "B1g",
                Irrep::B2g => "B2g",
                Irrep::B3g => "B3g",
                Irrep::Au => "Au",
                Irrep::B1u => "B1u",
                Irrep::B2u => "B2u",
                Irrep::B3u => "B3u",
            }
        ))
    }
}

impl FromStr for Irrep {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "A" | "a" => Ok(Irrep::A),
            "B" | "b" => Ok(Irrep::B),
            "A'" | "a'" => Ok(Irrep::Ap),
            "A''" | "a''" => Ok(Irrep::App),
            "A1" | "a1" => Ok(Irrep::A1),
            "B1" | "b1" => Ok(Irrep::B1),
            "B2" | "b2" => Ok(Irrep::B2),
            "A2" | "a2" => Ok(Irrep::A2),
            "Ag" | "ag" => Ok(Irrep::Ag),
            "B1g" | "b1g" => Ok(Irrep::B1g),
            "B2g" | "b2g" => Ok(Irrep::B2g),
            "B3g" | "b3g" => Ok(Irrep::B3g),
            "Au" | "au" => Ok(Irrep::Au),
            "B1u" | "b1u" => Ok(Irrep::B1u),
            "B2u" | "b2u" => Ok(Irrep::B2u),
            "B3u" | "b3u" => Ok(Irrep::B3u),
            _ => Err(()),
        }
    }
}
