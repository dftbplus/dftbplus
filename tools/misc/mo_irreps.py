#!/usr/bin/env python3
"""Assign molecular-orbital irreducible representations from DFTB+ output."""

from __future__ import annotations

import argparse
import dataclasses
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"
REAL_LINE = re.compile(rf"^(?P<label>.*?)\s+(?P<coef>{FLOAT})\s+(?P<frac>{FLOAT})\s*$")
CMPLX_LINE = re.compile(
    rf"^(?P<label>.*?)\s+\(\s*(?P<real>{FLOAT})\s*,\s*(?P<imag>{FLOAT})\s*\)\s+"
    rf"(?P<frac>{FLOAT})\s*$"
)
REAL_HEADER = re.compile(r"Eigenvector:\s*(?P<state>\d+)\s*\((?P<spin>[^)]*)\)")
K_HEADER = re.compile(
    r"K-point:\s*(?P<kpt>\d+)\s+Eigenvector:\s*(?P<state>\d+)"
    r"(?:\s*\((?P<spin>[^)]*)\))?"
)
BAND_HEADER = re.compile(
    r"KPT\s+(?P<kpt>\d+)(?:\s+SPIN\s+(?P<spin>\d+))?(?:\s+KWEIGHT\s+" + FLOAT + r")?"
)

ORBITAL_LABELS = [
    "y(3x2-y2)",
    "x2+y2+z2",
    "z(x2-y2)",
    "x(x2-3y2)",
    "x2-y2",
    "xy",
    "yz",
    "z2",
    "xz",
    "x",
    "y",
    "z",
    "s",
]

POLYS = {
    "s": {(0, 0, 0): 1.0},
    "x": {(1, 0, 0): 1.0},
    "y": {(0, 1, 0): 1.0},
    "z": {(0, 0, 1): 1.0},
    "xy": {(1, 1, 0): 1.0},
    "yz": {(0, 1, 1): 1.0},
    "z2": {(0, 0, 2): 2.0, (2, 0, 0): -1.0, (0, 2, 0): -1.0},
    "xz": {(1, 0, 1): 1.0},
    "x2-y2": {(2, 0, 0): 1.0, (0, 2, 0): -1.0},
    "y(3x2-y2)": {(2, 1, 0): 3.0, (0, 3, 0): -1.0},
    "x2+y2+z2": {(2, 0, 1): 1.0, (0, 2, 1): 1.0, (0, 0, 3): 1.0},
    "yz2": {(0, 1, 2): 1.0},
    "z3": {(0, 0, 3): 1.0},
    "xz2": {(1, 0, 2): 1.0},
    "z(x2-y2)": {(2, 0, 1): 1.0, (0, 2, 1): -1.0},
    "x(x2-3y2)": {(3, 0, 0): 1.0, (1, 2, 0): -3.0},
}

SHELLS = {
    "s": ["s"],
    "p": ["y", "z", "x"],
    "d": ["xy", "yz", "z2", "xz", "x2-y2"],
    "f": ["y(3x2-y2)", "x2+y2+z2", "yz2", "z3", "xz2", "z(x2-y2)", "x(x2-3y2)"],
}

D6H_CLASSES = [
    ("E", 1),
    ("2C6", 2),
    ("2C3", 2),
    ("C2", 1),
    ("3C2p", 3),
    ("3C2pp", 3),
    ("i", 1),
    ("2S3", 2),
    ("2S6", 2),
    ("sh", 1),
    ("3sd", 3),
    ("3sv", 3),
]

D6H_CHARACTERS = {
    "A1g": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    "A2g": [1, 1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1],
    "B1g": [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1],
    "B2g": [1, -1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1],
    "E1g": [2, 1, -1, -2, 0, 0, 2, 1, -1, -2, 0, 0],
    "E2g": [2, -1, -1, 2, 0, 0, 2, -1, -1, 2, 0, 0],
    "A1u": [1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1],
    "A2u": [1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1],
    "B1u": [1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1],
    "B2u": [1, -1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1],
    "E1u": [2, 1, -1, -2, 0, 0, -2, -1, 1, 2, 0, 0],
    "E2u": [2, -1, -1, 2, 0, 0, -2, 1, 1, -2, 0, 0],
}

POINT_GROUPS = {
    "C1": {
        "classes": [("E", 1)],
        "characters": {"A": [1]},
    },
    "Ci": {
        "classes": [("E", 1), ("i", 1)],
        "characters": {"Ag": [1, 1], "Au": [1, -1]},
    },
    "Cs": {
        "classes": [("E", 1), ("s", 1)],
        "characters": {"A'": [1, 1], 'A"': [1, -1]},
    },
    "C2": {
        "classes": [("E", 1), ("C2", 1)],
        "characters": {"A": [1, 1], "B": [1, -1]},
    },
    "C2v": {
        "classes": [("E", 1), ("C2", 1), ("sv", 1), ("svp", 1)],
        "characters": {
            "A1": [1, 1, 1, 1],
            "A2": [1, 1, -1, -1],
            "B1": [1, -1, 1, -1],
            "B2": [1, -1, -1, 1],
        },
    },
    "C2h": {
        "classes": [("E", 1), ("C2", 1), ("i", 1), ("sh", 1)],
        "characters": {
            "Ag": [1, 1, 1, 1],
            "Bg": [1, -1, 1, -1],
            "Au": [1, 1, -1, -1],
            "Bu": [1, -1, -1, 1],
        },
    },
    "D2": {
        "classes": [("E", 1), ("C2z", 1), ("C2y", 1), ("C2x", 1)],
        "characters": {
            "A": [1, 1, 1, 1],
            "B1": [1, 1, -1, -1],
            "B2": [1, -1, 1, -1],
            "B3": [1, -1, -1, 1],
        },
    },
    "D2h": {
        "classes": [
            ("E", 1),
            ("C2z", 1),
            ("C2y", 1),
            ("C2x", 1),
            ("i", 1),
            ("sxy", 1),
            ("sxz", 1),
            ("syz", 1),
        ],
        "characters": {
            "Ag": [1, 1, 1, 1, 1, 1, 1, 1],
            "B1g": [1, 1, -1, -1, 1, 1, -1, -1],
            "B2g": [1, -1, 1, -1, 1, -1, 1, -1],
            "B3g": [1, -1, -1, 1, 1, -1, -1, 1],
            "Au": [1, 1, 1, 1, -1, -1, -1, -1],
            "B1u": [1, 1, -1, -1, -1, -1, 1, 1],
            "B2u": [1, -1, 1, -1, -1, 1, -1, 1],
            "B3u": [1, -1, -1, 1, -1, 1, 1, -1],
        },
    },
    "C3v": {
        "classes": [("E", 1), ("2C3", 2), ("3sv", 3)],
        "characters": {"A1": [1, 1, 1], "A2": [1, 1, -1], "E": [2, -1, 0]},
    },
    "C4v": {
        "classes": [("E", 1), ("2C4", 2), ("C2", 1), ("2sv", 2), ("2sd", 2)],
        "characters": {
            "A1": [1, 1, 1, 1, 1],
            "A2": [1, 1, 1, -1, -1],
            "B1": [1, -1, 1, 1, -1],
            "B2": [1, -1, 1, -1, 1],
            "E": [2, 0, -2, 0, 0],
        },
    },
    "C6v": {
        "classes": [("E", 1), ("2C6", 2), ("2C3", 2), ("C2", 1), ("3sv", 3), ("3sd", 3)],
        "characters": {
            "A1": [1, 1, 1, 1, 1, 1],
            "A2": [1, 1, 1, 1, -1, -1],
            "B1": [1, -1, 1, -1, 1, -1],
            "B2": [1, -1, 1, -1, -1, 1],
            "E1": [2, 1, -1, -2, 0, 0],
            "E2": [2, -1, -1, 2, 0, 0],
        },
    },
    "C3h": {
        "classes": [("E", 1), ("C3+", 1), ("C3-", 1), ("sh", 1), ("S3+", 1), ("S3-", 1)],
        "characters": {
            "A'": [1, 1, 1, 1, 1, 1],
            "E1'": [1, np.exp(2j * np.pi / 3), np.exp(-2j * np.pi / 3), 1, np.exp(2j * np.pi / 3), np.exp(-2j * np.pi / 3)],
            "E2'": [1, np.exp(-2j * np.pi / 3), np.exp(2j * np.pi / 3), 1, np.exp(-2j * np.pi / 3), np.exp(2j * np.pi / 3)],
            'A"': [1, 1, 1, -1, -1, -1],
            'E1"': [1, np.exp(2j * np.pi / 3), np.exp(-2j * np.pi / 3), -1, -np.exp(2j * np.pi / 3), -np.exp(-2j * np.pi / 3)],
            'E2"': [1, np.exp(-2j * np.pi / 3), np.exp(2j * np.pi / 3), -1, -np.exp(-2j * np.pi / 3), -np.exp(2j * np.pi / 3)],
        },
    },
    "S4": {
        "classes": [("E", 1), ("S4+", 1), ("C2", 1), ("S4-", 1)],
        "characters": {
            "A": [1, 1, 1, 1],
            "E+": [1, 1j, -1, -1j],
            "B": [1, -1, 1, -1],
            "E-": [1, -1j, -1, 1j],
        },
    },
    "D2d": {
        "classes": [("E", 1), ("2S4", 2), ("C2", 1), ("2C2p", 2), ("2sd", 2)],
        "characters": {
            "A1": [1, 1, 1, 1, 1], "A2": [1, 1, 1, -1, -1],
            "B1": [1, -1, 1, 1, -1], "B2": [1, -1, 1, -1, 1],
            "E": [2, 0, -2, 0, 0],
        },
    },
    "D3h": {
        "classes": [
            ("E", 1),
            ("2C3", 2),
            ("3C2p", 3),
            ("sh", 1),
            ("2S3", 2),
            ("3sv", 3),
        ],
        "characters": {
            "A1'": [1, 1, 1, 1, 1, 1],
            "A2'": [1, 1, -1, 1, 1, -1],
            "E'": [2, -1, 0, 2, -1, 0],
            'A1"': [1, 1, 1, -1, -1, -1],
            'A2"': [1, 1, -1, -1, -1, 1],
            'E"': [2, -1, 0, -2, 1, 0],
        },
    },
    "D3d": {
        "classes": [("E", 1), ("2C3", 2), ("3C2p", 3), ("i", 1), ("2S6", 2), ("3sd", 3)],
        "characters": {
            "A1g": [1, 1, 1, 1, 1, 1], "A2g": [1, 1, -1, 1, 1, -1],
            "Eg": [2, -1, 0, 2, -1, 0],
            "A1u": [1, 1, 1, -1, -1, -1], "A2u": [1, 1, -1, -1, -1, 1],
            "Eu": [2, -1, 0, -2, 1, 0],
        },
    },
    "D4h": {
        "classes": [("E", 1), ("2C4", 2), ("C2", 1), ("2C2p", 2), ("2C2pp", 2), ("i", 1), ("2S4", 2), ("sh", 1), ("2sv", 2), ("2sd", 2)],
        "characters": {
            "A1g": [1,1,1,1,1, 1,1,1,1,1], "A2g": [1,1,1,-1,-1, 1,1,1,-1,-1],
            "B1g": [1,-1,1,1,-1, 1,-1,1,1,-1], "B2g": [1,-1,1,-1,1, 1,-1,1,-1,1],
            "Eg": [2,0,-2,0,0, 2,0,-2,0,0],
            "A1u": [1,1,1,1,1, -1,-1,-1,-1,-1], "A2u": [1,1,1,-1,-1, -1,-1,-1,1,1],
            "B1u": [1,-1,1,1,-1, -1,1,-1,-1,1], "B2u": [1,-1,1,-1,1, -1,1,-1,1,-1],
            "Eu": [2,0,-2,0,0, -2,0,2,0,0],
        },
    },
    "Td": {
        "classes": [("E", 1), ("8C3", 8), ("3C2", 3), ("6S4", 6), ("6sd", 6)],
        "characters": {
            "A1": [1,1,1,1,1], "A2": [1,1,1,-1,-1], "E": [2,-1,2,0,0],
            "T1": [3,0,-1,1,-1], "T2": [3,0,-1,-1,1],
        },
    },
    "Oh": {
        "classes": [("E",1),("8C3",8),("6C2p",6),("6C4",6),("3C2",3),("i",1),("8S6",8),("6sd",6),("6S4",6),("3sh",3)],
        "characters": {
            "A1g":[1,1,1,1,1, 1,1,1,1,1], "A2g":[1,1,-1,-1,1, 1,1,-1,-1,1],
            "Eg":[2,-1,0,0,2, 2,-1,0,0,2], "T1g":[3,0,-1,1,-1, 3,0,-1,1,-1],
            "T2g":[3,0,1,-1,-1, 3,0,1,-1,-1],
            "A1u":[1,1,1,1,1, -1,-1,-1,-1,-1], "A2u":[1,1,-1,-1,1, -1,-1,1,1,-1],
            "Eu":[2,-1,0,0,2, -2,1,0,0,-2], "T1u":[3,0,-1,1,-1, -3,0,1,-1,1],
            "T2u":[3,0,1,-1,-1, -3,0,-1,1,1],
        },
    },
    "D6h": {"classes": D6H_CLASSES, "characters": D6H_CHARACTERS},
}

ORIENTATION_WARNINGS = {
    "C2v": (
        "exchanging the x and y conventions swaps B1 <-> B2; "
        "A1 and A2 are unchanged"
    ),
    "D2": (
        "permuting the x, y, and z conventions permutes B1, B2, and B3; "
        "A is unchanged"
    ),
    "D2h": (
        "permuting the x, y, and z conventions permutes B1g/B2g/B3g and "
        "B1u/B2u/B3u separately; Ag, Au, and g/u parity are unchanged"
    ),
    "D2d": (
        "rotating the in-plane convention by 45 degrees swaps B1 <-> B2; "
        "A1, A2, and E are unchanged"
    ),
    "C4v": (
        "rotating the in-plane convention by 45 degrees swaps B1 <-> B2; "
        "A1, A2, and E are unchanged"
    ),
    "C6v": (
        "rotating the in-plane convention by 30 degrees swaps B1 <-> B2; "
        "the A and E labels are unchanged"
    ),
    "D4h": (
        "rotating the in-plane convention by 45 degrees swaps B1g <-> B2g "
        "and B1u <-> B2u; the A and E labels and g/u parity are unchanged"
    ),
    "D6h": (
        "rotating the in-plane convention by 30 degrees swaps B1g <-> B2g "
        "and B1u <-> B2u; the A and E labels and g/u parity are unchanged"
    ),
}

SUBGROUP_SUGGESTIONS = {
    "T": ["D2", "C2v", "C2", "Cs", "C1"],
    "Td": ["D2", "C2v", "C3v", "Cs", "C1"],
    "Th": ["D2h", "D2", "C2h", "C2v", "Ci", "Cs", "C1"],
    "O": ["D2", "C4v", "C2v", "C2", "C1"],
    "Oh": ["D2h", "C4v", "C2v", "C2h", "Ci", "Cs", "C1"],
    "I": ["D2", "C2", "C1"],
    "Ih": ["D2h", "D2", "C2h", "C2v", "Ci", "Cs", "C1"],
    "D3h": ["C3v", "C2v", "Cs", "C1"],
    "D4h": ["C4v", "D2h", "C2v", "C2h", "Cs", "Ci", "C1"],
    "D5h": ["C2v", "Cs", "C1"],
    "D6d": ["C6v", "D2", "C2v", "C1"],
}


@dataclasses.dataclass
class AO:
    atom: int
    species: str
    label: str
    shell: str


@dataclasses.dataclass
class Eigenvector:
    state: int
    spin: str
    kpt: int
    coeffs: np.ndarray
    fractions: np.ndarray


@dataclasses.dataclass
class BandState:
    kpt: int
    spin: int
    state: int
    energy: float
    occ: float


def parse_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    lines = path.read_text().splitlines()
    if not lines:
        raise ValueError(f"{path}: empty geometry file")
    natom = int(lines[0].strip())
    atom_lines = [line for line in lines[2:] if line.strip()][:natom]
    if len(atom_lines) != natom:
        raise ValueError(f"{path}: expected {natom} atoms, found {len(atom_lines)}")
    species = []
    coords = []
    for line in atom_lines:
        words = line.split()
        species.append(words[0])
        coords.append([float(words[1]), float(words[2]), float(words[3])])
    arr = np.array(coords, dtype=float)
    arr -= arr.mean(axis=0)
    return species, arr


def parse_ao_label(raw: str, cur_atom: int | None, cur_species: str | None) -> tuple[int, str, str]:
    text = raw.strip()
    words = text.split()
    if words and words[0].isdigit():
        cur_atom = int(words[0])
        cur_species = words[1]
        body = " ".join(words[2:])
    else:
        body = text
    if cur_atom is None or cur_species is None:
        raise ValueError(f"AO line does not start with an atom: {raw!r}")
    for label in ORBITAL_LABELS:
        if body.endswith(label):
            return cur_atom, cur_species, label
    raise ValueError(f"Could not infer AO label from {raw!r}")


def shell_for_label(label: str) -> str:
    for shell, labels in SHELLS.items():
        if label in labels:
            return shell
    raise ValueError(f"Unsupported AO label {label!r}")


def parse_eigenvec(path: Path) -> tuple[list[AO], list[Eigenvector]]:
    aos: list[AO] = []
    vecs: list[Eigenvector] = []
    cur_meta: list[AO] = []
    cur_coeffs: list[complex] = []
    cur_fracs: list[float] = []
    cur_state: int | None = None
    cur_spin = "1"
    cur_kpt = 1
    cur_atom: int | None = None
    cur_species: str | None = None

    def flush() -> None:
        nonlocal aos, cur_meta, cur_coeffs, cur_fracs, cur_state
        if cur_state is None:
            return
        if not aos:
            aos = cur_meta
        elif [(ao.atom, ao.species, ao.label) for ao in aos] != [
            (ao.atom, ao.species, ao.label) for ao in cur_meta
        ]:
            raise ValueError("AO order changes between eigenvectors")
        vecs.append(
            Eigenvector(
                cur_state,
                cur_spin,
                cur_kpt,
                np.array(cur_coeffs, dtype=complex),
                np.array(cur_fracs, dtype=float),
            )
        )
        cur_meta = []
        cur_coeffs = []
        cur_fracs = []
        cur_state = None

    for line in path.read_text().splitlines():
        match = K_HEADER.search(line)
        if match:
            flush()
            cur_kpt = int(match.group("kpt"))
            cur_state = int(match.group("state"))
            cur_spin = (match.group("spin") or "1").strip()
            cur_atom = None
            cur_species = None
            continue
        match = REAL_HEADER.search(line)
        if match:
            flush()
            cur_kpt = 1
            cur_state = int(match.group("state"))
            cur_spin = match.group("spin").strip()
            cur_atom = None
            cur_species = None
            continue
        if cur_state is None or not line.strip():
            continue
        match = CMPLX_LINE.match(line)
        if match:
            atom, species, label = parse_ao_label(match.group("label"), cur_atom, cur_species)
            cur_atom, cur_species = atom, species
            cur_meta.append(AO(atom, species, label, shell_for_label(label)))
            cur_coeffs.append(complex(float(match.group("real")), float(match.group("imag"))))
            cur_fracs.append(float(match.group("frac")))
            continue
        match = REAL_LINE.match(line)
        if match:
            atom, species, label = parse_ao_label(match.group("label"), cur_atom, cur_species)
            cur_atom, cur_species = atom, species
            cur_meta.append(AO(atom, species, label, shell_for_label(label)))
            cur_coeffs.append(complex(float(match.group("coef")), 0.0))
            cur_fracs.append(float(match.group("frac")))
    flush()
    if not vecs:
        raise ValueError(f"{path}: no eigenvectors found")
    return aos, vecs


def parse_band(path: Path) -> list[BandState]:
    states: list[BandState] = []
    kpt = 1
    spin = 1
    state = 0
    for line in path.read_text().splitlines():
        header = BAND_HEADER.search(line)
        if header:
            kpt = int(header.group("kpt"))
            spin = int(header.group("spin") or 1)
            state = 0
            continue
        words = line.split()
        if len(words) < 2:
            continue
        try:
            vals = [float(word) for word in words]
        except ValueError:
            continue
        if len(vals) >= 3 and abs(vals[0] - (state + 1)) < 1e-8:
            vals = vals[1:]
        state += 1
        states.append(BandState(kpt, spin, state, vals[0], vals[1]))
    return states


def parse_oversqr(path: Path) -> dict[int, np.ndarray]:
    tokens = path.read_text().split()
    if len(tokens) < 8 or tokens[0] != "#":
        raise ValueError(f"{path}: unsupported oversqr.dat header")
    is_real = tokens[4].upper().startswith("T")
    norb = int(tokens[5])
    nkpt = int(tokens[6])
    pos = 7
    matrices = {}
    for _ in range(nkpt):
        while pos < len(tokens) and tokens[pos] == "#":
            pos += 1
            while pos < len(tokens) and not re.fullmatch(r"[+-]?\d+", tokens[pos]):
                pos += 1
        ikpt = int(tokens[pos])
        pos += 2 if pos + 1 < len(tokens) and re.fullmatch(r"[+-]?\d+", tokens[pos + 1]) else 1
        while pos < len(tokens) and tokens[pos] == "#":
            pos += 2
        nvals = norb * norb if is_real else 2 * norb * norb
        vals = np.array([float(tok) for tok in tokens[pos : pos + nvals]])
        pos += nvals
        if is_real:
            mat = vals.reshape((norb, norb), order="F")
        else:
            mat = (vals[0::2] + 1j * vals[1::2]).reshape((norb, norb), order="F")
        matrices[ikpt] = mat
    return matrices


def unit_vector(vec: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(vec)
    if norm < 1e-12:
        raise ValueError("Cannot normalize a zero vector")
    return vec / norm


def rotation_angle(op: np.ndarray) -> float:
    arg = (np.trace(op) - 1.0) / 2.0
    return math.acos(float(np.clip(arg, -1.0, 1.0)))


def c2_axis(op: np.ndarray) -> np.ndarray:
    vals, vecs = np.linalg.eig(op)
    idx = int(np.argmin(np.abs(vals - 1.0)))
    return unit_vector(np.real(vecs[:, idx]))


def angle_about_axis(vec: np.ndarray, ref: np.ndarray, normal: np.ndarray) -> float:
    x = unit_vector(ref - np.dot(ref, normal) * normal)
    y = np.cross(normal, x)
    projected = unit_vector(vec - np.dot(vec, normal) * normal)
    return math.atan2(float(np.dot(projected, y)), float(np.dot(projected, x)))


def infer_d6_axis(coords: np.ndarray) -> np.ndarray:
    centered = coords - coords.mean(axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    return unit_vector(vh[-1])


def infer_c2_reference(coords: np.ndarray, axis: np.ndarray) -> np.ndarray:
    centered = coords - coords.mean(axis=0)
    projected = centered - np.outer(centered @ axis, axis)
    norms = np.linalg.norm(projected, axis=1)
    if np.max(norms) < 1e-12:
        raise ValueError("Could not infer an in-plane D6h reference axis from geometry")
    return unit_vector(projected[int(np.argmax(norms))])


def classify_proper_d6h_operation(op: np.ndarray, axis: np.ndarray, ref: np.ndarray) -> str:
    angle = rotation_angle(op)
    if angle < 1e-5:
        return "E"
    if abs(angle - math.pi / 3.0) < 1e-4:
        return "2C6"
    if abs(angle - 2.0 * math.pi / 3.0) < 1e-4:
        return "2C3"
    if abs(angle - math.pi) > 1e-4:
        raise ValueError(f"Cannot classify D6h proper rotation angle {angle:.8f}")
    rot_axis = c2_axis(op)
    if abs(float(np.dot(rot_axis, axis))) > 0.9:
        return "C2"
    phi = angle_about_axis(rot_axis, ref, axis)
    sector = (phi / (math.pi / 6.0)) % 2.0
    return "3C2p" if min(abs(sector), abs(sector - 2.0)) < 0.5 else "3C2pp"


def classify_d6h_operation(op: np.ndarray, axis: np.ndarray, ref: np.ndarray) -> str:
    det = np.linalg.det(op)
    if det > 0.0:
        return classify_proper_d6h_operation(op, axis, ref)
    proper_name = classify_proper_d6h_operation(-op, axis, ref)
    return {
        "E": "i",
        "2C6": "2S3",
        "2C3": "2S6",
        "C2": "sh",
        "3C2p": "3sd",
        "3C2pp": "3sv",
    }[proper_name]


def classify_d3h_operation(op: np.ndarray, axis: np.ndarray) -> str:
    """Classify a D3h operation relative to its principal C3 axis."""
    if np.linalg.det(op) > 0.0:
        angle = rotation_angle(op)
        if angle < 1e-5:
            return "E"
        if abs(angle - 2.0 * math.pi / 3.0) < 1e-4:
            return "2C3"
        if abs(angle - math.pi) < 1e-4:
            if abs(float(np.dot(operation_axis(op), axis))) < 0.1:
                return "3C2p"
        raise ValueError(f"Cannot classify D3h proper rotation angle {angle:.8f}")

    proper = -op
    angle = rotation_angle(proper)
    if abs(angle - math.pi / 3.0) < 1e-4:
        return "2S3"
    if abs(angle - math.pi) < 1e-4:
        rot_axis = operation_axis(proper)
        if abs(float(np.dot(rot_axis, axis))) > 0.9:
            return "sh"
        if abs(float(np.dot(rot_axis, axis))) < 0.1:
            return "3sv"
    raise ValueError(f"Cannot classify D3h improper operation (angle of -R {angle:.8f})")


def signed_planar_angle(op: np.ndarray, axis: np.ndarray, ref: np.ndarray) -> float:
    x = unit_vector(ref - np.dot(ref, axis) * axis)
    y = np.cross(axis, x)
    return math.atan2(float(np.dot(y, op @ x)), float(np.dot(x, op @ x)))


def improper_principal_axis(operations: list[np.ndarray], angle: float) -> np.ndarray:
    # For an S_n operation, -R is a proper rotation about the same principal
    # axis.  Project small numerical errors out before inspecting its angle.
    for op in operations:
        if np.linalg.det(op) >= 0.0:
            continue
        proper = -op
        left, _, right = np.linalg.svd(proper)
        proper = left @ right
        if np.linalg.det(proper) > 0.0 and abs(rotation_angle(proper) - angle) < 1e-3:
            return operation_axis(proper)

    # An S4 matrix has one real eigenvalue -1 along its axis and a complex
    # conjugate pair in the perpendicular plane.  This is independent of the
    # sign convention used for the rotation angle.
    for op in operations:
        if np.linalg.det(op) >= 0.0:
            continue
        values, vectors = np.linalg.eig(op)
        complex_pair = np.count_nonzero(np.abs(np.imag(values)) > 1e-5) == 2
        idx = int(np.argmin(np.abs(values + 1.0)))
        if complex_pair and abs(values[idx] + 1.0) < 1e-3:
            return unit_vector(np.real(vectors[:, idx]))

    determinants = [float(np.linalg.det(op)) for op in operations]
    raise ValueError(
        "Could not infer the principal improper-rotation axis from "
        f"{len(operations)} operations (determinants: "
        + ", ".join(f"{value:.3f}" for value in determinants)
        + ")"
    )


def classify_c3h_operation(op: np.ndarray, axis: np.ndarray, ref: np.ndarray) -> str:
    if np.linalg.det(op) > 0.0:
        angle = rotation_angle(op)
        if angle < 1e-5:
            return "E"
        signed = signed_planar_angle(op, axis, ref)
        return "C3+" if signed > 0.0 else "C3-"
    if np.linalg.norm(op - (np.eye(3) - 2.0 * np.outer(axis, axis))) < 1e-4:
        return "sh"
    signed = signed_planar_angle(op, axis, ref)
    return "S3+" if signed > 0.0 else "S3-"


def classify_s4_operation(op: np.ndarray, axis: np.ndarray, ref: np.ndarray) -> str:
    if np.linalg.det(op) > 0.0:
        return "E" if rotation_angle(op) < 1e-5 else "C2"
    signed = signed_planar_angle(op, axis, ref)
    return "S4+" if signed > 0.0 else "S4-"


def classify_d2d_operation(op: np.ndarray, axis: np.ndarray) -> str:
    det = np.linalg.det(op)
    angle = rotation_angle(op if det > 0.0 else -op)
    if det > 0.0:
        if angle < 1e-5:
            return "E"
        rot_axis = operation_axis(op)
        return "C2" if abs(float(np.dot(rot_axis, axis))) > 0.9 else "2C2p"
    if abs(angle - math.pi / 2.0) < 1e-4:
        return "2S4"
    return "2sd"


def classify_d3d_operation(op: np.ndarray, axis: np.ndarray) -> str:
    det = np.linalg.det(op)
    proper = op if det > 0.0 else -op
    angle = rotation_angle(proper)
    if angle < 1e-5:
        name = "E"
    elif abs(angle - 2.0 * math.pi / 3.0) < 1e-4:
        name = "2C3"
    elif abs(angle - math.pi) < 1e-4 and abs(float(np.dot(operation_axis(proper), axis))) < 0.1:
        name = "3C2p"
    else:
        raise ValueError(f"Cannot classify D3d operation angle {angle:.8f}")
    if det > 0.0:
        return name
    return {"E": "i", "2C3": "2S6", "3C2p": "3sd"}[name]


def classify_d4h_operation(op: np.ndarray, axis: np.ndarray, ref: np.ndarray) -> str:
    det = np.linalg.det(op)
    proper = op if det > 0.0 else -op
    angle = rotation_angle(proper)
    if angle < 1e-5:
        name = "E"
    elif abs(angle - math.pi / 2.0) < 1e-4:
        name = "2C4"
    elif abs(angle - math.pi) < 1e-4:
        rot_axis = operation_axis(proper)
        if abs(float(np.dot(rot_axis, axis))) > 0.9:
            name = "C2"
        else:
            phi = angle_about_axis(rot_axis, ref, axis)
            sector = (phi / (math.pi / 4.0)) % 2.0
            name = "2C2p" if min(abs(sector), abs(sector - 2.0)) < 0.5 else "2C2pp"
    else:
        raise ValueError(f"Cannot classify D4h operation angle {angle:.8f}")
    if det > 0.0:
        return name
    return {"E":"i", "2C4":"2S4", "C2":"sh", "2C2p":"2sv", "2C2pp":"2sd"}[name]


def classify_td_operation(op: np.ndarray) -> str:
    det = np.linalg.det(op)
    angle = rotation_angle(op if det > 0.0 else -op)
    if det > 0.0:
        if angle < 1e-5:
            return "E"
        if abs(angle - 2.0 * math.pi / 3.0) < 1e-4:
            return "8C3"
        if abs(angle - math.pi) < 1e-4:
            return "3C2"
    else:
        if abs(angle - math.pi / 2.0) < 1e-4:
            return "6S4"
        if abs(angle - math.pi) < 1e-4:
            return "6sd"
    raise ValueError(f"Cannot classify Td operation (angle {angle:.8f})")


def classify_oh_operation(op: np.ndarray, c4_axes: list[np.ndarray]) -> str:
    det = np.linalg.det(op)
    proper = op if det > 0.0 else -op
    angle = rotation_angle(proper)
    if angle < 1e-5:
        name = "E"
    elif abs(angle - 2.0 * math.pi / 3.0) < 1e-4:
        name = "8C3"
    elif abs(angle - math.pi / 2.0) < 1e-4:
        name = "6C4"
    elif abs(angle - math.pi) < 1e-4:
        axis = operation_axis(proper)
        name = "3C2" if any(abs(float(np.dot(axis, old))) > 0.9 for old in c4_axes) else "6C2p"
    else:
        raise ValueError(f"Cannot classify Oh operation angle {angle:.8f}")
    if det > 0.0:
        return name
    return {"E":"i", "8C3":"8S6", "6C2p":"6sd", "6C4":"6S4", "3C2":"3sh"}[name]


def axis_from_analyzer(analyzer, coords: np.ndarray, order: int | None = None) -> np.ndarray:
    if order is not None:
        for candidate, candidate_order in analyzer.rot_sym:
            if candidate_order == order:
                return unit_vector(np.array(candidate, dtype=float))
    if analyzer.rot_sym:
        candidate = max(analyzer.rot_sym, key=lambda item: item[1])[0]
        return unit_vector(np.array(candidate, dtype=float))
    return infer_d6_axis(coords)


def operation_axis(op: np.ndarray) -> np.ndarray:
    vals, vecs = np.linalg.eig(op)
    idx = int(np.argmin(np.abs(vals - 1.0)))
    return unit_vector(np.real(vecs[:, idx]))


def reflection_normal(op: np.ndarray) -> np.ndarray:
    vals, vecs = np.linalg.eig(op)
    idx = int(np.argmin(np.abs(vals + 1.0)))
    return unit_vector(np.real(vecs[:, idx]))


def classify_cn_proper(op: np.ndarray, nfold: int) -> str:
    angle = rotation_angle(op)
    if angle < 1e-5:
        return "E"
    if nfold == 2 and abs(angle - math.pi) < 1e-4:
        return "C2"
    if nfold == 3 and abs(angle - 2.0 * math.pi / 3.0) < 1e-4:
        return "2C3"
    if nfold == 4:
        if abs(angle - math.pi / 2.0) < 1e-4:
            return "2C4"
        if abs(angle - math.pi) < 1e-4:
            return "C2"
    if nfold == 6:
        if abs(angle - math.pi / 3.0) < 1e-4:
            return "2C6"
        if abs(angle - 2.0 * math.pi / 3.0) < 1e-4:
            return "2C3"
        if abs(angle - math.pi) < 1e-4:
            return "C2"
    raise ValueError(f"Cannot classify C{nfold}v proper rotation angle {angle:.8f}")


def classify_cnv_reflection(op: np.ndarray, pg: str, axis: np.ndarray, ref: np.ndarray) -> str:
    if pg == "C3v":
        return "3sv"
    normal = reflection_normal(op)
    phi = angle_about_axis(normal, ref, axis)
    if pg == "C2v":
        return "sv" if abs(math.cos(phi)) < 0.5 else "svp"
    sector = (phi / (math.pi / int(pg[1]))) % 2.0
    return f"{int(pg[1]) // 2}sv" if abs(math.sin(phi)) < 0.5 else f"{int(pg[1]) // 2}sd"


def analyzer_operation_matrices(analyzer) -> list[np.ndarray]:
    """Return the completed Cartesian point group from a pymatgen analyzer."""
    cached = getattr(analyzer, "_mo_irreps_operation_matrices", None)
    if cached is not None:
        return cached
    candidates = []
    for getter in (
        lambda: analyzer.get_pointgroup(),
        lambda: analyzer.get_symmetry_operations(),
        lambda: analyzer.symmops,
        lambda: analyzer._symmops,
    ):
        try:
            value = getter()
            value = getattr(value, "symmetry_ops", value)
            matrices = [np.array(item.rotation_matrix, dtype=float) for item in value]
            if matrices:
                candidates.append(matrices)
        except (AttributeError, TypeError):
            continue
    if not candidates:
        raise ValueError("pymatgen returned no molecular symmetry operations")

    # get_symmetry_operations() can return only generators for some pymatgen
    # versions. Complete each candidate set under matrix multiplication.
    completed = []
    for generators in candidates:
        group = [np.eye(3)]
        pending = list(generators)
        while pending:
            matrix = pending.pop()
            left, _, right = np.linalg.svd(matrix)
            matrix = left @ right
            if any(np.allclose(matrix, old, atol=1e-6) for old in group):
                continue
            group.append(matrix)
            if len(group) > 240:
                raise ValueError("Symmetry-operation closure unexpectedly exceeds 240 elements")
            pending.extend(matrix @ old for old in group)
            pending.extend(old @ matrix for old in group)
        completed.append(group)
    group = max(completed, key=len)

    # Some pymatgen releases label an S4 molecule correctly but expose only
    # its proper C2 subgroup.  Reconstruct the two S4 operations from that
    # unique C2 axis; their sense is immaterial because both are included.
    if analyzer.sch_symbol == "S4" and not any(np.linalg.det(op) < 0.0 for op in group):
        c2_ops = [
            op
            for op in group
            if np.linalg.det(op) > 0.0 and abs(rotation_angle(op) - math.pi) < 1e-4
        ]
        if len(c2_ops) == 1:
            axis = operation_axis(c2_ops[0])
            cross = np.array(
                [
                    [0.0, -axis[2], axis[1]],
                    [axis[2], 0.0, -axis[0]],
                    [-axis[1], axis[0], 0.0],
                ]
            )
            reflection = np.eye(3) - 2.0 * np.outer(axis, axis)
            for signed_angle in (math.pi / 2.0, -math.pi / 2.0):
                rotation = (
                    math.cos(signed_angle) * np.eye(3)
                    + math.sin(signed_angle) * cross
                    + (1.0 - math.cos(signed_angle)) * np.outer(axis, axis)
                )
                group.append(reflection @ rotation)
    try:
        analyzer._mo_irreps_operation_matrices = group
    except AttributeError:
        pass
    return group


def d2_frame(analyzer) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if hasattr(analyzer, "principal_axes"):
        axes = [unit_vector(np.array(axis, dtype=float)) for axis in analyzer.principal_axes]
        if len(axes) >= 3:
            return axes[2], axes[1], axes[0]
    proper_axes = []
    for mat in analyzer_operation_matrices(analyzer):
        if np.linalg.det(mat) > 0.0 and abs(rotation_angle(mat) - math.pi) < 1e-4:
            axis = operation_axis(mat)
            if all(abs(float(np.dot(axis, old))) < 0.9 for old in proper_axes):
                proper_axes.append(axis)
    if len(proper_axes) < 3:
        raise ValueError("Could not infer D2 axes from pymatgen operations")
    return proper_axes[0], proper_axes[1], proper_axes[2]


def closest_d2_axis_name(axis: np.ndarray, frame: tuple[np.ndarray, np.ndarray, np.ndarray]) -> str:
    names = ["C2z", "C2y", "C2x"]
    dots = [abs(float(np.dot(axis, candidate))) for candidate in frame]
    return names[int(np.argmax(dots))]


def classify_operation(pg: str, op: np.ndarray, analyzer, coords: np.ndarray) -> str:
    det = np.linalg.det(op)
    if pg == "C3h":
        axis = axis_from_analyzer(analyzer, coords, 3)
        return classify_c3h_operation(op, axis, infer_c2_reference(coords, axis))
    if pg in ("S4", "D2d"):
        matrices = analyzer_operation_matrices(analyzer)
        axis = improper_principal_axis(matrices, math.pi / 2.0)
        if pg == "S4":
            return classify_s4_operation(op, axis, infer_c2_reference(coords, axis))
        return classify_d2d_operation(op, axis)
    if pg == "D3d":
        return classify_d3d_operation(op, axis_from_analyzer(analyzer, coords, 3))
    if pg == "D3h":
        axis = axis_from_analyzer(analyzer, coords, 3)
        return classify_d3h_operation(op, axis)
    if pg == "D4h":
        axis = axis_from_analyzer(analyzer, coords, 4)
        return classify_d4h_operation(op, axis, infer_c2_reference(coords, axis))
    if pg == "D6h":
        axis = axis_from_analyzer(analyzer, coords, 6)
        ref = infer_c2_reference(coords, axis)
        return classify_d6h_operation(op, axis, ref)
    if pg == "Td":
        return classify_td_operation(op)
    if pg == "Oh":
        c4_axes = [
            operation_axis(matrix)
            for matrix in analyzer_operation_matrices(analyzer)
            if np.linalg.det(matrix) > 0.0
            and abs(rotation_angle(matrix) - math.pi / 2.0) < 1e-4
        ]
        return classify_oh_operation(op, c4_axes)
    if det > 0.0 and rotation_angle(op) < 1e-5:
        return "E"
    if pg == "C1":
        return "E"
    if pg == "Ci":
        return "i"
    if pg == "Cs":
        return "s"
    if pg in ("C2", "C2h"):
        if det > 0.0:
            return "C2"
        return "i" if np.linalg.norm(op + np.eye(3)) < 1e-5 else "sh"
    if pg in ("C2v", "C3v", "C4v", "C6v"):
        nfold = int(pg[1])
        axis = axis_from_analyzer(analyzer, coords, nfold)
        ref = infer_c2_reference(coords, axis)
        if det > 0.0:
            return classify_cn_proper(op, nfold)
        return classify_cnv_reflection(op, pg, axis, ref)
    if pg in ("D2", "D2h"):
        frame = d2_frame(analyzer)
        if det > 0.0:
            return closest_d2_axis_name(operation_axis(op), frame)
        if np.linalg.norm(op + np.eye(3)) < 1e-5:
            return "i"
        proper_name = closest_d2_axis_name(operation_axis(-op), frame)
        return {"C2z": "sxy", "C2y": "sxz", "C2x": "syz"}[proper_name]
    raise ValueError(f"Point group {pg} is not supported")


def subgroup_hint(pg: str) -> str:
    suggestions = [item for item in SUBGROUP_SUGGESTIONS.get(pg, []) if item in POINT_GROUPS]
    if not suggestions:
        suggestions = ["C1"]
    return (
        f"Point group {pg} is not supported yet. Supported groups: "
        f"{', '.join(sorted(POINT_GROUPS))}. "
        f"Possible subgroup analyses: {', '.join(suggestions)}. "
        "Use --point-group <subgroup> to request one explicitly."
    )


def select_subgroup_operations(
    requested_pg: str, analyzer, coords: np.ndarray
) -> list[tuple[str, np.ndarray]]:
    target_counts = dict(POINT_GROUPS[requested_pg]["classes"])
    selected: list[tuple[str, np.ndarray]] = []
    counts = defaultdict(int)

    for mat in analyzer_operation_matrices(analyzer):
        try:
            class_name = classify_operation(requested_pg, mat, analyzer, coords)
        except Exception:
            continue
        if class_name not in target_counts:
            continue
        if counts[class_name] >= target_counts[class_name]:
            continue
        selected.append((class_name, mat))
        counts[class_name] += 1

    missing = [
        f"{name}({target_counts[name] - counts[name]})"
        for name in target_counts
        if counts[name] != target_counts[name]
    ]
    if missing:
        raise ValueError(
            f"Could not construct requested {requested_pg} subgroup from detected symmetry "
            f"operations; missing classes: {', '.join(missing)}"
        )
    return selected


def pymatgen_point_group_operations(
    species: list[str], coords: np.ndarray, tolerance: float, requested_pg: str
) -> tuple[str, list[tuple[str, np.ndarray]]]:
    try:
        from pymatgen.core import Molecule
        from pymatgen.symmetry.analyzer import PointGroupAnalyzer
    except ImportError as exc:
        raise RuntimeError("pymatgen is not installed") from exc

    analyzer = PointGroupAnalyzer(Molecule(species, coords), tolerance=tolerance)
    pg = analyzer.sch_symbol
    if requested_pg != "auto" and pg != requested_pg:
        if requested_pg not in POINT_GROUPS:
            raise ValueError(f"Requested point group {requested_pg} is not supported")
        ops = select_subgroup_operations(requested_pg, analyzer, coords)
        print(
            f"mo_irreps.py: warning: pymatgen detected {pg}; analyzing requested "
            f"{requested_pg} subgroup instead.",
            file=sys.stderr,
        )
        return requested_pg, ops
    if pg not in POINT_GROUPS:
        raise ValueError(subgroup_hint(pg))

    ops = []
    for mat in analyzer_operation_matrices(analyzer):
        ops.append((classify_operation(pg, mat, analyzer, coords), mat))
    return pg, ops


def polynomial_value(poly: dict[tuple[int, int, int], float], points: np.ndarray) -> np.ndarray:
    out = np.zeros(len(points))
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    for (ix, iy, iz), coef in poly.items():
        out += coef * (x**ix) * (y**iy) * (z**iz)
    return out


def shell_transform(shell: str, op: np.ndarray) -> np.ndarray:
    labels = SHELLS[shell]
    if shell == "s":
        return np.ones((1, 1))
    grid = np.array(
        [
            [-1.3, -0.7, 0.2],
            [-0.4, 1.1, -1.2],
            [0.8, -1.5, 0.6],
            [1.4, 0.9, -0.3],
            [-1.1, 1.2, 1.0],
            [0.3, -0.2, -1.4],
            [1.7, -0.8, 1.3],
            [-0.6, -1.6, -0.9],
            [0.9, 1.5, 0.7],
        ],
        dtype=float,
    )
    basis = np.column_stack([polynomial_value(POLYS[label], grid) for label in labels])
    transformed_points = grid @ op
    transform = np.zeros((len(labels), len(labels)))
    for j, label in enumerate(labels):
        vals = polynomial_value(POLYS[label], transformed_points)
        transform[:, j], *_ = np.linalg.lstsq(basis, vals, rcond=None)
    return transform


def atom_permutation(species: list[str], coords: np.ndarray, op: np.ndarray, tol: float) -> list[int]:
    rotated = coords @ op.T
    used: set[int] = set()
    perm = [-1] * len(species)
    for old, point in enumerate(rotated):
        distances = np.linalg.norm(coords - point, axis=1)
        candidates = np.argsort(distances)
        for new in candidates:
            if new not in used and species[new] == species[old] and distances[new] <= tol:
                perm[old] = int(new)
                used.add(int(new))
                break
        if perm[old] < 0:
            raise ValueError(
                f"Could not map atom {old + 1} under symmetry operation; "
                f"nearest distance is {distances[candidates[0]]:.3e}"
            )
    return perm


def build_ao_operation(
    aos: list[AO], species: list[str], coords: np.ndarray, op: np.ndarray, tolerance: float
) -> np.ndarray:
    perm = atom_permutation(species, coords, op, tolerance)
    shell_blocks = {ao.shell: shell_transform(ao.shell, op) for ao in aos}
    by_key = {(ao.atom, ao.shell): [] for ao in aos}
    for idx, ao in enumerate(aos):
        by_key[(ao.atom, ao.shell)].append(idx)
    matrix = np.zeros((len(aos), len(aos)))
    for atom in range(1, len(species) + 1):
        new_atom = perm[atom - 1] + 1
        shells = {ao.shell for ao in aos if ao.atom == atom}
        for shell in shells:
            old_idxs = by_key[(atom, shell)]
            new_idxs = by_key[(new_atom, shell)]
            old_labels = [aos[i].label for i in old_idxs]
            new_labels = [aos[i].label for i in new_idxs]
            if old_labels != SHELLS[shell] or new_labels != SHELLS[shell]:
                raise ValueError(
                    f"AO order for atom {atom} shell {shell} is not the expected DFTB+ order"
                )
            block = shell_blocks[shell]
            matrix[np.ix_(new_idxs, old_idxs)] = block
    return matrix


def group_blocks_from_band(
    band: list[BandState], eigvecs: list[Eigenvector], energy_tol: float
) -> list[tuple[int, int, list[int], float | None, float | None]]:
    by_vec = {(vec.kpt, spin_number(vec.spin), vec.state): idx for idx, vec in enumerate(eigvecs)}
    blocks = []
    for key, states in groupby_key(band, lambda item: (item.kpt, item.spin)).items():
        current: list[int] = []
        current_energy: float | None = None
        current_occ = 0.0
        for state in states:
            idx = by_vec.get((state.kpt, state.spin, state.state))
            if idx is None:
                continue
            if current and abs(state.energy - (current_energy or state.energy)) > energy_tol:
                blocks.append((key[0], key[1], current, current_energy, current_occ))
                current = []
                current_occ = 0.0
            current.append(idx)
            current_energy = state.energy
            current_occ += state.occ
        if current:
            blocks.append((key[0], key[1], current, current_energy, current_occ))
    return blocks


def groupby_key(items, keyfunc):
    grouped = defaultdict(list)
    for item in items:
        grouped[keyfunc(item)].append(item)
    return grouped


def spin_number(spin: str) -> int:
    words = spin.lower().split()
    if words and words[-1].isdigit():
        return int(words[-1])
    if "down" in words:
        return 2
    return 1


def default_blocks(eigvecs: list[Eigenvector]) -> list[tuple[int, int, list[int], None, None]]:
    return [(vec.kpt, spin_number(vec.spin), [idx], None, None) for idx, vec in enumerate(eigvecs)]


def class_average(point_group: str, chars_by_op: list[tuple[str, complex]]) -> np.ndarray:
    values = []
    grouped = groupby_key(chars_by_op, lambda item: item[0])
    for name, count in POINT_GROUPS[point_group]["classes"]:
        vals = [char for _, char in grouped[name]]
        if len(vals) != count:
            raise ValueError(f"Internal error: class {name} has {len(vals)} operations, expected {count}")
        values.append(sum(vals) / count)
    return np.array(values, dtype=complex)


def decompose(point_group: str, characters: np.ndarray) -> dict[str, float]:
    classes = POINT_GROUPS[point_group]["classes"]
    irreps = POINT_GROUPS[point_group]["characters"]
    order = sum(count for _, count in classes)
    out = {}
    for name, chars in irreps.items():
        weighted = sum(
            count * np.conjugate(ref) * got
            for (_, count), ref, got in zip(classes, chars, characters)
        )
        mult = weighted / order
        if abs(mult.real) > 0.15:
            out[name] = mult.real
    return out


def normalize_coeffs(coeffs: np.ndarray, overlap: np.ndarray | None) -> tuple[np.ndarray, float]:
    if overlap is None:
        norm = np.vdot(coeffs, coeffs).real
    else:
        norm = np.vdot(coeffs, overlap @ coeffs).real
    if norm <= 0.0:
        return coeffs, norm
    return coeffs / math.sqrt(norm), norm


def format_mult(mult: dict[str, float], point_group: str) -> str:
    if not mult:
        return "unassigned"
    mult = dict(mult)
    conjugate_pairs = []
    if point_group == "C3h":
        conjugate_pairs = [("E1'", "E2'", "E'"), ('E1"', 'E2"', 'E"')]
    elif point_group == "S4":
        conjugate_pairs = [("E+", "E-", "E")]
    for first, second, combined in conjugate_pairs:
        if first in mult and second in mult and abs(mult[first] - mult[second]) < 0.2:
            mult[combined] = 0.5 * (mult.pop(first) + mult.pop(second))
    parts = []
    for name, value in sorted(mult.items()):
        rounded = int(round(value))
        if abs(value - rounded) < 0.2:
            if rounded == 1:
                parts.append(name)
            elif rounded:
                parts.append(f"{rounded}{name}")
        else:
            parts.append(f"{value:.2f}{name}")
    return " + ".join(parts) if parts else "unassigned"


def representation_characters(point_group: str, mult: dict[str, float]) -> np.ndarray:
    table = POINT_GROUPS[point_group]["characters"]
    chars = np.zeros(len(POINT_GROUPS[point_group]["classes"]), dtype=complex)
    for name, value in mult.items():
        chars += value * np.asarray(table[name], dtype=complex)
    return chars


def transition_product(
    point_group: str,
    occupied: dict[str, float],
    virtual: dict[str, float],
) -> dict[str, float]:
    occupied_chars = representation_characters(point_group, occupied)
    virtual_chars = representation_characters(point_group, virtual)
    return decompose(point_group, np.conjugate(occupied_chars) * virtual_chars)


def annotate_exc_dat(
    input_path: Path,
    output_path: Path,
    point_group: str,
    mo_irreps: dict[tuple[int, int, int], dict[str, float]],
) -> None:
    """Append dominant occupied-virtual product symmetries to EXC.DAT."""
    if input_path.resolve() == output_path.resolve():
        raise ValueError("Excitation input and annotated output paths must be different")

    transition = re.compile(r"(?P<occ>\d+)\s*->\s*(?P<virt>\d+)")
    output = []
    annotated = 0
    for line in input_path.read_text().splitlines():
        match = transition.search(line)
        if match is None:
            if "Transition" in line and "Spatial irrep" not in line:
                line += "      Spatial irrep (dominant MO product)"
            output.append(line)
            continue

        occ = int(match.group("occ"))
        virt = int(match.group("virt"))
        occ_irrep = mo_irreps.get((1, 1, occ))
        virt_irrep = mo_irreps.get((1, 1, virt))
        if occ_irrep is None or virt_irrep is None:
            label = "unassigned"
        else:
            label = format_mult(
                transition_product(point_group, occ_irrep, virt_irrep), point_group
            )
        output.append(f"{line}      {label}")
        annotated += 1

    if annotated == 0:
        raise ValueError(f"{input_path}: no occupied -> virtual transitions found")
    output_path.write_text("\n".join(output) + "\n")


def warn_if_orientation_convention_matters(point_group: str) -> None:
    detail = ORIENTATION_WARNINGS.get(point_group)
    if detail is None:
        return
    print(
        f"mo_irreps.py: warning: some {point_group} irrep names depend on the "
        f"symmetry-axis convention: {detail}. This is only a relabeling; the "
        "physical symmetry decomposition is unchanged.",
        file=sys.stderr,
    )


def orientation_description(
    point_group: str,
    operations: list[tuple[str, np.ndarray]],
    coords: np.ndarray,
) -> list[str]:
    """Describe the Cartesian directions used to name symmetry operations."""

    def fmt(vec: np.ndarray) -> str:
        vec = unit_vector(np.asarray(vec, dtype=float))
        # Axis signs are physically equivalent; choose one for stable output.
        significant = np.flatnonzero(np.abs(vec) > 1e-10)
        if len(significant) and vec[significant[0]] < 0.0:
            vec = -vec
        return "[" + ", ".join(f"{value:+.8f}" for value in vec) + "]"

    by_class = groupby_key(operations, lambda item: item[0])
    lines = [
        "orientation: molecule was not rotated; vectors below are in the input XYZ frame"
    ]

    if point_group in ("C2", "C2h", "C2v", "S4", "D2d"):
        c2 = by_class.get("C2", [])
        if c2:
            lines.append(f"orientation C2/z axis: {fmt(operation_axis(c2[0][1]))}")
        if point_group == "C2v" and c2:
            axis = operation_axis(c2[0][1])
            lines.append(
                "orientation geometry-selected in-plane reference: "
                f"{fmt(infer_c2_reference(coords, axis))}"
            )
    elif point_group in ("C3h", "C3v", "C4v", "C6v", "D3d", "D3h", "D4h", "D6h"):
        class_name = {
            "C3h": "C3+",
            "C3v": "2C3",
            "C4v": "2C4",
            "C6v": "2C6",
            "D3d": "2C3",
            "D3h": "2C3",
            "D4h": "2C4",
            "D6h": "2C6",
        }[point_group]
        rotations = by_class.get(class_name, [])
        if rotations:
            axis = operation_axis(rotations[0][1])
            lines.append(f"orientation principal/z axis: {fmt(axis)}")
            if point_group not in ("D3d", "D3h"):
                lines.append(
                    "orientation geometry-selected in-plane reference: "
                    f"{fmt(infer_c2_reference(coords, axis))}"
                )
    elif point_group in ("D2", "D2h"):
        for class_name, label in (
            ("C2x", "x"),
            ("C2y", "y"),
            ("C2z", "z"),
        ):
            rotations = by_class.get(class_name, [])
            if rotations:
                axis = operation_axis(rotations[0][1])
                lines.append(f"orientation {label} axis ({class_name}): {fmt(axis)}")
    elif point_group == "Cs":
        mirrors = by_class.get("s", [])
        if mirrors:
            normal = reflection_normal(mirrors[0][1])
            lines.append(f"orientation mirror-plane normal: {fmt(normal)}")

    if len(lines) == 1:
        lines.append("orientation no convention-dependent symmetry axis")
    return lines


def analyze(args: argparse.Namespace) -> int:
    for path, description in (
        (args.eigenvec, "eigenvector file"),
        (args.geometry, "geometry file"),
        (args.band, "band.out file"),
        (args.overlap, "oversqr.dat file"),
    ):
        if not path.exists():
            raise ValueError(f"Required {description} not found: {path}")

    species, coords = parse_xyz(args.geometry)
    aos, eigvecs = parse_eigenvec(args.eigenvec)
    if len(aos) != len(eigvecs[0].coeffs):
        raise ValueError("AO metadata length does not match coefficient length")
    if len(species) != max(ao.atom for ao in aos):
        raise ValueError("Geometry atom count does not match eigenvec.out atom numbering")

    overlaps = parse_oversqr(args.overlap)
    band = parse_band(args.band)
    blocks = group_blocks_from_band(band, eigvecs, args.energy_tol)

    point_group, operations = pymatgen_point_group_operations(
        species, coords, args.tol, args.point_group
    )
    warn_if_orientation_convention_matters(point_group)
    ao_ops = [
        (name, build_ao_operation(aos, species, coords, op, args.tol))
        for name, op in operations
    ]

    print(f"# MO irreps for point group {point_group}")
    for line in orientation_description(point_group, operations, coords):
        print(f"# {line}")
    print("# columns: kpt spin states energy occ irrep norm_error")
    metric_ops_by_kpt: dict[int, list[tuple[str, np.ndarray]]] = {}
    mo_irreps: dict[tuple[int, int, int], dict[str, float]] = {}
    for kpt, spin, indices, energy, occ in blocks:
        overlap = overlaps.get(kpt)
        coeff_cols = []
        norm_errors = []
        for idx in indices:
            coeff, norm = normalize_coeffs(eigvecs[idx].coeffs, overlap)
            coeff_cols.append(coeff)
            norm_errors.append(abs(norm - 1.0) if norm > 0.0 else float("nan"))
        coeffs = np.column_stack(coeff_cols)
        chars_by_op = []
        metric = overlap if overlap is not None else np.eye(len(aos))
        if kpt not in metric_ops_by_kpt:
            metric_ops_by_kpt[kpt] = [
                (class_name, metric @ ao_op) for class_name, ao_op in ao_ops
            ]
        for class_name, metric_op in metric_ops_by_kpt[kpt]:
            rep = coeffs.conj().T @ metric_op @ coeffs
            chars_by_op.append((class_name, np.trace(rep)))
        chars = class_average(point_group, chars_by_op)
        mult = decompose(point_group, chars)
        for idx in indices:
            vec = eigvecs[idx]
            mo_irreps[(vec.kpt, spin_number(vec.spin), vec.state)] = mult
        states = ",".join(str(eigvecs[idx].state) for idx in indices)
        energy_txt = f"{energy:.8f}" if energy is not None else "-"
        occ_txt = f"{occ:.5f}" if occ is not None else "-"
        norm_txt = f"{max(norm_errors):.2e}" if norm_errors else "-"
        print(
            f"{kpt:3d} {spin:3d} {states:>9s} {energy_txt:>14s} {occ_txt:>9s} "
            f"{format_mult(mult, point_group):>14s} {norm_txt:>10s}"
        )
    if args.exc_input.exists():
        annotate_exc_dat(
            args.exc_input, args.exc_output, point_group, mo_irreps
        )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Determine molecular-orbital irreps from DFTB+ eigenvec.out.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Required inputs:\n"
            "  eigenvec.out  MO coefficients and AO order, from WriteEigenvectors + EigenvectorsAsText.\n"
            "  geometry      XYZ geometry matching the atom order in eigenvec.out.\n"
            "  band.out      Energies/occupations; used to group degenerate orbitals.\n"
            "  oversqr.dat   AO overlap matrix; used for S-normalization and projections.\n"
            "  EXC.DAT       Optional excitation list; its dominant MO products are annotated.\n"
            "  pymatgen      Required for point-group detection and symmetry operations.\n\n"
            "Tolerance tips:\n"
            "  --tol controls both pymatgen point-group detection and the atom-mapping\n"
            "    check used to construct AO symmetry operations. Increase it for a\n"
            "    slightly distorted or rounded geometry.\n"
            "  --energy-tol controls grouping of nearly-degenerate levels from band.out.\n"
            "    Increase it if a true degenerate pair is printed as separate 0.5E-like pieces.\n\n"
            "Unsupported high-symmetry groups:\n"
            "  The script stops and suggests supported subgroup analyses. Use --point-group\n"
            "  explicitly to request one of those subgroups; no subgroup fallback is automatic."
        ),
    )
    parser.add_argument(
        "-e",
        "--eigenvec",
        type=Path,
        default=Path("eigenvec.out"),
        help="DFTB+ text eigenvector file with MO coefficients and AO order (required)",
    )
    parser.add_argument(
        "-g",
        "--geometry",
        type=Path,
        default=Path("geo.xyz"),
        help="XYZ geometry matching eigenvec.out atom order (required)",
    )
    parser.add_argument(
        "-b",
        "--band",
        type=Path,
        default=Path("band.out"),
        help="DFTB+ band.out file for energies, occupations, and degeneracy grouping (required)",
    )
    parser.add_argument(
        "-s",
        "--overlap",
        type=Path,
        default=Path("oversqr.dat"),
        help="DFTB+ square overlap matrix for normalization/projections (required)",
    )
    parser.add_argument(
        "--exc-input",
        type=Path,
        default=Path("EXC.DAT"),
        help="DFTB+ excitation file to annotate when present (default: EXC.DAT)",
    )
    parser.add_argument(
        "--exc-output",
        type=Path,
        default=Path("EXC.DAT.SYM"),
        help="annotated excitation output file (default: EXC.DAT.SYM)",
    )
    parser.add_argument(
        "--point-group",
        choices=["auto"] + sorted(POINT_GROUPS),
        default="auto",
        help="Point group to use; auto detects with pymatgen",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=0.1,
        help="geometry tolerance for point-group detection and atom mapping",
    )
    parser.add_argument(
        "--energy-tol",
        type=float,
        default=1e-5,
        help="energy tolerance for grouping degenerate states from band.out",
    )
    try:
        return analyze(parser.parse_args(argv))
    except Exception as exc:
        print(f"mo_irreps.py: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
