from random import randint

import sympy as sp
from dtransform import Spectrum


EQUATIONS: tuple = (
    ("x + y", "1 + x * y", {'x': 1, 'y': 2}, None),
    ("2 + x + 3 * y", "1 + x * y", {'x': 1, 'y': 2}, None),
    ("2 + x + y", "1 + x * y", {'x': 1, 'y': 2}, None),
    ("pi + x + y", "1 + x * y", {'x': 1, 'y': 2}, None),
    ("pi + x + y", "1 + x * y", {'x': 0, 'y': 0}, None),
    ("pi + x + y", "1 + x * y", {'x': 1, 'y': 2}, {'x': 0.5, 'y': 1.5}),
    ("(1 + x) / y", "1 + x - y", {'x': 1, 'y': 1}, None),
    ("sin(x) / y", "1 - x - y", {'x': 1, 'y': 2}, None),
)


def test_addition() -> None:
    global EQUATIONS

    for f1, f2, center, scaling in EQUATIONS:
        s1 = Spectrum(f1, order=3, center=center, scaling=scaling)
        s2 = Spectrum(f2, order=3, center=center, scaling=scaling)
        assert (s1 + s2).inverse().evalf(subs=center) == (
            sp.sympify(f1) + sp.sympify(f2)
        ).evalf(subs=center)


def test_subtraction() -> None:
    global EQUATIONS

    for f1, f2, center, scaling in EQUATIONS:
        s1 = Spectrum(f1, order=3, center=center, scaling=scaling)
        s2 = Spectrum(f2, order=3, center=center, scaling=scaling)
        assert (s1 - s2).inverse().evalf(subs=center) == (
            sp.sympify(f1) - sp.sympify(f2)
        ).evalf(subs=center)


def test_multiplication() -> None:
    global EQUATIONS

    for f1, f2, center, scaling in EQUATIONS:
        s1 = Spectrum(f1, order=3, center=center, scaling=scaling)
        s2 = Spectrum(f2, order=3, center=center, scaling=scaling)
        assert (s1 * s2).inverse().evalf(subs=center) == (
            sp.sympify(f1) * sp.sympify(f2)
        ).evalf(subs=center)


def test_scalar_multiplication() -> None:
    global EQUATIONS

    for f1, f2, center, scaling in EQUATIONS:
        s1 = Spectrum(f1, order=3, center=center, scaling=scaling)
        s2 = Spectrum(f2, order=3, center=center, scaling=scaling)

        scalar = randint(2, 9)

        assert (scalar * s1).inverse().evalf(subs=center) == (
            scalar * sp.sympify(f1)
        ).evalf(subs=center)

        assert (scalar * s2).inverse().evalf(subs=center) == (
            scalar * sp.sympify(f2)
        ).evalf(subs=center)


def test_division() -> None:
    global EQUATIONS

    for f1, f2, center, scaling in EQUATIONS:
        s1 = Spectrum(f1, order=3, center=center, scaling=scaling)
        s2 = Spectrum(f2, order=3, center=center, scaling=scaling)

        assert round((s1 / s2).inverse().evalf(subs=center), 14) == round(
            (sp.sympify(f1) / sp.sympify(f2)).evalf(subs=center),
            14
        )


def test_scalar_division() -> None:
    global EQUATIONS

    for f1, f2, center, scaling in EQUATIONS:
        s1 = Spectrum(f1, order=3, center=center, scaling=scaling)
        s2 = Spectrum(f2, order=3, center=center, scaling=scaling)

        scalar = randint(2, 9)

        assert (s1 / scalar).inverse().evalf(subs=center) == (
            sp.sympify(f1) / scalar
        ).evalf(subs=center)

        assert (s2 / scalar).inverse().evalf(subs=center) == (
            sp.sympify(f2) / scalar
        ).evalf(subs=center)
