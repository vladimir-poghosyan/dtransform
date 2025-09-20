from itertools import product
from typing import Union

import sympy as sp


class Spectrum:

    __slots__ = (
        '__expr',
        '__order',
        '__center',
        '__scaling',
        '__coeffs',
        '__variables',
        '__order_prod',
    )

    def __init__(
        self,
        expr: str,
        order: int = 4,
        center: dict[str, int | float] | None = None,
        scaling: dict[str, int | float] | None = None,
        **kwargs: int | float
    ) -> None:
        if not isinstance(order, int) or order <= 0:
            raise ValueError("Order must be a positive integer")

        self.__order = order
        self.__expr = sp.sympify(expr)
        self.__center: dict = {}
        self.__scaling: dict = {}
        self.__variables = tuple(sorted(
            self.__expr.free_symbols,
            key=lambda s: s.name
        ))

        # populate values or use defaults for center and scaling
        _center = (center or {}) | kwargs
        _scaling = scaling or {}

        for var in self.__variables:
            str_var = str(var)

            self.__center[var] = _center.get(str_var, 0)

            if (scale := _scaling.get(str_var, 1)) <= 0:
                raise ValueError(
                    f"Scaling constant for variable '{var}'"
                    " must be a positive integer"
                )

            self.__scaling[var] = scale

        del _center, _scaling

        # Compute DTM coefficients around center (transformation spectrum)
        coeffs = {}

        for multi_idx in product(*[range(order)] * len(self.__variables)):
            deriv = self.__expr
            for var, k in zip(self.__variables, multi_idx):
                deriv = sp.diff(deriv, var, k)

            denom = sp.prod([sp.factorial(k) for k in multi_idx])
            H_terms = sp.prod([
                self.__scaling[var] ** k
                for var, k in zip(self.__variables, multi_idx)
            ])
            value = deriv.subs(self.__center) * H_terms / denom
            coeffs[multi_idx] = value

        self.__coeffs: dict[tuple[int, ...], sp.Expr] = coeffs

    @property
    def order(self) -> int:
        return self.__order

    @property
    def scaling(self) -> dict:
        return self.__scaling

    @property
    def variables(self) -> tuple:
        return self.__variables

    @property
    def coeffs(self) -> dict:
        return self.__coeffs

    @property
    def center(self) -> dict:
        return self.__center

    def inverse(self) -> sp.Expr:
        """Reconstruct function from transformation spectrum
        (Taylor expansion).

        """
        reconstructed = 0
        for multi_idx, coeff in self.__coeffs.items():
            term = coeff
            for var, power in zip(self.__variables, multi_idx):
                a = self.__center.get(var, 0)
                term *= ((var - a) / self.__scaling[var]) ** power
            reconstructed += term

        return sp.simplify(reconstructed)

    def clone(self) -> 'Spectrum':
        new = object.__new__(Spectrum)
        new.__expr = self.__expr
        new.__order = self.__order
        new.__scaling = self.__scaling
        new.__variables = self.__variables
        new.__center = self.__center
        new.__coeffs = self.__coeffs

        return new

    def display_coefficients(self) -> None:
        """Print non-zero transformation spectrum coefficients."""
        for idx, val in sorted(self.coeffs.items()):
            print(f"Spectrum[{idx}] = {val}")

    def _check_compatibility(self, other: 'Spectrum') -> None:
        if self.__variables != other.variables:
            raise ValueError("Variables do not match.")
        if self.__order != other.order:
            raise ValueError("Expansion order mismatch.")
        if self.__scaling != other.scaling:
            raise ValueError("Scaling constants mismatch.")
        if self.__center != other.center:
            raise ValueError("Expansion center mismatch.")

    def __repr__(self) -> str:
        return (
            f"Spectrum(expr='{str(self.__expr)}', order={self.__order},"
            f" center={self.__center}, scaling={self.__scaling})"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Spectrum):
            return False

        return (
            self.__order == other.order
            and self.__scaling == other.scaling
            and self.__variables == other.variables
            and self.__center == other.center
            and self.__coeffs == other.coeffs
        )

    def __neg__(self) -> 'Spectrum':
        new = self.clone()
        new.__coeffs = {k: -v for k, v in self.__coeffs.items()}
        return new

    def __add__(self, other: 'Spectrum') -> 'Spectrum':
        self._check_compatibility(other)

        new = self.clone()
        new.__coeffs = {
            k: self.__coeffs.get(k, 0) + other.coeffs.get(k, 0)
            for k in (frozenset(self.__coeffs) | frozenset(other.coeffs))
        }

        return new

    def __sub__(self, other: 'Spectrum') -> 'Spectrum':
        self._check_compatibility(other)

        new = self.clone()
        new.__coeffs = {
            k: self.__coeffs.get(k, 0) - other.coeffs.get(k, 0)
            for k in (frozenset(self.__coeffs) | frozenset(other.coeffs))
        }

        return new

    def __mul__(
        self,
        other: Union['Spectrum', int, float, sp.Basic]
    ) -> 'Spectrum':
        if isinstance(other, Spectrum):
            self._check_compatibility(other)

            new = self.clone()
            coeffs = {}

            for idx in product(*[range(self.__order)] * len(self.variables)):
                val = 0
                for i in product(*[range(k + 1) for k in idx]):
                    j = tuple(k - l for k, l in zip(idx, i))
                    val += self.__coeffs.get(i, 0) * other.coeffs.get(j, 0)
                coeffs[idx] = val

            new.__coeffs = coeffs
            return new
        elif isinstance(other, (int, float, sp.Basic)):
            new = self.clone()
            new.__coeffs = {
                k: other * v
                for k, v in self.__coeffs.items()
            }
            return new

        else:
            raise TypeError(
                f'"Spectrum" can\'t be multiplied by {type(other)}'
            )

    def __rmul__(
        self,
        other: Union['Spectrum', int, float, sp.Basic]
    ) -> 'Spectrum':
        return self.__mul__(other)

    def __truediv__(
        self,
        other: Union['Spectrum', int, float, sp.Basic]
    ) -> 'Spectrum':
        if isinstance(other, Spectrum):
            self._check_compatibility(other)

            new = self.clone()
            n_vars = len(self.__variables)
            coeffs: dict[tuple[int, ...], sp.Expr] = {}
            zeros = (0,) * n_vars

            if (zeros_coeff := other.coeffs.get(zeros)) == 0:
                raise ZeroDivisionError(
                    "Leading coefficient of denominator is zero."
                )

            for idx in product(*[range(self.__order)] * n_vars):
                val = self.__coeffs.get(idx, 0)

                for i in product(*[range(k + 1) for k in idx]):
                    j = tuple(k - l for k, l in zip(idx, i))
                    if j != idx:
                        val -= other.coeffs.get(i, 0) * coeffs.get(j, 0)

                coeffs[idx] = val / zeros_coeff

            new.__coeffs = coeffs
            return new
        elif isinstance(other, (int, float, sp.Basic)):
            if other == 0:
                raise ZeroDivisionError("Division by zero.")

            new = self.clone()
            new.__coeffs = {
                k: v / other
                for k, v in self.__coeffs.items()
            }
            return new
        else:
            raise TypeError(
                '"Spectrum" can only be divided by'
                f' another Spectrum or a scalar, not {type(other)}.'
            )
