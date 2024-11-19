from typing import Union, Optional
from visualisation import _plot_curve, _plot_point_addition
from finite_field import FiniteFieldElement

class EllipticCurve:
    """
    A class representing an elliptic curve over a finite field or regular integers.
    
    Attributes:
        a (Union[int, FiniteFieldElement]): The a coefficient of the elliptic curve equation.
        b (Union[int, FiniteFieldElement]): The b coefficient of the elliptic curve equation.
    """

    def __init__(self, a: Union[int, FiniteFieldElement], b: Union[int, FiniteFieldElement]):
        """
        Initialize an EllipticCurve.

        Args:
            a (Union[int, FiniteFieldElement]): The 'a' coefficient.
            b (Union[int, FiniteFieldElement]): The 'b' coefficient.

        Raises:
            ValueError: If a or b is not of type int or FiniteFieldElement.
        """
        if not isinstance(a, (int, FiniteFieldElement)) or not isinstance(b, (int, FiniteFieldElement)):
            raise ValueError("Coefficients 'a' and 'b' must be either integers or FiniteFieldElement instances.")
        self.a = a
        self.b = b

    def __repr__(self) -> str:
        """Return a string representation of the elliptic curve."""
        if isinstance(self.a, FiniteFieldElement):
            return f'EllipticCurve(a={self.a.num}, b={self.b.num})'
        else:
            return f'EllipticCurve(a={self.a}, b={self.b})'
        

    def is_on_curve(self, x: Union[int, FiniteFieldElement], y: Union[int, FiniteFieldElement]) -> bool:
        """
        Check if the given point (x, y) lies on the elliptic curve.

        Args:
            x (Union[int, FiniteFieldElement]): The x-coordinate of the point.
            y (Union[int, FiniteFieldElement]): The y-coordinate of the point.

        Returns:
            bool: True if the point lies on the curve, False otherwise.
        """
        # IMPORTANT: the rounding here is done for the sole purpose of making the plotting
        # of point addition work visually and should never be used in practice!
        # In practice we want precision here so never use rounding when using this in an 
        # actual application!!!!
        return round(y**2, 5) == round(x**3 + self.a * x + self.b, 5)

    def plot_curve(self, x_min: float, x_max: float, linspace: int = 10000, point: Optional['ECPoint'] = False,
                   inverse: bool = False, lims: tuple = (5, 10), figsize: tuple = (10, 10), line: Optional[tuple] = False,
                   vline_x: Optional[float] = False, tangent_point: Optional['ECPoint'] = False) -> None:
        """
        Plot the elliptic curve between the given x-values.

        See visualisation.py for more detail.
        """
        _plot_curve(self, x_min, x_max, linspace=linspace, point=point, inverse=inverse, lims=lims, figsize=figsize, line=line,
                    vline_x=vline_x, tangent_point=tangent_point)

    def plot_point_addition(self, x_min: float, x_max: float, p1: 'ECPoint', p2: 'ECPoint', linspace: int = 10000,
                            lims: tuple = (5, 10), figsize: tuple = (10, 10), point_names: Optional[list[str]] = False) -> None:
        """
        Plot the addition of two points on the elliptic curve.

        See visualisation.py for more detail.
        """
        _plot_point_addition(self, x_min, x_max, p1, p2, linspace=linspace, lims=lims, figsize=figsize, point_names=point_names)


class ECPoint:
    """
    A class representing a point on an elliptic curve over a finite field or regular integers.

    Attributes:
        x (Union[int, FiniteFieldElement]): The x-coordinate of the point.
        y (Union[int, FiniteFieldElement]): The y-coordinate of the point.
        curve (EllipticCurve): The elliptic curve the point lies on.
    """

    def __init__(self, x: Optional[Union[int, FiniteFieldElement]], y: Optional[Union[int, FiniteFieldElement]], curve: EllipticCurve):
        """
        Initialize an ECPoint.

        Args:
            x (Optional[Union[int, FiniteFieldElement]]): The x-coordinate of the point. None if the point is at infinity.
            y (Optional[Union[int, FiniteFieldElement]]): The y-coordinate of the point. None if the point is at infinity.
            curve (EllipticCurve): The elliptic curve the point lies on.

        Raises:
            ValueError: If the point (x, y) does not lie on the elliptic curve.
        """
        self.x = x
        self.y = y
        self.curve = curve

        if self.x is not None and self.y is not None and not self.curve.is_on_curve(self.x, self.y):
            raise ValueError(f'({x}, {y}) is not on the curve {curve}')

    def __repr__(self) -> str:
        """Return a string representation of the point."""
        if self.x is None:
            return 'ECPoint(infinity)'
        elif isinstance(self.x, FiniteFieldElement):
            return f'ECPoint({self.x.num}, {self.y.num})_on_{self.curve}'
        else:
            return f'ECPoint({self.x}, {self.y})_on_{self.curve}'

    def __eq__(self, other: 'ECPoint') -> bool:
        """
        Check if this point is equal to another point.

        Args:
            other (ECPoint): The other point to compare with.

        Returns:
            bool: True if the points are equal (i.e., same x, y coordinates and the same curve), False otherwise.

        Raises:
            TypeError: If 'other' is not an instance of ECPoint.
        """
        if not isinstance(other, ECPoint):
            raise TypeError("Comparison must be with another ECPoint.")
        # IMPORTANT: This uses rounding solely or purposes of making the notebook work again...
        return (round(self.x, 5) == round(other.x, 5) and round(self.y, 5) == round(other.y, 5) 
                and self.curve == other.curve)

    def __ne__(self, other: 'ECPoint') -> bool:
        """
        Check if this point is not equal to another point.

        Args:
            other (ECPoint): The other point to compare with.

        Returns:
            bool: True if the points are not equal (i.e., different x, y coordinates or different curves), False otherwise.

        Raises:
            TypeError: If 'other' is not an instance of ECPoint.
        """
        return not self.__eq__(other)

    def __add__(self, other: 'ECPoint') -> 'ECPoint':
        """
        Add this point to another point.

        Args:
            other (ECPoint): The other point to add.

        Returns:
            ECPoint: The resulting point after addition.

        Raises:
            TypeError: If the points are not on the same curve.
        """
        if self.curve != other.curve:
            raise TypeError('Points are not on the same curve')

        # Handle the case of the point at infinity
        if self.x is None:
            return other
        if other.x is None:
            return self

        # Handle the case of P + -P = 0 (the point at infinity)
        if self.x == other.x and self.y != other.y:
            return self.__class__(None, None, self.curve)

        # Handle the case of P + P = 2P (point doubling)
        if self == other:
            if self.y == 0 * self.x:
                return self.__class__(None, None, self.curve)
            else:
                s = (3 * self.x**2 + self.curve.a) / (2 * self.y)
                new_x = s**2 - 2 * self.x
                new_y = s * (self.x - new_x) - self.y
                return self.__class__(new_x, new_y, self.curve)

        # Handle the case of P + Q (point addition)
        if self.x != other.x:
            s = (other.y - self.y) / (other.x - self.x)
            new_x = s**2 - self.x - other.x
            new_y = s * (self.x - new_x) - self.y
            return self.__class__(new_x, new_y, self.curve)

        raise ValueError('Addition case not handled')

    def __rmul__(self, coefficient: int) -> 'ECPoint':
        """
        Perform scalar multiplication using the double-and-add algorithm.

        Args:
            coefficient (int): The integer coefficient for scalar multiplication.

        Returns:
            ECPoint: The resulting point after performing scalar multiplication.

        Raises:
            TypeError: If the coefficient is not an integer.
        """
        if not isinstance(coefficient, int):
            raise TypeError("The coefficient must be an integer.")
        
        coef = coefficient
        current = self
        result = self.__class__(None, None, self.curve)
        while coef:
            if coef & 1:
                result += current
            current += current
            coef >>= 1
        return result