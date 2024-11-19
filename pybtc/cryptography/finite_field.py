from typing import Union

class FiniteFieldElement:
    """
    A class representing an element of a finite field.

    Attributes:
        num (int): The value of the element.
        prime (int): The prime order of the finite field.
    """

    def __init__(self, num: int, prime: int):
        """
        Initialize a FiniteFieldElement.

        Args:
            num (int): The value of the element in the field.
            prime (int): The prime number defining the field.

        Raises:
            ValueError: If num is not in the range [0, prime) or if prime is not positive.
        """
        if not isinstance(prime, int) or prime <= 1:
            raise ValueError("Prime must be a positive integer greater than 1.")
        if not isinstance(num, int) or num < 0 or num >= prime:
            raise ValueError(f"Num {num} not in field range 0 to {prime - 1}")
        self.num = num
        self.prime = prime

    def __repr__(self):
        """Return a string representation of the element."""
        return f'FiniteFieldElement_{self.prime}({self.num})'

    def __eq__(self, other: Union['FiniteFieldElement', None]) -> bool:
        """
        Check if this element is equal to another.

        Args:
            other (Union[FiniteFieldElement, None]): Another element to compare with.

        Returns:
            bool: True if the elements are equal, False otherwise.
        """
        if other is None:
            return False
        return self.num == other.num and self.prime == other.prime

    def __ne__(self, other: Union['FiniteFieldElement', None]) -> bool:
        """
        Check if this element is not equal to another.

        Args:
            other (Union[FiniteFieldElement, None]): Another element to compare with.

        Returns:
            bool: True if the elements are not equal, False otherwise.
        """
        return not (self == other)

    def __add__(self, other: Union['FiniteFieldElement', int]) -> 'FiniteFieldElement':
        """
        Add this element to another.

        Args:
            other (Union[FiniteFieldElement, int]): Another element or an integer to add.

        Returns:
            FiniteFieldElement: The result of the addition.

        Raises:
            TypeError: If the fields are not the same or if other is not an integer or FiniteFieldElement.
        """
        if isinstance(other, int):
            other = self.__class__(other, self.prime)
        elif not isinstance(other, FiniteFieldElement):
            raise TypeError("Operand must be an integer or FiniteFieldElement")
        
        if self.prime != other.prime:
            raise TypeError('Cannot add two numbers in different fields')
        
        num = (self.num + other.num) % self.prime
        return self.__class__(num, self.prime)

    def __sub__(self, other: Union['FiniteFieldElement', int]) -> 'FiniteFieldElement':
        """
        Subtract another element from this one.

        Args:
            other (Union[FiniteFieldElement, int]): Another element or an integer to subtract.

        Returns:
            FiniteFieldElement: The result of the subtraction.

        Raises:
            TypeError: If the fields are not the same or if other is not an integer or FiniteFieldElement.
        """
        if isinstance(other, int):
            other = self.__class__(other, self.prime)
        elif not isinstance(other, FiniteFieldElement):
            raise TypeError("Operand must be an integer or FiniteFieldElement")
        
        if self.prime != other.prime:
            raise TypeError('Cannot subtract two numbers in different fields')
        
        num = (self.num - other.num) % self.prime
        return self.__class__(num, self.prime)

    def __mul__(self, other: Union['FiniteFieldElement', int]) -> 'FiniteFieldElement':
        """
        Multiply this element by another.

        Args:
            other (Union[FiniteFieldElement, int]): Another element or an integer to multiply by.

        Returns:
            FiniteFieldElement: The result of the multiplication.

        Raises:
            TypeError: If the fields are not the same or if other is not an integer or FiniteFieldElement.
        """
        if isinstance(other, int):
            other = self.__class__(other, self.prime)
        elif not isinstance(other, FiniteFieldElement):
            raise TypeError("Operand must be an integer or FiniteFieldElement")
        
        if self.prime != other.prime:
            raise TypeError('Cannot multiply two numbers in different fields')
        
        num = (self.num * other.num) % self.prime
        return self.__class__(num, self.prime)

    def __pow__(self, exponent: int) -> 'FiniteFieldElement':
        """
        Raise this element to a power.

        Args:
            exponent (int): The power to raise the element to.

        Returns:
            FiniteFieldElement: The result of the exponentiation.
        """        
        # Fermat's Little Theorem: a^(p-1) ≡ 1 (mod p)
        # To find a^b mod p, we can reduce b mod (p-1)
        n = exponent % (self.prime - 1)
        num = pow(self.num, n, self.prime)
        return self.__class__(num, self.prime)

    def __truediv__(self, other: Union['FiniteFieldElement', int]) -> 'FiniteFieldElement':
        """
        Divide this element by another.

        Args:
            other (Union[FiniteFieldElement, int]): Another element or an integer to divide by.

        Returns:
            FiniteFieldElement: The result of the division.

        Raises:
            TypeError: If the fields are not the same or if other is not an integer or FiniteFieldElement.
            ZeroDivisionError: If division by zero is attempted.
        """
        if isinstance(other, int):
            other = self.__class__(other, self.prime)
        elif not isinstance(other, FiniteFieldElement):
            raise TypeError("Operand must be an integer or FiniteFieldElement")
        
        if self.prime != other.prime:
            raise TypeError('Cannot divide two numbers in different fields')
        if other.num == 0:
            raise ZeroDivisionError("Cannot divide by zero in finite fields.")
        
        # Division in a finite field is multiplication by the modular inverse
        # Fermat's Little Theorem: a^(p-2) ≡ a^(-1) (mod p)
        num = (self.num * pow(other.num, self.prime - 2, self.prime)) % self.prime
        return self.__class__(num, self.prime)

    def __rmul__(self, coefficient: Union[int, 'FiniteFieldElement']) -> 'FiniteFieldElement':
        """
        Right multiplication of this element by a scalar or another FiniteFieldElement.

        Args:
            coefficient (Union[int, FiniteFieldElement]): The scalar or element to multiply by.

        Returns:
            FiniteFieldElement: The result of the multiplication.

        Raises:
            TypeError: If coefficient is neither an integer nor a FiniteFieldElement.
        """
        if isinstance(coefficient, int):
            num = (self.num * coefficient) % self.prime
        elif isinstance(coefficient, FiniteFieldElement):
            if self.prime != coefficient.prime:
                raise TypeError('Cannot multiply two numbers in different fields')
            num = (self.num * coefficient.num) % self.prime
        else:
            raise TypeError("Operand must be an integer or FiniteFieldElement")
        
        return self.__class__(num, self.prime)

    def mul_inverse(self) -> 'FiniteFieldElement':
        """
        Compute the multiplicative inverse of this element.

        Returns:
            FiniteFieldElement: The multiplicative inverse.

        Raises:
            ZeroDivisionError: If the element is zero (no inverse).
        """
        if self.num == 0:
            raise ZeroDivisionError("The multiplicative inverse of zero does not exist.")
        
        num = pow(self.num, self.prime - 2, self.prime)
        return self.__class__(num, self.prime)
    
    def add_inverse(self) -> 'FiniteFieldElement':
        """
        Compute the additive inverse of this element.

        Returns:
            FiniteFieldElement: The additive inverse of this element.

        Note:
            In a finite field, the additive inverse of a is -a mod p, where p is the prime modulus.
        """
        num = (-self.num) % self.prime
        return self.__class__(num, self.prime)