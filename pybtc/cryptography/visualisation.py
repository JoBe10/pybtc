from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from typing import Optional, Union, TYPE_CHECKING
from finite_field import FiniteFieldElement
import math

def _plot_curve(curve, x_min: float, x_max: float, linspace: int = 10000, point = False,
                inverse: bool = False, lims: tuple = (5, 10), figsize: tuple = (10, 10), line: Optional[tuple] = False,
                vline_x: Optional[float] = False, tangent_point = False) -> None:
    """
    Plot the elliptic curve between the given x-values.

    Args:
        curve (EllipticCurve): The elliptic curve to plot.
        x_min (float): The minimum x-value.
        x_max (float): The maximum x-value.
        linspace (int): The number of points to generate between x_min and x_max. Default is 10000.
        point (Optional[ECPoint]): A point to plot on the curve. Default is False.
        inverse (bool): Whether to plot the inverse of the point. Default is False.
        lims (tuple): Limits for the plot axes as (x_limit, y_limit). Default is (5, 10).
        figsize (tuple): Size of the figure as (width, height). Default is (10, 10).
        line (Optional[tuple]): Coefficients (m, c) of a straight line y = mx + c to plot. Default is False.
        vline_x (Optional[float]): x-coordinate of a vertical line to plot. Default is False.
        tangent_point (Optional[ECPoint]): Point at which to plot the tangent line. Default is False.

    Raises:
        ValueError: If any input is invalid.
    """
    # Validate input arguments
    if not isinstance(lims, tuple) or len(lims) != 2:
        raise ValueError('lims must be a tuple of length 2')
    if not isinstance(figsize, tuple) or len(figsize) != 2:
        raise ValueError('figsize must be a tuple of length 2')
    if not isinstance(linspace, int) or linspace <= 0:
        raise ValueError('linspace must be a positive integer')
    if line and (not isinstance(line, tuple) or len(line) != 2):
        raise ValueError('line must be a tuple of length 2')
    if vline_x and not isinstance(vline_x, (int, float)):
        raise ValueError('vline_x must be a number')

    # Generate x values and corresponding y values for the elliptic curve
    x = np.linspace(x_min, x_max, linspace)
    y_elliptic_positive = [np.sqrt(xi**3 + curve.a * xi + curve.b) if xi**3 + curve.a * xi + curve.b >= 0 else np.nan for xi in x]
    y_elliptic_negative = [-yi for yi in y_elliptic_positive]

    # Set up the plot
    plt.figure(figsize=figsize)

    sign_a = f'{curve.a:+}x' if curve.a != 0 else ''
    sign_b = f'{curve.b:+}' if curve.b != 0 else ''
    label = f'Elliptic Curve: $y^2 = x^3 {sign_a} {sign_b}$'

    plt.plot(x, y_elliptic_positive, label=label, color='red')
    plt.plot(x, y_elliptic_negative, color='red')

    if point:
        plt.plot(point.x, point.y, 'go')
        if inverse:
            plt.plot(point.x, -point.y, 'go')
            plt.text(point.x-0.1, point.y-0.1, '$P$', fontsize=12, ha='right')
            plt.text(point.x-0.1, -point.y-0.2, '$P^{-1}$', fontsize=12, ha='right')

    if line:
        y_line = [line[0] * xi + line[1] for xi in x]
        sign_m = f'{line[0]}x' if line[0] != 0 else ''
        sign_b_line = (f'{line[1]:+}' if line[1] != 0 else '') if line[0] != 0 else (f'{line[1]}' if line[1] != 0 else '')
        label_line = f'Straight Line: $y = {sign_m} {sign_b_line}$'
        plt.plot(x, y_line, color='blue', label=label_line)

    if vline_x:
        label_vline = f'Vertical Line: $x = {vline_x}$'
        plt.axvline(vline_x, color='blue', label=label_vline)

    if tangent_point:
        # Plot the tangent line at the point
        s = (3 * tangent_point.x**2 + curve.a) / (2 * tangent_point.y)
        y_tangent = s * (x - tangent_point.x) + tangent_point.y
        intercept = s * (0 - tangent_point.x) + tangent_point.y
        sign_m_tangent = f'{s}x' if s != 0 else ''
        sign_b_tangent = f'{intercept:+}' if intercept != 0 else ''
        label_tline = f'Tangent Line: $y = {sign_m_tangent} {sign_b_tangent}$'
        plt.plot(x, y_tangent, 'blue', label=label_tline)

    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.legend()

    # Set axis limits
    plt.xlim(-lims[0], lims[0])
    plt.ylim(-lims[1], lims[1])

    # Set the aspect ratio
    plt.gca().set_aspect(lims[0] / lims[1], adjustable='datalim')

    plt.show()


def _plot_point_addition(curve, x_min: float, x_max: float, p1, p2, linspace: int = 10000,
                        lims: tuple = (5, 10), figsize: tuple = (10, 10), point_names: Optional[list[str]] = False) -> None:
    """
    Plot the addition of two points on the elliptic curve.

    Args:
        curve (EllipticCurve): The elliptic curve to plot.
        x_min (float): The minimum x-value.
        x_max (float): The maximum x-value.
        p1 (ECPoint): The first point.
        p2 (ECPoint): The second point.
        linspace (int): The number of points to generate between x_min and x_max. Default is 10000.
        lims (tuple): Limits for the plot axes as (x_limit, y_limit). Default is (5, 10).
        figsize (tuple): Size of the figure as (width, height). Default is (10, 10).
        point_names (Optional[list[str]]): Custom names for the points. Default is False.

    Raises:
        ValueError: If the points are not on the curve or any input is invalid.
    """
    if p1.curve != curve or p2.curve != curve:
        raise ValueError('Both points must be on the same curve')

    # Compute the resulting point
    result = p1 + p2

    # Generate x values and corresponding y values for the elliptic curve
    x = np.linspace(x_min, x_max, linspace)
    y_elliptic_positive = [np.sqrt(xi**3 + curve.a * xi + curve.b) if xi**3 + curve.a * xi + curve.b >= 0 else np.nan for xi in x]
    y_elliptic_negative = [-yi for yi in y_elliptic_positive]

    # Set up the plot
    plt.figure(figsize=figsize)

    sign_a = f'{curve.a:+}x' if curve.a != 0 else ''
    sign_b = f'{curve.b:+}' if curve.b != 0 else ''
    label = f'Elliptic Curve: $y^2 = x^3 {sign_a} {sign_b}$'

    plt.plot(x, y_elliptic_positive, label=label, color='red')
    plt.plot(x, y_elliptic_negative, color='red')

    # Set label position adjustment
    text_adjustment_x = -0.1
    text_adjustment_y = 0

    if p1 == p2:
        # Plot the point and its label
        plt.plot(p1.x, p1.y, 'go')
        plt.text(p1.x + text_adjustment_x, p1.y + text_adjustment_y, 'A', fontsize=12, ha='right')

        # Plot the tangent line at the point
        s = (3 * p1.x**2 + curve.a) / (2 * p1.y)
        y_tangent = s * (x - p1.x) + p1.y
        plt.plot(x, y_tangent, 'blue', linestyle='--')

        # Plot the resulting point and label it as 2A
        plt.plot(result.x, result.y, 'go')
        plt.text(result.x + text_adjustment_x, result.y + text_adjustment_y, '2A', fontsize=12, ha='right')

        # Flip over the x-axis to get the correct point as a result of addition
        plt.plot(result.x, -result.y, 'go')
        plt.axvline(result.x, color='blue', linestyle='--')

    elif p1.x == p2.x:
        # Plot the points and their labels
        plt.plot(p1.x, p1.y, 'go')
        plt.text(p1.x + text_adjustment_x, p1.y + text_adjustment_y, 'A', fontsize=12, ha='right')
        plt.plot(p2.x, p2.y, 'go')
        plt.text(p2.x + text_adjustment_x, p2.y + text_adjustment_y, 'B', fontsize=12, ha='right')

        plt.axvline(p1.x, color='blue', linestyle='--')
    
    else:
        if point_names:
            p1_name, p2_name, p3_name, p4_name = point_names[0], point_names[1], point_names[2], point_names[3]
        else:
            p1_name, p2_name, p3_name, p4_name = 'A', 'B', 'C', 'A+B'

        # Plot the points and their labels
        plt.plot(p1.x, p1.y, 'go')
        plt.text(p1.x + text_adjustment_x, p1.y + text_adjustment_y, p1_name, fontsize=12, ha='right')
        plt.plot(p2.x, p2.y, 'go')
        plt.text(p2.x + text_adjustment_x, p2.y + text_adjustment_y, p2_name, fontsize=12, ha='right')

        # Plot the line connecting the points and the vertical line for the result
        slope = (p2.y - p1.y) / (p2.x - p1.x)
        y_line = slope * (x - p1.x) + p1.y
        plt.plot(x, y_line, 'blue', linestyle='--')

        # Plot the resulting point and label it as A+B
        plt.plot(result.x, result.y, 'go')
        plt.text(result.x + text_adjustment_x, result.y + text_adjustment_y, p4_name, fontsize=12, ha='right')

        # Flip over the x-axis to get the correct point as a result of addition
        plt.plot(result.x, -result.y, 'go')
        plt.text(result.x + text_adjustment_x, -result.y + text_adjustment_y, p3_name, fontsize=12, ha='right')
        plt.axvline(result.x, color='blue', linestyle='--')

    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.legend()

    # Set axis limits
    plt.xlim(-lims[0], lims[0])
    plt.ylim(-lims[1], lims[1])

    # Set the aspect ratio
    plt.gca().set_aspect(lims[0] / lims[1], adjustable='datalim')

    plt.show()


def _plot_curve_over_field(curve, figsize: tuple = (10, 6), point = False, p_label: str = 'G') -> None:
    """
    Plot the points on the elliptic curve over a finite field.

    Args:
        curve (EllipticCurve): The elliptic curve to plot.
        figsize (tuple): Size of the figure as (width, height). Default is (10, 6).
        point (Optional[ECPoint]): A point to plot on the curve. Default is False.
        p_label (str): Label for the point. Default is 'G'.

    Raises:
        ValueError: If the elliptic curve is not defined over a finite field or any input is invalid.
    """
    if not isinstance(curve.a, FiniteFieldElement) or not isinstance(curve.b, FiniteFieldElement):
        raise ValueError('Elliptic curve coefficients must be FiniteFieldElement instances for this function')

    field_order = curve.a.prime

    if field_order > 547:
        raise ValueError('Field order must be at most 547 in order to reduce plotting time')

    if not isinstance(figsize, tuple) or len(figsize) != 2:
        raise ValueError('figsize must be a tuple of length 2')

    # Generate all points on the elliptic curve over the finite field
    points = []
    for x in range(field_order):
        for y in range(field_order):
            x_field = FiniteFieldElement(x, field_order)
            y_field = FiniteFieldElement(y, field_order)
            if curve.is_on_curve(x_field, y_field):
                points.append((x, y))
    
    # Elliptic curve equation label
    sign_a = f'{curve.a.num:+}x' if curve.a.num != 0 else ''
    sign_b = f'{curve.b.num:+}' if curve.b.num != 0 else ''

    # Plotting the points
    plt.figure(figsize=figsize)
    x_coords, y_coords = zip(*points)
    plt.scatter(x_coords, y_coords, c='blue', marker='.')
    if point:
        plt.scatter([point.x.num], [point.y.num], c='red')
        plt.text(point.x.num-0.4, point.y.num-0.4, p_label, fontsize=12, ha='right')
    #plt.title(f'Points on elliptic curve $y^2 = x^3 {sign_a} {sign_b} over finite field of order {field_order}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()


def _plot_point_addition_over_field(curve, p1, p2, figsize: tuple = (10, 6), max_wraps: int = 20) -> None:
    """
    Plot point addition on the elliptic curve over a finite field.

    Args:
        curve (EllipticCurve): The elliptic curve to plot.
        p1 (ECPoint): The first point.
        p2 (ECPoint): The second point.
        figsize (tuple): Size of the figure as (width, height). Default is (10, 6).
        max_wraps (int): Maximum number of times to wrap around the field order to find the result. Default is 20.

    Raises:
        ValueError: If the elliptic curve is not defined over a finite field or any input is invalid.
    """
    if not isinstance(curve.a, FiniteFieldElement) or not isinstance(curve.b, FiniteFieldElement):
        raise ValueError('Elliptic curve coefficients must be FiniteFieldElement instances for this function')
    if p1.curve != curve or p2.curve != curve:
        raise ValueError('Both points must be on the same curve')
    if not isinstance(figsize, tuple) or len(figsize) != 2:
        raise ValueError('figsize must be a tuple of length 2')

    field_order = curve.a.prime

    if field_order > 547:
        raise ValueError('Field order must be at most 547 in order to reduce plotting time')

    # Generate all points on the elliptic curve over the finite field
    points = []
    for x in range(field_order):
        for y in range(field_order):
            x_field = FiniteFieldElement(x, field_order)
            y_field = FiniteFieldElement(y, field_order)
            if curve.is_on_curve(x_field, y_field):
                points.append((x, y))
    
    # Elliptic curve equation label
    sign_a = f'{curve.a.num:+}x' if curve.a.num != 0 else ''
    sign_b = f'{curve.b.num:+}' if curve.b.num != 0 else ''

    # Plotting the points
    plt.figure(figsize=figsize)
    x_coords, y_coords = zip(*points)
    plt.scatter(x_coords, y_coords, c='blue', marker='.')
    #plt.title(f'Point addition on elliptic curve $y^2 = x^3 {sign_a} {sign_b}$ over finite field of order {field_order}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)

    # Plot points A and B
    plt.scatter([p1.x.num, p2.x.num], [p1.y.num, p2.y.num], c='red')
    if p1 == p2:
        plt.text(p1.x.num-0.1, p1.y.num, 'P', fontsize=12, ha='right')
    else:
        plt.text(p1.x.num-0.1, p1.y.num, 'A', fontsize=12, ha='right')
    if p1 != p2:
        plt.text(p2.x.num-0.1, p2.y.num, 'B', fontsize=12, ha='right')

    # Compute the resulting point
    result = p1 + p2
    result_inverse_y = field_order - result.y.num

    # Determine the slope and intercept of the line
    if p1 == p2:
        s = (3 * p1.x.num ** 2 + curve.a.num) / (2 * p1.y.num)
    else:
        s = (p2.y.num - p1.y.num) / (p2.x.num - p1.x.num)

    wraparound_y = p1.y.num
    wraparound_x = p1.x.num

    point_reached = False
    i = 1
    while not point_reached and i < max_wraps + 1:
        y_int = wraparound_y - wraparound_x * s
        if round(result.x.num * s + y_int) == result_inverse_y or i == max_wraps:
            plt.plot([wraparound_x, result.x.num], [wraparound_y, result_inverse_y], 'orange')
            point_reached = True
            break

        # Plot the line connecting A and B
        x_vals = []
        y_vals = []

        x_vals.append(wraparound_x)
        y_vals.append(wraparound_y)
        for x_i in range(math.ceil(wraparound_x), field_order):
            y_i = s * x_i + y_int
            if y_i > (field_order - 1):
                last_x = (field_order - 1 - y_int) / s
                x_vals.append(last_x)
                y_vals.append(field_order - 1)
                wraparound_x = last_x
                wraparound_y = 0
                break
            if x_i == (field_order - 1):
                last_y = (field_order - 1) * s + y_int
                x_vals.append(field_order - 1)
                y_vals.append(last_y)
                wraparound_x = 0
                wraparound_y = last_y
                break
            else:
                x_vals.append(x_i)
                y_vals.append(y_i)

        plt.plot(x_vals, y_vals, 'orange')
        i += 1

    # Plot the resulting point and label it
    plt.scatter([result.x.num], [result.y.num], c='green')
    if p1 == p2:
        plt.text(result.x.num-0.1, result.y.num, '2P', fontsize=12, ha='right')
    else:
        plt.text(result.x.num-0.1, result.y.num, 'A+B', fontsize=12, ha='right')

    # Plot the flip over the y = p/2 axis
    plt.scatter([result.x.num], [result_inverse_y], c='green')
    if p1 != p2:
        plt.text(result.x.num-0.1, result_inverse_y, 'C', fontsize=12, ha='right')

    # Plot the vertical dotted line from C to A+B
    plt.plot([result.x.num, result.x.num], [result_inverse_y, result.y.num], 'orange', linestyle='dotted')

    plt.show()


def _plot_cyclic_subgroup(curve, generator, figsize: tuple = (10, 6), max_loops: int = 1) -> None:
    """
    Plot the cyclic subgroup generated by repeatedly adding the generator point to itself.

    Args:
        curve (EllipticCurve): The elliptic curve to plot.
        generator (ECPoint): The generator point to use for the cyclic group.
        figsize (tuple): Size of the figure as (width, height). Default is (10, 6).
        max_loops (int): Maximum number of loops to complete. Default is 1.

    Raises:
        ValueError: If the elliptic curve is not defined over a finite field or any input is invalid.
    """
    from ec import ECPoint
    if not isinstance(curve.a, FiniteFieldElement) or not isinstance(curve.b, FiniteFieldElement):
        raise ValueError('Elliptic curve coefficients must be FiniteFieldElement instances for this function')

    field_order = curve.a.prime

    if field_order > 547:
        raise ValueError('Field order must be at most 547 in order to reduce plotting time')

    if not isinstance(figsize, tuple) or len(figsize) != 2:
        raise ValueError('figsize must be a tuple of length 2')

    points = [(generator.x.num, generator.y.num)]
    nums = [str(1)]
    result = generator
    loops = 0

    for i in range(2, field_order * 2 * max_loops):
        result += generator
        if result == ECPoint(None, None, curve):  # Handle point at infinity
            loops += 1
            if loops >= max_loops:
                break
            else:
                continue
        if loops >= 1:
            nums[i % (len(points)) - 2] += ', ' + str(i)
        else:
            points.append((result.x.num, result.y.num))
            nums.append(str(i))

    plt.figure(figsize=figsize)
    x_coords, y_coords = zip(*points)
    plt.scatter(x_coords, y_coords, c='blue', marker='.')

    for i in range(len(nums)):
        plt.text(points[i][0]+0.1, points[i][1]-0.1, nums[i], c='red')

    sign_a = f'{curve.a.num:+}x' if curve.a.num != 0 else ''
    sign_b = f'{curve.b.num:+}' if curve.b.num != 0 else ''

    #plt.title(f'Cyclic Subgroup on Elliptic Curve $y^2 = x^3 {sign_a} {sign_b}$ over Finite Field of order {field_order}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()