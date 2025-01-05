import numpy as np
from scipy.optimize import minimize

def calculate_roll_center(upper, lower, cop):
    """
    Calculate the roll center (RC) based on 3D suspension geometry and tire CoP.
    Assumes suspension symmetry about the vehicle centerline.

    Parameters:
        upper (np.ndarray): 3x3 matrix of upper control arm pickup points.
                        [[outer], [inner front], [inner rear]] points.
        lower (np.ndarray): 3x3 matrix of lower control arm pickup points.
                    [[outer], [inner front], [inner rear]] points.
        cop (np.ndarray): 1x3 array for tire center of pressure (CoP).

    Returns:
        np.ndarray: Roll center coordinates [x, y, z] in the global coordinate system.
        dict: Additional geometry data including instant centers and angles.
    """
    # Input validation
    if not all(isinstance(x, np.ndarray) for x in [upper, lower, cop]):
        raise ValueError("Inputs must be numpy arrays")
    if upper.shape != (3, 3) or lower.shape != (3, 3):
        raise ValueError("Control arm matrices must be 3x3")
    if cop.shape != (3,):
        raise ValueError("CoP must be 1x3 array")

    # Calculate control arm planes normal vectors
    def get_plane_normal(points):
        """Calculate normal vector of plane defined by three points."""
        v1 = points[1] - points[0]
        v2 = points[2] - points[0]
        return np.cross(v1, v2)

    upper_normal = get_plane_normal(upper)
    lower_normal = get_plane_normal(lower)

    # Project arms onto front view (XZ plane)
    upper_outer_2d = upper[0, [0, 2]]  # [x, z]
    upper_inner_2d = upper[1, [0, 2]]  # Using front inner point
    lower_outer_2d = lower[0, [0, 2]]
    lower_inner_2d = lower[1, [0, 2]]
    tire_cop_2d = cop[[0, 2]]

    # Calculate control arm angles in front view
    def get_arm_angle(outer, inner):
        """Calculate arm angle from horizontal in front view."""
        delta = inner - outer
        return np.degrees(np.arctan2(delta[1], delta[0]))

    upper_angle = get_arm_angle(upper_outer_2d, upper_inner_2d)
    lower_angle = get_arm_angle(lower_outer_2d, lower_inner_2d)

    # Find instant center (IC) in front view
    def line_equation(p1, p2):
        """Return slope and intercept of line through two points."""
        m = (p2[1] - p1[1]) / (p2[0] - p1[0])
        c = p1[1] - m * p1[0]
        return m, c

    m_upper, c_upper = line_equation(upper_outer_2d, upper_inner_2d)
    m_lower, c_lower = line_equation(lower_outer_2d, lower_inner_2d)
    
    ic_x = (c_lower - c_upper) / (m_upper - m_lower)
    ic_z = m_upper * ic_x + c_upper
    instant_center_2d = np.array([ic_x, ic_z])

    # Calculate roll center height
    # RC is where line from IC to CoP intersects vehicle centerline (x=0)
    m_ic_to_cop = (tire_cop_2d[1] - ic_z) / (tire_cop_2d[0] - ic_x)
    c_ic_to_cop = ic_z - m_ic_to_cop * ic_x
    rc_z = c_ic_to_cop  # At x=0

    # Consider 3D effects
    # Calculate swing arm length (SAL) considering y-axis position
    sal = np.sqrt(ic_x**2 + (upper[1, 1] - upper[0, 1])**2)  # Using y-distance
    
    # Calculate roll center including y-coordinate
    # Using average y-position of inner pickup points
    rc_y = (np.mean(upper[:, 1]) + np.mean(lower[:, 1])) / 2
    
    roll_center = np.array([0, rc_y, rc_z])

    # Package additional geometry data
    geometry_data = {
        'instant_center_2d': instant_center_2d,
        'swing_arm_length': sal,
        'upper_arm_angle': upper_angle,
        'lower_arm_angle': lower_angle,
        'upper_normal': upper_normal,
        'lower_normal': lower_normal
    }

    return roll_center, geometry_data

def rotate_points(points, axis_point, rotation_axis, angle, height):
    """
    Rotate points around a given axis in 3D space.

    Parameters:
        points (np.array): Nx3 array of points to be rotated.
        axis_point (np.array): A point on the rotation axis (3D).
        rotation_axis (np.array): Direction vector of the rotation axis (3D).
        angle (float): Rotation angle in radians.

    Returns:
        np.array: Rotated points.
    """
    
    points = jounce_innerpoints(points, height)
    
    # Normalize the rotation axis
    rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
    
    # Translate points to the rotation axis origin
    translated_points = points[1:] - axis_point
    
    # Compute the rotation matrix using the Rodrigues' rotation formula
    K = np.array([[0, -rotation_axis[2], rotation_axis[1]],
                  [rotation_axis[2], 0, -rotation_axis[0]],
                  [-rotation_axis[1], rotation_axis[0], 0]])
    
    I = np.eye(3)
    R = I + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)
    
    # Apply rotation and translate back
    rotated_points = (R @ translated_points.T).T + axis_point
    
    # Combine the unrotated first point with the rotated points
    modified_points = np.vstack((points[0], rotated_points))
    
    return modified_points

def jounce_innerpoints(points, height):

    """
    Move the inner points of the suspension system by teh prescribed distance
    
    Parameters:
        points (np.array): 3x3 array of suspension joint locations
        height (float): movement of the sprung mass (m)
    
    Returns:
        np.array: moved inner points
    
    """
    
    inner_points_z = points[1:, 2:]
    inner_points_z = inner_points_z + height
    
    points[1:,2:] = inner_points_z
    
    return points

def solve_outer_position(points, target_lengths):
    """
    Solve for the outer point position while preserving distances to the inner points.

    Parameters:
        inner_points (np.array): Inner points (2x3), extracted from 'points' matrix input (4x3)
        initial_outer (np.array): Initial guess for the outer point (1x3), extracted from 'points' matrix input (4x3)
        target_lengths (np.array): Target lengths to maintain between the outer and inner points (2x1)

    Returns:
        np.array: The updated outer point position.
    """
    
    inner_points = points[1:3]
    initial_outer = points[0]
    constraint_point = points[3]
    link_lengths = target_lengths[:2]
    third_length = target_lengths[2]
    
    def objective(outer_point):
        # Distances to inner points (leading and trailing arms)
        distances = np.linalg.norm(inner_points - outer_point, axis=1)
        inner_error = np.sum((distances - link_lengths) ** 2)

        # Constraint for the fixed distance to the reference point
        constraint_dist = np.linalg.norm(constraint_point - outer_point)
        constraint_error = (constraint_dist - third_length) ** 2

        return inner_error + constraint_error

    # Solve using the optimization library
    result = minimize(objective, initial_outer, method='SLSQP')
    if result.success:
        points[0] = result.x  # Update the first row of the points matrix
        return result.x
    else:
        raise ValueError("Optimization failed:", result.message)