import numpy as np
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation as R

def calculate_rc_simple(upper, lower, cop):
    """
    Calculates the roll centre based on 3D geometry:
    1. Project a-arms to YZ plane (X=0) 
    2. Find equation of upper & lower lines
    3. Find SV IC
    4. Equation of line from tire_cop==>SV IC
    5. Find z when y=0
    6. Return rc_z
    """
    
    # Find equation of line in form y = m*x + c
    def line_equation(p1, p2):
        """Return slope and intercept of line through two points."""
        m = (p2[1] - p1[1]) / (p2[0] - p1[0])       # Find slope of line
        c = p1[1] - m * p1[0]                       # Find intercept of line
        return m, c
    
    # Step 1: Project arms onto front view (XZ plane)
    upper_outer_2d = upper[0, [1,2]]  # [y, z]
    upper_inner_2d = upper[1, [1,2]]  # Using front inner point
    lower_outer_2d = lower[0, [1,2]]
    lower_inner_2d = lower[1, [1,2]]
    tire_cop_2d = cop[[1, 2]]
    
    # Step 2: Find equation of line in XY plane
    m_upper, c_upper = line_equation(upper_outer_2d, upper_inner_2d)
    m_lower, c_lower = line_equation(lower_outer_2d, lower_inner_2d)
    
    # Step 3: Find SV IC
    ic_y = (c_lower - c_upper) / (m_upper - m_lower)
    ic_z = m_upper * ic_y + c_upper
    SV_IC = [ic_y, ic_z]
    
    # Step 4: Find line of action of force
    m_load, rc_z = line_equation(tire_cop_2d, SV_IC)
    output_array = np.array([0,0,rc_z,ic_y,ic_z])

    return output_array

def rotate_points(points, axis_point, rotation_axis, angle):
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
# Smoothly adjust heights before rotation

    # Translate points to the axis origin
    translated_points = points - axis_point

    # Create the rotation object
    rotation = R.from_rotvec(angle * rotation_axis)

    # Apply the rotation
    rotated_translated_points = rotation.apply(translated_points)

    # Translate back to the original position
    rotated_points = rotated_translated_points + axis_point

    return rotated_points

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

def trail_calc(upper_points, lower_points, cop):
    
    # Import the outer joint points of the suspension system
    
    # Lower outer point
    p0 = lower_points[0,:]
    # Upper outer point
    p1 = upper_points[0,:]
    
    
    # Calculate the vector
    tz = -p0[2]/(p1[2]-p0[2])
    x0 = p0[0] + tz*(p1[0]-p0[0])
    y0 = p0[1] + tz*(p1[1]-p0[1])
    
    trail = (x0-cop[0])
    scrub = (y0-cop[1])
    
    
    trails = np.array([[trail, scrub]])
    
    return trails

def suspension_geometry_calc(upper_points, lower_points, cop):
    """
    Calculate trail, scrub, kingpin inclination, and camber angle.
    
    Args:
        upper_points: np.array, upper suspension points [outer, inner].
        lower_points: np.array, lower suspension points [outer, inner].
        cop: np.array, center of the contact patch (x, y, z).
    
    Returns:
        dict: Dictionary containing trail, scrub, kingpin inclination (degrees), and camber angle (degrees).
    """
    # Lower outer point
    p0 = lower_points[0, :]
    # Upper outer point
    p1 = upper_points[0, :]
    
    # Calculate the vector for kingpin axis
    kingpin_vector = p1 - p0
    
    # Trail and Scrub
    tz = -p0[2] / (p1[2] - p0[2])  # Intersection with the ground (z = 0)
    x0 = p0[0] + tz * (p1[0] - p0[0])
    y0 = p0[1] + tz * (p1[1] - p0[1])
    trail = x0 - cop[0]
    scrub = y0 - cop[1]
    
    # Kingpin Inclination (XZ-plane)
    kpi_vector_xz = np.array([kingpin_vector[0], kingpin_vector[2]])  # Project to XZ-plane
    kpi_angle = np.arctan2(kpi_vector_xz[0], kpi_vector_xz[1])  # Angle with vertical (z-axis)
    kpi_angle_deg = np.degrees(kpi_angle)
    
    # Camber Angle (YZ-plane)
    camber_vector_yz = np.array([kingpin_vector[1], kingpin_vector[2]])  # Project to YZ-plane
    camber_angle = np.arctan2(camber_vector_yz[0], camber_vector_yz[1])  # Angle with vertical (z-axis)
    camber_angle_deg = np.degrees(camber_angle)
    
    return {
        "trail": trail,
        "scrub": scrub,
        "kingpin_inclination_deg": kpi_angle_deg,
        "camber_angle_deg": camber_angle_deg,
    }

def bump_steer_calc(upper, lower, tie, previous_angle=None):
    """
    Calculate the bump steer angle change for a given suspension iteration.
    
    Parameters:
        outer_point (np.array): Current position of the outer point (3D vector).
        steering_axis_points (np.array): Two points defining the steering axis (2x3 array).
        previous_angle (float, optional): The previous steering angle in radians. If None, only returns the current angle.

    Returns:
        float: The change in steering angle (if previous_angle is provided).
        float: The current steering angle in radians.
    """
    
    # Define the steering axis direction vector
    axis_point1 = upper[0,:]
    axis_point2 = lower[0,:]
    steering_axis = axis_point2 - axis_point1
    steering_axis /= np.linalg.norm(steering_axis)  # Normalize the vector

    # Vector from the axis to the outer point
    vector_to_outer = tie[0,:] - axis_point1

    # Project the outer point vector onto the steering axis (normalized direction)
    projection = np.dot(vector_to_outer, steering_axis) * steering_axis

    # Perpendicular vector from the projection to the outer point
    perpendicular_vector = vector_to_outer - projection

    # Calculate the steering angle (angle between perpendicular vector and projection)
    angle = np.arctan2(
        np.linalg.norm(np.cross(steering_axis, vector_to_outer)),
        np.dot(steering_axis, vector_to_outer)
    )
    
    angle = np.rad2deg(angle)

    return angle

def force_calc(upper,lower,tie, force, location):
    """
    Calculate the force reaction at the outer points in xyz
    Assumptions:
    * All joints are revolute
    * Lower a-arm reacts contact force z
    * Upper & tie rod only react load in polar form (along their line of action)

    Args:
        upper (_matrix_): a-arm joint position matrix
        lower (_matrix_): a-arm joint position matrix
        tie (_matrix_): tie rod joint position matrix
        force (_vector_): Force vector
        location (_vector_): tire force cop location
    """
    #Create matrix of outter points
    outer = np.zeros((3,3))
    outer[0] = upper[0,:]
    outer[1] = lower[0,:]
    outer[2] = tie[0,:]
    
    # Moment arms of locations
    r = outer - location
    
    # Force equilibrium equations
    A = np.array([
        # Force equilibrium (ΣFx = 0)
        [1, 0,  1, 0,  1, 0,  0],  # Fx_upper, Fy_upper, Fx_lower, Fy_lower, Fx_tie, Fy_tie, Fz_lower
        [0, 1,  0, 1,  0, 1,  0],  # Fx_upper, Fy_upper, Fx_lower, Fy_lower, Fx_tie, Fy_tie, Fz_lower
        [0, 0,  0, 0,  0, 0,  1],  # Only Fz_lower (Z reaction assumed only at lower A-arm)
        # Moment equilibrium (ΣMx = 0, ΣMy = 0, ΣMz = 0)
        [0,  r[0,2],  0, r[1,2],  0, r[2,2], -r[1,1]],  # Mx = r × F (y * Fz - z * Fy)
        [-r[0,2],  0, -r[1,2], 0, -r[2,2], 0,  r[1,0]],  # My = r × F (z * Fx - x * Fz)
        [r[0,1], -r[0,0],  r[1,1], -r[1,0],  r[2,1], -r[2,0], 0]  # Mz = r × F (x * Fy - y * Fx)
    ], dtype=np.float64)

    # Right-hand side vector (force & moment equilibrium)
    b = np.array([force[0], force[1], force[2], 0, 0, 0], dtype=np.float64)

    # Solve the system Ax = b
    reactions, residuals, rank, s_values = np.linalg.lstsq(A, b, rcond=None)
    
    """print("Reaction Forces:")
    force_labels = [
        "Fx (Upper A-arm)", 
        "Fy (Upper A-arm)", 
        "Fx (Lower A-arm)", 
        "Fy (Lower A-arm)", 
        "Fx (Tie Rod)", 
        "Fy (Tie Rod)", 
        "Fz (Lower A-arm)"
    ]
    
    for label, value in zip(force_labels, reactions):
        print(f"{label}: {value}")
        """
    
    return reactions