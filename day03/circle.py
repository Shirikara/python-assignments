# ----- This program prints out the area and the circumference of the circle in an interactive way with the user -----#

# Libraries
import argparse
from math import pi

# Welcome message
print("Welcome to the Circle Program!")
print("This program calculates the area and circumference of a circle based on the radius you provide.\n")

# Set up the argument parser
parser = argparse.ArgumentParser(description="Calculate the area and circumference of a circle.")
parser.add_argument('--radius', help="Radius of the circle in cm (must be a positive number)", required=True, type=float)

try:
    # Parse the arguments
    args = parser.parse_args()
    radius = args.radius

    # Validate the radius
    if radius <= 0:
        raise ValueError("Radius must be a positive number greater than zero.")

    # Calculate the area of the circle
    circle_area = pi * radius ** 2

    # Calculate the circumference of the circle
    circle_circumference = 2 * pi * radius

    # Print the results
    print(f"\nResults:")
    print(f"The radius you entered is: {radius:.2f} cm")
    print(f"The area of the circle is: {circle_area:.2f} cm²")
    print(f"The circumference of the circle is: {circle_circumference:.2f} cm")

except ValueError as e:
    # Handle invalid radius input
    print(f"\nError: {e}")
except Exception as e:
    # Handle any unexpected errors
    print(f"\nAn unexpected error occurred: {e}")