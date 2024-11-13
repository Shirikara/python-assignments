# ----- This program is about a circle -----#

# Libraries
from math import pi

# Prompt the user to input the radius of the circle
radius = float(input("Please type the radius of the circle: "))

# Calculate the area of the circle (formula = π * r^2)
circle_area = pi * radius ** 2

# Print the result, including the radius and calculated area
print("The area of the circle is: " + str(round(circle_area,2)))

# Calculate the area of the circle (formula = 2πr)

# Print the result, including the radius and calculated circumference
circle_circumference = 2 * pi * radius
print("The circumference of the circle is: " + str(round(circle_circumference,2)))


