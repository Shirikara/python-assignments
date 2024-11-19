# ----- This program prints out the area and the circumference of the rectangle in an interactive way with the user -----#

# Prompt the user to define width
rec_width = float(input("Please type the width of the rectangle: "))

# Prompt the user to define length
rec_length = float(input("Please type the length of the rectangle: "))

# Prinet the area of the rectangle
rec_area = rec_width * rec_length
print("The area of the rectangle is: ", rec_area)

# Pritne the circumference of the rectangle
rec_circumference = 2*rec_width + 2*rec_length
print("The circumference of the rectangle is: ", rec_circumference)
