'''
 This code predicts the protein concentration given 
a spectrophotometer readouts in Fluorescamine Assay using tkinter (Tk) GUI

Note: this code is based on the assignment in Day03 in this Python course
'''

import tkinter as tk
from tkinter import messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

# Dummy calibration function
def calibrate(rfu_values):
    try:
        # Example: Linear relationship for simplicity
        rfu_values = [float(val) for val in rfu_values.split(",")]
        concentrations = [0.5 * val for val in rfu_values]  # Example calculation
        return concentrations
    except ValueError:
        messagebox.showerror("Input Error", "Please enter numeric values separated by commas.")
        return None

# Function to plot the calibration curve
def plot_curve():
    # Example RFU and concentration data
    rfu = np.linspace(0, 100, 50)
    concentration = 0.5 * rfu  # Example linear relationship

    # Create the plot
    plt.figure(figsize=(6, 4))
    plt.plot(rfu, concentration, label="Calibration Curve", color="blue")
    plt.title("Protein Quantification Calibration Curve")
    plt.xlabel("RFU")
    plt.ylabel("Concentration")
    plt.grid(True)
    plt.legend()

    # Display the plot in the Tkinter app
    fig = plt.gcf()
    canvas = FigureCanvasTkAgg(fig, master=frame_plot)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

# Function to predict concentrations from input RFU values
def predict():
    input_rfu = entry_rfu.get()
    concentrations = calibrate(input_rfu)
    if concentrations:
        output = "\n".join([f"RFU: {rfu}, Concentration: {conc}" for rfu, conc in zip(input_rfu.split(","), concentrations)])
        text_output.delete(1.0, tk.END)
        text_output.insert(tk.END, output)

# Initialize Tkinter App
app = tk.Tk()
app.title("Protein Quantification")

# Frame for the plot
frame_plot = tk.Frame(app, bd=2, relief=tk.SUNKEN)
frame_plot.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

# Input Frame
frame_input = tk.Frame(app)
frame_input.pack(side=tk.TOP, pady=10)

# Input label and text entry
label_rfu = tk.Label(frame_input, text="Enter RFU values (comma-separated):")
label_rfu.pack(side=tk.LEFT, padx=5)
entry_rfu = tk.Entry(frame_input, width=50)
entry_rfu.pack(side=tk.LEFT, padx=5)

# Predict Button
btn_predict = tk.Button(frame_input, text="Predict", command=predict)
btn_predict.pack(side=tk.LEFT, padx=5)

# Output Frame
frame_output = tk.Frame(app)
frame_output.pack(side=tk.TOP, pady=10)

# Output Text Box
text_output = tk.Text(frame_output, height=10, width=60)
text_output.pack()

# Plot Button
btn_plot = tk.Button(app, text="Show Calibration Curve", command=plot_curve)
btn_plot.pack(side=tk.TOP, pady=5)

# Run the Tkinter app
app.mainloop()
