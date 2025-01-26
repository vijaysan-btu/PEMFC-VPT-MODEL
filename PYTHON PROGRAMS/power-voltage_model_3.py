import sys
import tkinter as tk
from tkinter import messagebox, filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd
import math
import numpy as np

# Base Fuel Cell Model Class
class FuelCellModel:
    def __init__(self, name):
        self.name = name
        self.parameters = {}

    def set_parameters(self, params):
        """Set model parameters."""
        self.parameters.update(params)

    def simulate(self):
        """Simulate the model. To be implemented by subclasses."""
        raise NotImplementedError("This method should be implemented in subclasses.")

# Loss-Less Model
class LossLessModel(FuelCellModel):
    def __init__(self):
        super().__init__("Loss-Less")

    def simulate(self):
        """Simulate the Loss-Less model."""
        T = self.parameters.get("temperature", 298.15)  # Kelvin
        Ph2 = self.parameters.get("pressure_hydrogen", 1.0)  # atm
        PO2 = self.parameters.get("pressure_oxygen", 1.0)  # atm
        Reference_Voltage = 1.229
        Nernst_Voltage = Reference_Voltage - ((8.5e-4) * (T - 298.15)) + ((4.308e-5) * T * (np.log10(Ph2) + np.log10(PO2)))
        return {"Voltage": Nernst_Voltage}

# Loss Model
class LossModel(FuelCellModel):
    def __init__(self, lossless_model):
        super().__init__("Loss")
        self.lossless_model = lossless_model

    def simulate(self):
        """Simulate the Loss model."""
        T = self.parameters.get("temperature", 298.15)  # Kelvin
        Ph2 = self.parameters.get("pressure_hydrogen", 1.0)  # atm
        PO2 = self.parameters.get("pressure_oxygen", 1.0)  # atm
        A = self.parameters.get("active_area", 1.0)  # cm²
        i = self.parameters.get("current_density", 1.0)  # A/cm²
        lamb = 23  # Always use lambda as 23
        membrane_thickness = self.parameters.get("membrane_thickness", 0.1)  # cm
        num_cells = self.parameters.get("number_of_cells", 1)

        # Loss-less voltage calculation
        lossless_voltage = self.lossless_model.simulate()["Voltage"]

        # Activation loss calculation
        Ch2 = Ph2 / (1.09e6 * np.exp(77 / T))
        Co2 = PO2 / (5.08e6 * np.exp(-498 / T))
        yeta_1 = -0.948
        yeta_2 = 0.00286 + 0.002 * np.log10(A) + (4.3e-5) * np.log10(Ch2)
        yeta_3 = 7.6e-5
        yeta_4 = -1.93e-4
        activation_loss = yeta_1 + (yeta_2 * T) + (yeta_3 * T * np.log10(Co2)) + (yeta_4 * T * np.log10(i))

        # Ohmic loss calculation
        c_a = i / A
        ep = (4.18 * ((T - 303) / T))
        Rho = (181.6 * (1 + (0.03 * c_a) + (0.062 * ((T / 303) ** 2) * (c_a ** 2.5)))) / (lamb - 0.634 - (3 * c_a * np.exp(ep)))
        R_proton = (Rho * membrane_thickness) / A
        ohmic_loss = i * R_proton

        # Concentration loss calculation
        B = 0.015
        J_max = 1.5
        J = c_a
        concentration_loss = -B * np.log10(1 - (J / J_max))

        # Total loss
        total_loss = activation_loss + ohmic_loss + concentration_loss

        # Voltage calculation
        voltage = lossless_voltage - total_loss
        power_cell = voltage * i
        power_stack = power_cell * num_cells
        return {"Voltage": voltage, "Power_Cell": power_cell, "Power_Stack": power_stack}

# Combined Fuel Cell Model
class CombinedFuelCellModel:
    def __init__(self):
        self.lossless_model = LossLessModel()
        self.models = {
            "Loss": LossModel(self.lossless_model),
        }

    def run_simulation(self, params):
        """Run simulation for models and return results."""
        results = {}
        for name, model in self.models.items():
            model.set_parameters(params)
            results[name] = model.simulate()
        return results

# GUI Application using tkinter
class FuelCellApp:
    def __init__(self, master):
        self.master = master
        self.master.title("Fuel Cell Model Simulator")

        # Combined model instance
        self.combined_model = CombinedFuelCellModel()

        # Frame for input fields
        self.input_frame = tk.Frame(master)
        self.input_frame.pack(pady=10)

        # Input fields with equal spacing
        self.current_density_input = self.add_input("Current Density (A/cm²):", 1.0)
        self.temperature_input = self.add_input("Temperature (K):", 298.15)
        self.pressure_hydrogen_input = self.add_input("Hydrogen Pressure (atm):", 1.0)
        self.pressure_oxygen_input = self.add_input("Oxygen Pressure (atm):", 1.0)
        self.active_area_input = self.add_input("Active Area (cm²):", 1.0)
        self.membrane_thickness_input = self.add_input("Membrane Thickness (cm):", 0.1)
        self.number_of_cells_input = self.add_input("Number of Cells:", 1)

        # Run button
        self.run_button = tk.Button(master, text="Run Simulation", command=self.run_simulation)
        self.run_button.pack(pady=5)

        # Save graph button
        self.save_button = tk.Button(master, text="Save Graph", command=self.save_graph)
        self.save_button.pack(pady=5)

        # Frame for plots
        self.plot_frame = tk.Frame(master)
        self.plot_frame.pack(pady=10)

        # Create initial empty figures
        self.fig, (self.ax_voltage, self.ax_power) = plt.subplots(1, 2, figsize=(10, 5))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack()

    def add_input(self, label_text, default_value):
        """Helper method to add input fields with equal spacing."""
        frame = tk.Frame(self.input_frame)
        label = tk.Label(frame, text=label_text, width=25, anchor='w')  # Fixed width for labels
        entry = tk.Entry(frame, width=10)  # Fixed width for entries
        entry.insert(0, str(default_value))
        label.pack(side=tk.LEFT, padx=5)
        entry.pack(side=tk.RIGHT, padx=5)
        frame.pack(side=tk.TOP, pady=5)
        return entry

    def run_simulation(self):
        """Run the simulation and update plots."""
        try:
            # Gather input parameters
            params = {
                "current_density": float(self.current_density_input.get()),
                "temperature": float(self.temperature_input.get()),
                "pressure_hydrogen": float(self.pressure_hydrogen_input.get()),
                "pressure_oxygen": float(self.pressure_oxygen_input.get()),
                "active_area": float(self.active_area_input.get()),
                "membrane_thickness": float(self.membrane_thickness_input.get()),
                "number_of_cells": int(self.number_of_cells_input.get()),
            }

            # Simulate at multiple current densities
            current_densities = [i / 10 for i in range(1, 101)]
            model_results = {"Loss": []}

            for cd in current_densities:
                params["current_density"] = cd
                results = self.combined_model.run_simulation(params)
                for name, result in results.items():
                    model_results[name].append((cd, result["Voltage"], result["Power_Stack"]))

            # Update plots
            self.update_plots(model_results)

        except Exception as e:
            messagebox.showerror("Error", str(e))

    def update_plots(self, model_results):
        self.ax_voltage.clear()
        self.ax_power.clear()

        # Voltage vs Current Density for Loss Model
        results = model_results["Loss"]
        current_densities = [res[0] for res in results]
        voltages = [res[1] for res in results]
        powers = [res[2] for res in results]

        self.ax_voltage.plot(current_densities, voltages, label="Loss Model")
        self.ax_voltage.set_title("Voltage vs Current Density")
        self.ax_voltage.set_xlabel("Current Density (A/cm²)")
        self.ax_voltage.set_ylabel("Voltage (V)")
        self.ax_voltage.legend()

        # Power vs Voltage for Loss Model
        self.ax_power.plot(voltages, powers, label="Loss Model")
        self.ax_power.set_title("Power vs Voltage")
        self.ax_power.set_xlabel("Voltage (V)")
        self.ax_power.set_ylabel("Power (W)")
        self.ax_power.legend()

        self.canvas.draw()

    def save_graph(self):
        """Save the graph data and image to an Excel file."""
        file_path = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")])
        if file_path:
            try:
                # Prepare data for saving
                current_densities = [i / 10 for i in range(1, 101)]
                data = []

                for cd in current_densities:
                    params = {
                        "current_density": cd,
                        "temperature": float(self.temperature_input.get()),
                        "pressure_hydrogen": float(self.pressure_hydrogen_input.get()),
                        "pressure_oxygen": float(self.pressure_oxygen_input.get()),
                        "active_area": float(self.active_area_input.get()),
                        "membrane_thickness": float(self.membrane_thickness_input.get()),
                        "number_of_cells": int(self.number_of_cells_input.get()),
                    }
                    result = self.combined_model.run_simulation(params)["Loss"]
                    data.append([cd, result["Voltage"], result["Power_Stack"]])

                # Create a DataFrame and save to Excel
                df = pd.DataFrame(data, columns=["Current Density (A/cm²)", "Voltage (V)", "Power (W)"])
                df.to_excel(file_path, index=False)

                # Save the figure
                self.fig.savefig(file_path.replace('.xlsx', '.png'))

                messagebox.showinfo("Success", "Graph and data saved successfully.")
            except Exception as e:
                messagebox.showerror("Error", str(e))

if __name__ == "__main__":
    root = tk.Tk()
    app = FuelCellApp(root)
    root.mainloop()
