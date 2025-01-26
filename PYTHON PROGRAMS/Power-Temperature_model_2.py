import sys
import tkinter as tk
from tkinter import messagebox, filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd
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

# Temperature Model
class TemperatureModel:
    def __init__(self):
        pass

    def calculate_temperature(self, T_initial, Q_generated, Q_dissipated, C_p, dt):
        """Calculate the new temperature based on heat generation and dissipation."""
        return T_initial + (Q_generated - Q_dissipated) * dt / C_p

# Combined Fuel Cell Model
class CombinedFuelCellModel:
    def __init__(self):
        self.lossless_model = LossLessModel()
        self.loss_model = LossModel(self.lossless_model)
        self.temperature_model = TemperatureModel()

    def run_simulation(self, params):
        """Run simulation for models and return results."""
        results = []

        # Initialize parameters
        T = params.get("temperature", 298.15)  # Initial temperature
        C_p = params.get("heat_capacity", 1000.0)  # Heat capacity (J/K)
        k = params.get("heat_transfer_coefficient", 10.0)  # Heat transfer coefficient (W/K)
        dt = params.get("time_step", 1.0)  # Time step (s)

        for time_step in range(params.get("simulation_steps", 100)):
            params["temperature"] = T

            # Run loss model simulation
            self.loss_model.set_parameters(params)
            loss_results = self.loss_model.simulate()

            # Heat generation and dissipation
            Q_generated = params["current_density"] * loss_results["Voltage"] + loss_results["Power_Cell"]
            Q_dissipated = k * (T - params.get("ambient_temperature", 298.15))

            # Update temperature
            T = self.temperature_model.calculate_temperature(T, Q_generated, Q_dissipated, C_p, dt)

            # Append results
            results.append({
                "Time": time_step * dt,
                "Temperature": T,
                "Voltage": loss_results["Voltage"],
                "Power_Stack": loss_results["Power_Stack"],
                "Current_Density": params["current_density"]
            })

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
        self.add_input("Fuel Flow Rate (L/min):", 10)
        self.add_input("Air Flow Rate (L/min):", 20)
        self.add_input("Pressure (atm):", 1.0)
        self.add_input("Temperature (K):", 300)
        self.add_input("Relative Humidity (%):", 50)

        # New input fields for graph size
        self.fig_width_input = self.add_input("Figure Width (inches):", 14)
        self.fig_height_input = self.add_input("Figure Height (inches):", 12)

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
        self.fig, self.axes = plt.subplots(2, 2, figsize=(14, 12))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack()

        # Add control panel for adjusting grid layout
        self.grid_controls = tk.Frame(master)
        self.grid_controls.pack(pady=10)

        tk.Label(self.grid_controls, text="Adjust Graph Layout:").grid(row=0, column=0, columnspan=2)
        tk.Label(self.grid_controls, text="Rows:").grid(row=1, column=0)
        self.grid_rows = tk.Entry(self.grid_controls, width=5)
        self.grid_rows.insert(0, "2")
        self.grid_rows.grid(row=1, column=1)

        tk.Label(self.grid_controls, text="Columns:").grid(row=2, column=0)
        self.grid_cols = tk.Entry(self.grid_controls, width=5)
        self.grid_cols.insert(0, "2")
        self.grid_cols.grid(row=2, column=1)

        tk.Button(self.grid_controls, text="Update Layout", command=self.update_layout).grid(row=3, column=0, columnspan=2)

    def add_input(self, label_text, default_value):
        """Helper method to add input fields with equal spacing."""
        label = tk.Label(self.input_frame, text=label_text)
        label.pack(side=tk.LEFT, padx=5)
        entry = tk.Entry(self.input_frame, width=10)
        entry.insert(0, str(default_value))
        entry.pack(side=tk.LEFT, padx=5)
        return entry

    def run_simulation(self):
        """Run the simulation and update plots."""
        # Gather input parameters
        params = {
            "Fuel_Flow_Rate": float(self.input_frame.children["!entry"].get()),
            "Air_Flow_Rate": float(self.input_frame.children["!entry2"].get()),
            "Pressure": float(self.input_frame.children["!entry3"].get()),
            "Temperature": float(self.input_frame.children["!entry4"].get()),
            "Relative_Humidity": float(self.input_frame.children["!entry5"].get()),
        }

        # Run simulation
        simulation_results = self.combined_model.run_simulation(params)

        # Update plots
        self.update_plots(simulation_results)

    def update_plots(self, simulation_results):
        """Update the plots based on simulation results."""
        self.fig.clear()
        rows = int(self.grid_rows.get())
        cols = int(self.grid_cols.get())
        
        # Get figure size from inputs
        fig_width = float(self.fig_width_input.get())
        fig_height = float(self.fig_height_input.get())
        self.fig, self.axes = plt.subplots(rows, cols, figsize=(fig_width, fig_height))

        # Ensure axes is iterable
        if rows * cols == 1:
            self.axes = [self.axes]
        else:
            self.axes = self.axes.flatten()

        # Extract results
        temperatures = [res["Temperature"] for res in simulation_results]
        powers = [res["Power_Stack"] for res in simulation_results]
        voltages = [res["Voltage"] for res in simulation_results]
        current_densities = [res["Current_Density"] for res in simulation_results]

        # Plot data
        plot_data = [
            (current_densities, voltages, "Voltage vs Current Density", "Current Density (A/cm²)", "Voltage (V)"),
            (temperatures, powers, "Power vs Temperature", "Temperature (K)", "Power (W)"),
            (voltages, powers, "Power vs Voltage", "Voltage (V)", "Power (W)"),
            (temperatures, voltages, "Voltage vs Temperature", "Temperature (K)", "Voltage (V)"),
        ]

        for ax, data in zip(self.axes, plot_data):
            x, y, title, xlabel, ylabel = data
            ax.plot(x, y)
            ax.set_title(title)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

        self.fig.tight_layout()
        self.canvas.draw()

    def update_layout(self):
        """Update the grid layout of the plots."""
        self.update_plots([])  # Refresh the plots

    def save_graph(self):
        """Save the graph data and image to an Excel file."""
        # Get the current figure
        fig = self.fig

        # Ask the user for a file name
        file_name = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel files", "*.xlsx")])
        if file_name:
            # Create a new Excel workbook
            workbook = pd.ExcelWriter(file_name, engine="xlsxwriter")

            # Save the figure to a BytesIO object
            fig.savefig(workbook.path.replace(".xlsx", ".png"))

            # Create a DataFrame with the plot data
            plot_data = [
                (current_densities, voltages, "Voltage vs Current Density"),
                (temperatures, powers, "Power vs Temperature"),
                (voltages, powers, "Power vs Voltage"),
                (temperatures, voltages, "Voltage vs Temperature"),
            ]
            for x, y, title in plot_data:
                df = pd.DataFrame({"X": x, "Y": y})
                df.to_excel(workbook, sheet_name=title, index=False)

            # Save the workbook
            workbook.save()
            messagebox.showinfo("Graph Saved", f"Graph data and image saved to {file_name}")

# Main execution
if __name__ == "__main__":
    root = tk.Tk()
    app = FuelCellApp(root)
    root.mainloop()