import sys
import tkinter as tk
from tkinter import messagebox, filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd

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
        current_density = self.parameters.get("current_density", 1.0)
        temperature = self.parameters.get("temperature", 298.15)  # Kelvin (Typocal Room temperature of 25 degC)
        pressure = self.parameters.get("pressure", 1.0)  # atm (Atmospheric Pressure)
        voltage = (0.6 - 0.1 * current_density + 0.001 * (temperature - 298.15) - 0.01 * (pressure - 1.0))  # Simplified
        power = voltage * current_density
        return {"Voltage": voltage, "Power": power}

# Loss Model
class LossModel(FuelCellModel):
    def __init__(self):
        super().__init__("Loss")

    def simulate(self):
        """Simulate the Loss model."""
        current_density = self.parameters.get("current_density", 1.0)
        activation_loss = self.parameters.get("activation_loss", 0.05) * current_density
        ohmic_loss = self.parameters.get("ohmic_loss", 0.02) * current_density
        concentration_loss = self.parameters.get("concentration_loss", 0.01) * current_density
        voltage = 0.7 - (activation_loss + ohmic_loss + concentration_loss)
        power = voltage * current_density
        return {"Voltage": voltage, "Power": power}

# Combined Fuel Cell Model
class CombinedFuelCellModel:
    def __init__(self):
        self.models = {
            "Loss-Less": LossLessModel(),
            "Loss": LossModel(),
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
        self.pressure_input = self.add_input("Pressure (atm):", 1.0)
        self.activation_loss_input = self.add_input("Activation Loss (V):", 0.05)
        self.ohmic_loss_input = self.add_input("Ohmic Loss (V):", 0.02)
        self.concentration_loss_input = self.add_input("Concentration Loss (V):", 0.01)

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
        self.fig, (self.ax_power, self.ax_voltage) = plt.subplots(1, 2, figsize=(10, 5))
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
                "pressure": float(self.pressure_input.get()),
                "activation_loss": float(self.activation_loss_input.get()),
                "ohmic_loss": float(self.ohmic_loss_input.get()),
                "concentration_loss": float(self.concentration_loss_input.get()),
            }

            # Simulate at multiple current densities
            current_densities = [i / 10 for i in range(1, 101)]
            model_results = {name: [] for name in self.combined_model.models.keys()}

            for cd in current_densities:
                params["current_density"] = cd
                results = self.combined_model.run_simulation(params)
                for name, result in results.items():
                    model_results[name].append((cd, result["Voltage"], result["Power"]))

            # Update plots
            self.update_plots(model_results)

        except Exception as e:
            messagebox.showerror("Error", str(e))

    def update_plots(self, model_results):
        self.ax_power.clear()
        self.ax_voltage.clear()

        # Power vs Voltage
        for name, results in model_results.items():
            voltages = [res[1] for res in results]
            powers = [res[2] for res in results]
            self.ax_power.plot(voltages, powers, label=f"{name} Model")
        self.ax_power.set_title("Power vs Voltage")
        self.ax_power.set_xlabel("Voltage (V)")
        self.ax_power.set_ylabel("Power (W)")
        self.ax_power.legend()

        # Voltage vs Current Density
        for name, results in model_results.items():
            current_densities = [res[0] for res in results]
            voltages = [res[1] for res in results]
            self.ax_voltage.plot(current_densities, voltages, label=f"{name} Model")
        self.ax_voltage.set_title("Voltage vs Current Density")
        self.ax_voltage.set_xlabel("Current Density (A/cm²)")
        self.ax_voltage.set_ylabel("Voltage (V)")
        self.ax_voltage.legend()

        self.canvas.draw()

    def save_graph(self):
        """Save the graph data and image to an Excel file."""
        file_path = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")])
        if file_path:
            try:
                # Prepare data for saving
                current_densities = [i / 10 for i in range(1, 101)]
                data = []

                for name, model in self.combined_model.models.items():
                    for cd in current_densities:
                        params = {
                            "current_density": cd,
                            "temperature": float(self.temperature_input.get()),
                            "pressure": float(self.pressure_input.get()),
                            "activation_loss": float(self.activation_loss_input.get()),
                            "ohmic_loss": float(self.ohmic_loss_input.get()),
                            "concentration_loss": float(self.concentration_loss_input.get()),
                        }
                        model.set_parameters(params)
                        result = model.simulate()
                        data.append([name, cd, result["Voltage"], result["Power"]])

                # Create a DataFrame and save to Excel
                df = pd.DataFrame(data, columns=["Model", "Current Density (A/cm²)", "Voltage (V)", "Power (W)"])
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
