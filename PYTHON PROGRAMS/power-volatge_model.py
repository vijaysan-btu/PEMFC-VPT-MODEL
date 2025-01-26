import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

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
        temperature = self.parameters.get("temperature", 298.15)  # Kelvin
        pressure = self.parameters.get("pressure", 1.0)  # atm
        voltage = 0.6 - 0.1 * current_density + 0.001 * (temperature - 298.15) - 0.01 * (pressure - 1.0)  # Simplified
        power = voltage * current_density
        return {"Voltage": voltage, "Power": power, "Current Density": current_density, "Temperature": temperature, "Pressure": pressure}

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
        return {"Voltage": voltage, "Power": power, "Current Density": current_density, "Activation Loss": activation_loss, "Ohmic Loss": ohmic_loss, "Concentration Loss": concentration_loss}

# Combined Fuel Cell Model
class CombinedFuelCellModel:
    def __init__(self):
        self.models = {
            "Loss-Less": LossLessModel(),
            "Loss": LossModel(),
        }

    def run_simulation(self, params):
        """Run simulation for all models and return results."""
        results = {}
        for name, model in self.models.items():
            model.set_parameters(params)
            results[name] = model.simulate()
        return results

# GUI Application
class FuelCellApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Fuel Cell Model Simulator")

        # Combined model instance
        self.combined_model = CombinedFuelCellModel()

        # GUI Variables
        self.current_density_var = tk.DoubleVar(value=1.0)
        self.temperature_var = tk.DoubleVar(value=298.15)  # Default: 298.15 K
        self.pressure_var = tk.DoubleVar(value=1.0)  # Default: 1 atm
        self.activation_loss_var = tk.DoubleVar(value=0.05)
        self.ohmic_loss_var = tk.DoubleVar(value=0.02)
        self.concentration_loss_var = tk.DoubleVar(value=0.01)
        self.num_cells_var = tk.IntVar(value=1)  # Number of cells in the stack

        # Current density input
        ttk.Label(root, text="Current Density (A/cm^2):").grid(row=0, column=0, padx=10, pady=10)
        current_density_entry = ttk.Entry(root, textvariable=self.current_density_var)
        current_density_entry.grid(row=0, column=1, padx=10, pady=10)

        # Temperature input
        ttk.Label(root, text="Temperature (K):").grid(row=1, column=0, padx=10, pady=10)
        temperature_entry = ttk.Entry(root, textvariable=self.temperature_var)
        temperature_entry.grid(row=1, column=1, padx=10, pady=10)

        # Pressure input
        ttk.Label(root, text="Pressure (atm):").grid(row=2, column=0, padx=10, pady=10)
        pressure_entry = ttk.Entry(root, textvariable=self.pressure_var)
        pressure_entry.grid(row=2, column=1, padx=10, pady=10)

        # Number of cells input
        ttk.Label(root, text="Number of Cells:").grid(row=3, column=0, padx=10, pady=10)
        num_cells_entry = ttk.Entry(root, textvariable=self.num_cells_var)
        num_cells_entry.grid(row=3, column=1, padx=10, pady=10)

        # Activation loss input (Loss specific)
        ttk.Label(root, text="Activation Loss (V):").grid(row=4, column=0, padx=10, pady=10)
        activation_loss_entry = ttk.Entry(root, textvariable=self.activation_loss_var)
        activation_loss_entry.grid(row=4, column=1, padx=10, pady=10)

        # Ohmic loss input (Loss specific)
        ttk.Label(root, text="Ohmic Loss (V):").grid(row=5, column=0, padx=10, pady=10)
        ohmic_loss_entry = ttk.Entry(root, textvariable=self.ohmic_loss_var)
        ohmic_loss_entry.grid(row=5, column=1, padx=10, pady=10)

        # Concentration loss input (Loss specific)
        ttk.Label(root, text="Concentration Loss (V):").grid(row=6, column=0, padx=10, pady=10)
        concentration_loss_entry = ttk.Entry(root, textvariable=self.concentration_loss_var)
        concentration_loss_entry.grid(row=6, column=1, padx=10, pady=10)

        # Run button
        run_button = ttk.Button(root, text="Run Simulation", command=self.run_simulation)
        run_button.grid(row=7, column=0, columnspan=2, pady=10)

        # Result display frames
        self.fig1, self.ax1 = plt.subplots(figsize=(5, 4))
        self.ax1.set_title("Power vs Voltage")
        self.ax1.set_xlabel("Voltage (V)")
        self.ax1.set_ylabel("Power (W)")
        self.canvas1 = FigureCanvasTkAgg(self.fig1, master=root)
        self.canvas1_widget = self.canvas1.get_tk_widget()
        self.canvas1_widget.grid(row=8, column=0, padx=10, pady=10)

        self.fig2, self.ax2 = plt.subplots(figsize=(5, 4))
        self.ax2.set_title("Voltage vs Current Density")
        self.ax2.set_xlabel("Current Density (A/cm^2)")
        self.ax2.set_ylabel("Voltage (V)")
        self.canvas2 = FigureCanvasTkAgg(self.fig2, master=root)
        self.canvas2_widget = self.canvas2.get_tk_widget()
        self.canvas2_widget.grid(row=8, column=1, padx=10, pady=10)

    def run_simulation(self):
        """Run the simulation for all models and display results."""
        params = {
            "current_density": self.current_density_var.get(),
            "temperature": self.temperature_var.get(),
            "pressure": self.pressure_var.get(),
            "activation_loss": self.activation_loss_var.get(),
            "ohmic_loss": self.ohmic_loss_var.get(),
            "concentration_loss": self.concentration_loss_var.get(),
        }
        num_cells = self.num_cells_var.get()

        try:
            # Simulate at multiple current densities for graphs
            current_densities = [i / 10 for i in range(1, 101)]
            model_results = {name: [] for name in self.combined_model.models.keys()}

            for cd in current_densities:
                params["current_density"] = cd
                results = self.combined_model.run_simulation(params)
                for name, result in results.items():
                    stack_voltage = result["Voltage"] * num_cells
                    stack_power = result["Power"] * num_cells
                    model_results[name].append((cd, result["Voltage"], result["Power"], stack_voltage, stack_power))

            # Update Power vs Voltage graph
            self.ax1.clear()
            for name, results in model_results.items():
                voltages = [res[1] for res in results]
                powers = [res[2] for res in results]
                self.ax1.plot(voltages, powers, label=f"{name} Model")
            self.ax1.set_title("Power vs Voltage")
            self.ax1.set_xlabel("Voltage (V)")
            self.ax1.set_ylabel("Power (W)")
            self.ax1.legend()
            self.canvas1.draw()

            # Update Voltage vs Current Density graph
            self.ax2.clear()
            for name, results in model_results.items():
                current_densities = [res[0] for res in results]
                voltages = [res[1] for res in results]
                self.ax2.plot(current_densities, voltages, label=f"{name} Model")
            self.ax2.set_title("Voltage vs Current Density")
            self.ax2.set_xlabel("Current Density (A/cm^2)")
            self.ax2.set_ylabel("Voltage (V)")
            self.ax2.legend()
            self.canvas2.draw()

        except Exception as e:
            self.ax1.clear()
            self.ax1.text(0.5, 0.5, f"Error: {e}", fontsize=12, ha='center')
            self.canvas1.draw()
            self.ax2.clear()
            self.ax2.text(0.5, 0.5, f"Error: {e}", fontsize=12, ha='center')
            self.canvas2.draw()

# Run the application
if __name__ == "__main__":
    root = tk.Tk()
    app = FuelCellApp(root)
    root.mainloop()