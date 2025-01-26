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
        return {"Voltage": voltage, "Current Density": current_density, "Temperature": temperature, "Pressure": pressure}

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
        return {"Voltage": voltage, "Current Density": current_density, "Activation Loss": activation_loss, "Ohmic Loss": ohmic_loss, "Concentration Loss": concentration_loss}

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

        # Activation loss input (Loss specific)
        ttk.Label(root, text="Activation Loss (V):").grid(row=3, column=0, padx=10, pady=10)
        activation_loss_entry = ttk.Entry(root, textvariable=self.activation_loss_var)
        activation_loss_entry.grid(row=3, column=1, padx=10, pady=10)

        # Ohmic loss input (Loss specific)
        ttk.Label(root, text="Ohmic Loss (V):").grid(row=4, column=0, padx=10, pady=10)
        ohmic_loss_entry = ttk.Entry(root, textvariable=self.ohmic_loss_var)
        ohmic_loss_entry.grid(row=4, column=1, padx=10, pady=10)

        # Concentration loss input (Loss specific)
        ttk.Label(root, text="Concentration Loss (V):").grid(row=5, column=0, padx=10, pady=10)
        concentration_loss_entry = ttk.Entry(root, textvariable=self.concentration_loss_var)
        concentration_loss_entry.grid(row=5, column=1, padx=10, pady=10)

        # Run button
        run_button = ttk.Button(root, text="Run Simulation", command=self.run_simulation)
        run_button.grid(row=6, column=0, columnspan=2, pady=10)

        # Result display frame
        self.fig, self.ax = plt.subplots(figsize=(5, 4))
        self.ax.set_title("Fuel Cell Simulation Results")
        self.ax.set_xlabel("Current Density (A/cm^2)")
        self.ax.set_ylabel("Voltage (V)")
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.grid(row=7, column=0, columnspan=2, pady=10)

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

        try:
            # Simulate at multiple current densities for graph
            current_densities = [i / 10 for i in range(1, 101)]
            model_results = {name: [] for name in self.combined_model.models.keys()}

            for cd in current_densities:
                params["current_density"] = cd
                results = self.combined_model.run_simulation(params)
                for name, result in results.items():
                    model_results[name].append(result["Voltage"])

            # Update graph
            self.ax.clear()
            for name, voltages in model_results.items():
                self.ax.plot(current_densities, voltages, label=f"{name} Model")
            self.ax.set_title("Fuel Cell Simulation Results")
            self.ax.set_xlabel("Current Density (A/cm^2)")
            self.ax.set_ylabel("Voltage (V)")
            self.ax.legend()
            self.canvas.draw()
        except Exception as e:
            self.ax.clear()
            self.ax.text(0.5, 0.5, f"Error: {e}", fontsize=12, ha='center')
            self.canvas.draw()

# Run the application
if __name__ == "__main__":
    root = tk.Tk()
    app = FuelCellApp(root)
    root.mainloop()
