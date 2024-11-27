########################
#Script to install depdencies 
########################
#Author: Alejandro Cardenas-Avendano

# Import the Pkg module
import Pkg

# List of required packages
packages = [
    "HDF5",          # HDF5 file handling
    "ArgParse",      # Argument parsing
    "GSL",           # GNU Scientific Library bindings
    "Plots",         # Plotting library
    "LaTeXStrings"   # LaTeX formatting for strings in plots
]

# Install missing packages
for pkg in packages
    if !(pkg in keys(Pkg.installed()))
        println("Installing ", pkg, "...")
        Pkg.add(pkg)
    else
        println(pkg, " is already installed.")
    end
end

println("All dependencies are installed.")

# Optional: Update all packages to the latest version
Pkg.update()
