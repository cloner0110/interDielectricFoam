# interDielectricFoam: An Electrohydrodynamics Solver for OpenFOAM

![OpenFOAM Version](https://img.shields.io/badge/OpenFOAM-V2312-blue)

## Overview

**interDielectricFoam** is a custom solver based on OpenFOAM's `interFoam`, designed to simulate **electrohydrodynamic (EHD)** phenomena in multiphase flows. The solver accounts for the coupling between fluid dynamics, electric fields, and interfacial phenomena, enabling the study of dielectric fluids, electrically induced flows, and other EHD-related processes.

## Features

- **Multiphase Flow Solver**: Built upon `interFoam`, providing capabilities for tracking immiscible fluid interfaces using the Volume of Fluid (VOF) method.
- **Electrostatic Effects**: Incorporates electric field solvers to handle the interaction between fluids and electric charges.
- **Dielectric and Conductive Materials**: Supports simulations involving dielectric and weakly conductive fluids.
- **Electrocapillarity**: Models surface tension modifications due to electric fields.

## Applications

interDielectricFoam can be used for:

- Electrospraying and electrospinning processes.
- Taylor cone formation.
- Dielectrophoretic manipulation of fluids.
- Electrically driven breakup and coalescence of droplets.
- Electrocapillary phenomena at fluid interfaces.

## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/cloner0110/interDielectricFoam.git
   ```

2. **Compile the Solver**:
   Navigate to the solver directory and compile:
   ```bash
   cd src/interDielectricFoam
   wmake
   ```


<!-- ## Examples

Example cases are provided in the `tutorials/` directory. These include:

- **Electrospraying**: Demonstrates the formation of a Taylor cone.
- **Droplet Dynamics**: Simulates droplet coalescence under electric fields.
- **Dielectrophoresis**: Shows particle motion in response to non-uniform electric fields. -->

## Dependencies

- OpenFOAM (v2312)

## Contributing

Contributions are welcome! If you encounter bugs or have feature requests, feel free to open an issue or submit a pull request. Alternatively, you can contact me at: [amirhossein.taran@ucdconnect.ie](mailto:amirhossein.taran@ucdconnect.ie)

## License

This project is licensed under the [MIT License](LICENSE).

---

Enjoy simulating with **interDielectricFoam**!