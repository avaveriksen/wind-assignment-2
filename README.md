# wind-assignment-2
## Electro assignment - Wind Farm

### Code structure
***main.ipynb*** is the main notebook.

***toolbox.py*** hold class **Tools** which is a collection of functions we can call and use throughout the main script.

***equations.ipynb*** has the derivation of the three equations from the grid (right) side of the system.

### Tasks
- [x] Code how to solve for power on generator side
  - [ ] Rewrite into callable function, inputs are windspeed and/or omega.
- [x] Derive three equations on the grid side.
  - [ ] Find out a way to verify the three functions, maybe speak with Nenad?
- [ ] Write function in toolbox.py to output all impedances in system **Tools.get_impedances()**
- [ ] Implement solver
- [ ] Put the whole dang thing into a for loop to compute over range of wind speeds.

