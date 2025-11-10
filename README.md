# wind-assignment-2
## Electro assignment - Wind Farm

### Code structure
***main.ipynb*** is the main notebook.

***toolbox.py*** hold class **Tools** which is a collection of functions we can call and use throughout the main script.

***equations.ipynb*** has the derivation of the three equations from the grid (right) side of the system.

### Tasks
- [x] Code how to solve for power on generator side
  - [ ] Rewrite into callable function, inputs are windspeed and/or omega, put function in toolbox.
- [ ] Code how to solve for power, second method.
- [x] Derive three equations on the grid side.
  - [ ] Find out a way to verify the three functions, maybe speak with Nenad?
- [ ] Write function in toolbox.py to output all impedances in system **Tools.get_impedances()**
- [ ] Implement solver (somewhat implemented, verify if it works)
- [ ] Put the whole dang thing into a for loop to compute over range of wind speeds.
- [ ] Also, write a report

### Questions
##### *Should we use single-phase power as input to the right hand side of the system?*
AVE: I think we solve the whole thing in single phase and multiply as the very last thing to get three phase values, I'm not sure though.