# PHEV-TSPS
Codes and instances used for Paper 'Modeling and Solving the Traveling Salesman Problem for a Plug-in Hybrid Electric Vehicle Considering Speed Optimization"

# Codes
All codes were implemented in Python using the Gurobi solver. 

para.py : parameters used for the model; 

phev_tsps.py : the solution method for the original PHEV-TSPS without valid inequalities;

phev_tsps_vi.py : the solution method for the PHEV-TSPS with valid inequalities;

phev_tspsd_vi.py : the solution method for the PHEV-TSPSD with valid inequalities;

phev_tsps_cha_vi.py : the solution method for the PHEV-TSPS-CS with valid inequalities.

# Inst
The instances with 8, 20, 50 customers were directly from 'Doppstadt C, Koberstein A, Vigo D, 2019 The hybrid electric vehicle—traveling salesman problem with time windows instances. Mendeley Dataset, URL doi:10.17632/9j3tt84hyx.1.'.

The instances with 30 and 40 customers were chosen from the instances with 50 customers.

The instances with 60, 70, and 80 customers are derived from the instances with 50 customers by adding additional customers who are randomly generated within the delivery area.

Note: The elevation of each customer location is randomly generated between 0 and 100 meters

# Inst2
These instances are similar to Inst but differ in terms of customer elevations. 

For the instances with the last number 2, customer elevations are randomly generated between 0 and 200 meters, while for the instances with the last number 3, customer elevations are randomly generated between 0 and 300 meters.
