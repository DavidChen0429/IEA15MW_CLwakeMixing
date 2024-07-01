This folder stores the windspeed data from LiDAR sampling in .mat format. the naming of each data file looks like the framework below:


(turbine name)_(control method)_(average inflow speed)_(turbulence situation)_(simulation time in seconds)_(LiDAR sampling range)_windspeed

- Turbine Name: IEA15MW mainly for this thesis work
- Control Method: TC(Torque Control), Helix(Helix + TC)
  - When Helix is activated, additional properties like CCW(Counterclockwise) and CW(Clockwise), as well as Strhoual number(Str) are recorded
- Average Inflow Speed: Very straightforward, in m/s
- Turbulence Situation:
  - Uni: Uniform
  - NTM: Normal Turbulence Model
  - ETM: Extreme Turbulence Model
- Simulation Time in Seconds: As the name implies
- LiDAR Sampling Range: Snapshot where the continuous-wave LiDAR samples data, "d" stands for downwind position. "D" stands for 1 rotor diameter 
