# NREL5MW_Helix_LiDAR (!!!)
Code for master thesis: LiDAR enhanced Closed-Loop Wake Mixing Control. Please note that the object wind turbine has been changed from IEA15MW to NREL5MW for the sack of faster simulation.

## File Description
### Materials
Reference paper, slides, etc

### Data
#### Figures&Video
Figures and videos are used for presentations and demonstrations.
#### NREL5MW
Data for system identification, system dynamic study, and so on.

### Functions
Matlab function
* BetaCenter_Comparison_Visualization: Data processing for visualizing the LiDAR sampling; Comparison of visualization.
* BetaCenter_Visualization: Visualization code without the comparison functionality.
* calculateBandwidth: Calculate the bandwidth (estimated) of the given system.
* Circle_LiDAR: The assumed perfect LiDAR (seeing the information of the entire plant) for sampling.
* Circle_LiDAR_Parallel: Optimized version of Circle_LiDAR that integrates GPU computing, very fast.
* Circle_LiDAR_Parallel_WakeCetner: Optimized version of Circle_LiDAR_Parallel with real-time wake center calculation functionality.
* FFT_func: Implemented Fast Fourier Transform for a given signal with visualization in the frequency domain.
* HelixCenter: Function that calculates the helix center by brutally filtering out the highspeed region.
* lowpassFilter: As the name implies, stop the high-frequency signal from passing.
* ringVisualization: As the name implies already.
* slideWindow: Create a sliding window that is useful when implementing the low-pass filter in real-time.
* spa_avf: Spectral analysis with frequency averaging, is useful for comparing the frequency response of the identified model and the real model.
* testCoupling: Compute the RGA of a given MIMO system at steady-state frequency and bandwidth frequency.
* videoCompare_func: Visualize the LiDAR sampling as a video and also give the option to save the video for demonstration.
* wakeCenterTraj: Study the trajectory of the wake center.
* writeDisconIPC: This file has been here forever, it is so old that I can't remember what the fuck does it do.
* ZXTM_LiDAR: LiDAR ring sampling.
* ZXTM_LiDAR_Parallel: Optimized version of ZXTM_LiDAR that integrates GPU computing, very fast.

### Source
This file contains the project and simulation files for NREL5MW and IEA15MW wind turbines
* .qpr: Wind turbine project file
* .sim: Simulation file for turbine project with simulation definitions
* Folders: Supporting files for simulation

### Matlab.m Files
* CLctrl_controllerDesign_MIMO: MIMO controller design based on the identified model. A fundamental step before implementing in the loop.
* CLctrl_controllerDesign_SISO: SISO controller design based on the identified model. A fundamental step before implementing in the loop.
* CLctrl_debug: Another file that made me want to crash my computer many many many times.
* CLctrl_debug_oppositeModel_MIMO_Hinf: H infinity controller for controlling the single channel.
* CLctrl_debug_oppositeModel_SISO_I: I controller for controlling the single channel.
* CLctrl_debug_oppositeModel_SISO_PI: PI controller for controlling the single channel.
* coordinateTransform: Study and explore the different coordinate transforms including MBC transform, new transform, and so on.
* correlationStudy: A file that I spent a lot of time exploring to find the control input and output.
* debug: The file that made me want to crash my computer many many many times.
* decouple: Decouple the MIMO system for easier control (I wish).
* decouple_RealApproximation: Conduct a real approximation of the system's response at a certain frequency for building a decouple compensator.
* delay: Add a delay factor to the system for building the Smith predictor.
* ForMarion: Marion is from Belgian.
* ForMarionFigure: Marion is a postdoc.
* Helix_IEA15MW_fixedFrame: Simulation of IEA15MW wind turbine in the fixed frame.
* Helix_IEA15MW_helixFrame: Simulation of IEA15MW wind turbine in the helix frame.
* Helix_NREL5MW_fixedFrame: Simulation of NREL5MW wind turbine in the fixed frame.
* Helix_NREL5MW_helixFrame: Simulation of NREL5MW wind turbine in the helix frame.
* Helix_NREL5MW_helixFrame_CLctrl: Closed-loop control of the Helix for NREL5MW wind turbine in the helix frame.
* HelixCenter_study: Study the helix center data sampled from the LiDAR.
* sysIDE: As the name implies, the file that conducts system identification.
* sysIDE_dataAcquision_NREL5MW: Simulation for acquiring the data for system identification.

### .qpr Files
Files that store the simulation result

### QBladeFunctionHelp
Doc for using q-blade functions in Matlab
