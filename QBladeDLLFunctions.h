/****************************************************************************

    QBladeDLLInclude Class
        Copyright (C) 2020 David Marten david.marten@qblade.org

*****************************************************************************/

///    PSEUDOCODE
///
///    the common order of execution should be:
///
///    createInstance()                     mandatory, create an instance of the QBladeEE DLL
///    loadProject() / loadSimDefinition()  mandatory, loads a .qpr project file or a .sim simulation definition file
///    addTurbulentWind()                   optional, if you want to create a turbulent wind field on the fly, overwrites sim definition
///    loadTurbulentWindBinary()            optional, instead of creating a windfield on the fly this loads an existing .bts file, overwrites sim definition
///    setPowerLawWind()                    optional, if you want to use laminar wind, overwrites sim definition
///    setInitialConditions_at_num(i)       optional, sets initial collective pitch, yaw and azimuthal position, overwrites sim definition, applied to the SELECTED turbine i
///    setTurbinePosition_at_num(i)         optional, if the turbine should "start" at a position/orientation other than (0,0,0,0,0,0), overwrites sim definition, applied to the SELECTED turbine i
///    setRPMPrescribeType_at_num(i)        optional, sets the prescribe rpm type, overwrites sim definition, applied to the SELECTED turbine i
///    setRampupTime()                      optional, sets the ramp-up time, overwrites sim definition
///    setTimestepSize()                    optional, overwrites timestep from sim definition
///
///    initializeSimulation()               mandatory!
///
///    for (int i=0;i<num;i++)
///    {
///       advanceTurbineSimulation()        mandatory, advances the aerodynamic and strucural simulation for ALL turbines
///       setPowerLawWind()                 optional, the wind input can be changed during the simulation to model gusts, wind ramps, etc...
///       advanceAero()                     mandatory, advances the aerodynamic simulation for all turbines
///       getCustomData_at_num(i)           optional, get the current value of an arbitrary data channel of the SeLECTED turbine i
///       advanceController_at_num(i)       optional, advances the controller for the SELECTED turbine i
///       getTurbineOperation_at_num(i)     optional, get turbine operational characteristics/loads of the SELECTED turbine i
///       setControlVars_at_num(i)          optional, overwrite the controller variables for the SELECTED turbine i, if not called the control vars are passed by the controller if a controller is part of the simulation
///       setTurbinePosition_at_num(i)      optional, if the turbine position/orientation should change, for cosimulation with hydrodynamics or other software, applied to the SELECTED turbine i
///    }
///
///    storeProject()                       (optional, stores the project as the specified file)


 void  loadProject(char *str);
//This function loads a qblade (.qpr) project into the dll

 void  loadSimDefinition(char *str);
//This function loads a simulation definition (.sim) file into the dll

 void  storeProject(char *str);
//This functions stores a project file

 void  setDebugInfo(bool isDebug);
//This function enables the debug info output if set to true

 void  createInstance(int clDevice, int groupSize);
//This function creates a new instance of QBlade

 void  closeInstance();
//This function creates a new instance of QBlade

 void  loadTurbulentWindBinary(char *str);
//This function allows to load a turbulent windfield that is stored in binary format

 void  addTurbulentWind(double windspeed, double refheight, double hubheight, double dimensions, int gridPoints, double length, double dT, char *turbulenceClass, char *turbulenceType, int seed, double vertInf, double horInf, bool removeFiles);
//This function allows to define and add a turbulent windfield to the simulations,m if a turbulent windfield is used the function 'setPowerLawWind' has no effect
//windspeed: the mean windspeed at the reference height [m/s]
//refheight: the reference height [m]
//hubheight: the hubheight, more specifically the height of the windfield center [m]
//dimensions: the y- and z- dimensions of the windfield in meters [m]
//gridpoints: the number of points in the y and z direction for which velocities are evaluated [-]
//length: the simulated length of the windfield in seconds [s]
//dT: the temporal resolution of the windfield [s]
//turbulenceClass: the turbulence class, can be "A", "B" or "C"
//turbulenceType: the turbulence type, can be "NTM", "ETM", "xEWM1" or "xEWM50" - where x is the turbine class (1,2 or 3)
//seed: the random seed for the turbulent windfield
//vertInf: vertical inflow angle in degrees [deg]
//horInf: horizontal inflow angle in degrees [deg]

 void  initializeSimulation();
//the simulation is reset and structural precomp is carried out

 void  setTimestepSize(double timestep);
//set the timestep size (in [s]), needs to be called befor initializeSimulation()

 void  setRPMPrescribeType_at_num(int type, int num);
//0 - RPM prescribed during ramp-up only, 1 - RPM prescribed for the whole simulation, 3 - no prescribed RPM, call before initializeSimulation()

 void  setRampupTime(double time);
//set the ramp-up time, call before initializeSimulation()

 void  getWindspeed(double x, double y, double z, double *velocity);
//gets the current windspeed at the chosen position (x,y,z), returns the windspeed vector in *velocity

 double  getCustomData_at_num(char *str, double pos, int num);
//get the current value of a custom channel (specify the data name as str) for the turbine instance i

 void  setInitialConditions_at_num(double yaw, double pitch, double azimuth, double rpm, int num);
//set the turbine initial yaw [deg], collective pitch [deg], azimuthal angle [deg] and initial rotSpeed [rpm], needs to be called before initializeSimulation()

 void  setTurbinePosition_at_num(double x, double y, double z, double rotx, double roty, double rotz, int num);
//set the turbine initial tower bottom x, y and z position [m], and xrot, yrot zrot rotation [deg], if other than (0,0,0,0,0,0) initially needs to be called befor initializeSimulation()
//can be used for cosimulation with a hydrodynamics software that models the floater

 void  setPowerLawWind(double windspeed, double horAngle, double vertAngle, double shearExponent, double referenceHeight);
//This function can be called at any time after the simulation has been initialized with initializeSimulation().
//This function defines a power law wind profile (https://en.wikipedia.org/wiki/Wind_profile_power_law) and the inflow direction.
//Using this function the windspeed can also be dynamically changed during a simulation.
//Brief description of the inputs:
//windspeed: constant windspeed in m/s [m/s]
//horAngle: the horizontal inflow angle in degrees [deg]
//vertAngle: the vertical inflow angle in degrees [deg]
//shearExponent: this is the exponent for the power law boundary layer profile, if this is set to 0 the windspeed is constant with height [-]
//referenceHeight: this is the height at which the velocity in the boundary layer is the defined windspeed, usually set to the hubheight [m]
//exemplary call: addTurbulentWind(12,115,115,220,20,60,0.1,"A","NTM",1000000,0,0);


 void  setControlVars_at_num(double *vars, int num);
//this sets the turbine control variables of the selected turbine (num) for torque and pitch
//if called after advanceController() and before advanceStructure(), the controller outputs can be overridden
//vars[0] = generator torque [Nm];
//vars[1] = yaw angle [deg];
//vars[2] = pitch blade 1 [deg];
//vars[3] = pitch blade 2 [deg];
//vars[4] = pitch blade 3 [deg];

 void  getTurbineOperation_at_num(double *vars, int num);
//this returns the turbines operational parameters of the selected turbine (num) to use as control input
//vars[0] = rotational speed [rad/s]
//vars[1] = power [W]
//vars[2] = HH wind velocity [m/s]
//vars[3] = yaw angle [deg]
//vars[4] = pitch blade 1 [deg]
//vars[5] = pitch blade 2 [deg]
//vars[6] = pitch blade 3 [deg]
//vars[7] = oop blade root bending moment blade 1 [Nm]
//vars[8] = oop blade root bending moment blade 2 [Nm]
//vars[9] = oop blade root bending moment blade 3 [Nm]
//vars[10] = ip blade root bending moment blade 1 [Nm]
//vars[11] = ip blade root bending moment blade 2 [Nm]
//vars[12] = ip blade root bending moment blade 3 [Nm]
//vars[13] = tor blade root bending moment blade 1 [Nm]
//vars[14] = tor blade root bending moment blade 2 [Nm]
//vars[15] = tor blade root bending moment blade 3 [Nm]
//vars[16] = oop tip deflection blade 1 [m]
//vars[17] = oop tip deflection blade 2 [m]
//vars[18] = oop tip deflection blade 3 [m]
//vars[19] = ip tip deflection blade 1 [m]
//vars[20] = ip tip deflection blade 2 [m]
//vars[21] = ip tip deflection blade 3 [m]
//vars[22] = tower top acceleration in global X [m/s^2]
//vars[23] = tower top acceleration in global Y [m/s^2]
//vars[24] = tower top acceleration in global Z [m/s^2]
//vars[25] = tower top fore aft acceleration [m/s^2]
//vars[26] = tower top side side acceleration [m/s^2]
//vars[27] = tower top X position [m]
//vars[28] = tower top Y position [m]
//vars[29] = tower bottom force along global X [Nm]
//vars[30] = tower bottom force along global Y [Nm]
//vars[31] = tower bottom force along global Z [Nm]
//vars[32] = tower bottom bending moment along global X [Nm]
//vars[33] = tower bottom bending moment along global Y [Nm]
//vars[34] = tower bottom bending moment along global Z [Nm]
//vars[35] = current time [s]
//vars[36] = azimuthal position of the LSS [deg]
//vars[37] = azimuthal position of the HSS [deg]
//vars[38] = HSS torque [Nm]
//vars[39] = wind speed at hub height [m/s]
//vars[40] = horizontal inflow angle [deg]

 void advanceController_at_num(double *vars, int num);
//this advancess the controller dll of the selected turbine (num) and controller outputs are applied to the turbine
//the controller ouputs are written to the vars array:
//vars[0] = generator torque [Nm];
//vars[1] = yaw angle [deg];
//vars[2] = pitch blade 1 [deg];
//vars[3] = pitch blade 2 [deg];
//vars[4] = pitch blade 3 [deg];

 void advanceTurbineSimulation();
//this advances the aerodynamic and structural model for all turbines and finishes the timestep
//the generator torque and pitch angles are applied during this step

 void runFullSimulation();
//this 

 void setLibraryPath(char *str);
// this tells where the dll is located
// this function is mandatory

 void getTowerBottomLoads_at_num(double *loads, int num);
 // this allows to export load at the bottom of the turbine tower



