# tskf_estimator_pkg
Estimator based on Two Stage Kalman Filter [1]

### Inputs
- Commanded vertical thrust and angular velocities
- Drone's linear positions and angular velocities 


### Outputs
- /fault_estim: estimation of the fault value (no very accurate)
- /f_detection: boolen vector which indicates the faulty motors

### Tested
- Configuration files for Hummingbird (+ configuration) and Iris (x configuration)

## References
<a id="1">[1]</a> 
M. Hadi Amoozgar, A.Chamseddine,Y.Zhang, (2013).
Experimental Test of a Two-Stage Kalman Filter for Actuator Fault Detection and Diagnosis of an Unmanned Quadrotor Helicopter.
J Intell Robot Syst, 70:107–117.
