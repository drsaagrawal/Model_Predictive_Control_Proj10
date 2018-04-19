
# Project Intro

This project is done in C++ using Model Predictive Control (MPC). Aim is to drive a simulated car around the track using waypoints. The simulator has latency of 100ms that has to be accounted. IPOPT and CPPAD libraries used to calculate an optimal trajectory and its associated actuation commands in order to minimize error with a third-degree polynomial fit to the given waypoints. 

# Rubric Ponits

### The Model: Student describes their model in detail. This includes the state, actuators and update equations.

The kinematic model includes the vehicle's x and y coordinates, orientation angle (psi), and velocity, as well as the cross-track error and psi error (epsi). Actuator outputs are acceleration and delta (steering angle).The model calculate current timestamp's state using previous actuation and state. Equations are:![eqns.png](attachment:eqns.png)

### Timestep Length and Elapsed Duration (N & dt): Student discusses the reasoning behind the chosen N (timestep length) and dt (elapsed duration between timesteps) values. Additionally the student details the previous values tried.

The values chosen for N and dt are 10 and 0.1 respectively. N is number of steps for the trajectory and dt is time elapsed between each steps. Other values tried are 25/0.1, 15/0.05 and many more.

### Polynomial Fitting and MPC Preprocessing: A polynomial is fitted to waypoints. If the student preprocesses waypoints, the vehicle state, and/or actuators prior to the MPC procedure it is described.

The waypoints are transformed into vehicle's perspective (main.cpp line 103-106). This simplifies the process to fit a polynomial to the waypoints because the vehicle's x and y coordinates are now at the origin (0, 0) and the orientation angle is also zero.

### Model Predictive Control with Latency: The student implements Model Predictive Control that handles a 100 millisecond latency. Student provides details on how they deal with latency.

I have taken care this in main.cpp as setting the dt with 0.1 and then prdicting the values taken into the latency. These values are then feeded in the state vector. This is in line 118-123 (main.cpp)
