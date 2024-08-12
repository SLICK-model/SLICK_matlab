## The Stochastic Low-dimensional Inflated Convolutional Koopman model (SLICK)


## Model-order reduction

![Alt text](Modal_analysis_schematic.png)

The choice of model order reduction usually goes hand-in-hand with the modal basis.
We categorize the approaches for identifying flow instabilities using three perspectives:
* Physical aspects (statistics or dynamics),
* State vector (linear $\to$ forced linear $\to$ nonlinear),
* Time-delay horizon (instantaneous $\to$ finite $\to$ long-time limit).
  
Note that all the methods here are applied to analyze general nonlinear fluid flows. The terminology 'linear,' 'forced linear,' and 'nonlinear' pertains to the state vector under analysis.



## SLICK

![SLICK overview](Koopman_SPOD_schematic.png)


A graphical summary of the stochastic low-dimensional inflated convolutional Koopman (SLICK) model is shown. It has the following steps: 
* Collecting data as equally sampled snapshots,
* Performing model-order reduction and Fourier convolution,
* Training the SLICK model.

The turbulent open cavity flow is demonstrated as an example. 
While the model encompasses complex concepts, it can be formulated using a simple discrete-time algorithm.

## Stochastic complex Ginzburg-Landau equation

![SCGL overview](SCGL_nonlinear_correlated_full.png)

Here, we consider the SCGL equation, which is widely used as a model to study instabilities in spatially evolving flows.


