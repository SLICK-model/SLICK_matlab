# The Stochastic Low-dimensional Inflated Convolutional Koopman model (SLICK)

The complex and chaotic dynamics, coupled with a wide range of length and time scales, remain an open challenge for the reduced-order modeling of turbulent flows.
Despite the understanding of the Navier-Stokes equations and the availability of accurate high-fidelity methods to solve them, real-time prediction of turbulent flows remains elusive across various fields, primarily due to their immense computational costs.
Available applications, therefore, have often been demonstrated using simplified dynamical systems, such as the Ginzburg-Landau equation and the Lorenz system, or canonical problems like cylinder wakes. This work aims to develop a reduced-order model (ROM) applicable to real-world data.

# Model-order reduction

![Alt text](Modal_analysis_schematic.png)

The choice of model order reduction usually goes hand-in-hand with the modal basis.
We categorize the approaches for identifying flow instabilities using three perspectives:
* Physical aspects (statistics or dynamics),
* State vector (linear $\to$ forced linear $\to$ nonlinear),
* Time-delay horizon (instantaneous $\to$ finite $\to$ long-time limit).
  
Note that all the methods here are applied to analyze general nonlinear fluid flows. The terminology 'linear,' 'forced linear,' and 'nonlinear' pertains to the state vector under analysis.

## Spectral proper orthogonal decomposition (SPOD)
SPOD identifies large coherent structures in stationary flows. Analogous to POD, SPOD is computed as the eigendecomposition of the cross-spectral density tensor and decomposes the flow into energy-ranked structures that evolve coherently in both space and time. Besides optimally accounting for the second-order statistics, SPOD modes are also dynamically significant as they are optimally averaged ensemble dynamical mode decomposition (DMD) modes and are formally equivalent to the spectral expansion of the stochastic Koopman operator. The SPOD expansion coefficients, $\mathbf{a}$, can be alternatively viewed as time-delay observables of the flow states, specifically through the Fourier convolution. This perspective can be understood from the inherent relationship between the Hankel singular vectors and SPOD modes. Based on these properties, our proposed model uses a finite-dimensional approximation of the linear Koopman operator to propagate the convolutional coordinates of SPOD in time.

# SLICK

![SLICK overview](Koopman_SPOD_schematic.png)


A graphical summary of the stochastic low-dimensional inflated convolutional Koopman (SLICK) model is shown. It has the following steps: 
* Collecting data as equally sampled snapshots,
* Performing model-order reduction and Fourier convolution,
* Training the SLICK model.

The turbulent open cavity flow is demonstrated as an example. 
While the model encompasses complex concepts, it can be formulated using a simple discrete-time algorithm.

# Stochastic complex Ginzburg-Landau equation

![SCGL overview](SCGL_nonlinear_correlated_full.png)

Here, we consider the SCGL equation, which is widely used as a model to study instabilities in spatially evolving flows.


