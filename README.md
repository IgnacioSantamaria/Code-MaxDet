# Code-MaxDet

This code package contains a simulation environment, based on Matlab, that reproduces the numerical results and figures in the paper

Ignacio Santamaria, Mohammad Soleymani, Jesus Gutierrez, and  Eduard Jorswieck, and "Optimal symmetric low-rank BD-RIS configuration maximizing the determinant of a MIMO link", submitted to IEEE Signal Processing, 2026.

Abstract.- Beyond-diagonal reconfigurable intelligent surfaces (BD-RISs) significantly improve wireless performance by allowing tunable interconnections among elements, but their design in multiple-input multiple-output (MIMO) systems has so far relied on complex iterative algorithms or suboptimal approximations.  This work introduces a simple yet powerful approach: instead of directly maximizing the achievable rate, we maximize the absolute value of the determinant of the equivalent MIMO channel. We derive a closed-form symmetric unitary scattering matrix whose rank is exactly twice the channel’s degrees of freedom \(2r\). Remarkably, this low-rank solution achieves the same determinant value as the optimal unitary BD-RIS. Using log-majorization theory, we prove that the rate loss relative to the optimal unitary BD-RIS vanishes at high signal-to-noise ratio (SNR) or when the number of BD-RIS elements becomes large. Moreover, the proposed solution can be perfectly implemented using a $q$-stem BD-RIS architecture with only \(q = 2r-1\) stems, requiring a minimum number of reconfigurable circuits. The resulting Max-Det solution is orders of magnitude faster to compute than existing iterative methods while achieving near-optimal rates in practical scenarios. This makes high-performance BD-RIS deployment feasible even with large surfaces and limited computational resources.  

Acknowledgements: This work is supported by the European Commission’s Horizon Europe, Smart Networks and Services Joint Undertaking, research and innovation program under grant agreement 101139282, 6G-SENSES project. The work of I. Santamaria was also partly supported under grant PID2022-137099NB-C43 (MADDIE) funded by MICIU/AEI /10.13039/501100011033 and FEDER, UE.


