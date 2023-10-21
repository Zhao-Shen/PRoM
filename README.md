## Venc Design and Velocity Estimation for Phase Contrast MRI
* Shen Zhao, The Ohio State University (zhao.1758@osu.edu, shenzhao@stanford.edu)
* Rizwan Ahmad, The Ohio State University (ahmad.46@osu.edu)
* Lee C. Potter, The Ohio State University (potter.36@osu.edu)

In phase-contrast magnetic resonance imaging (PC-MRI), spin velocity contributes to the phase measured at each voxel. Therefore, estimating velocity from potentially wrapped phase measurements is the task of solving a system of noisy congruence equations. We propose <em> Phase Recovery from Multiple Wrapped Measurements </em> (<strong>PRoM</strong>) as a fast, approximate maximum likelihood estimator of velocity from multi-coil data with possible amplitude attenuation due to dephasing. The estimator can recover the fullest possible extent of unambiguous velocities, which can greatly exceed twice the highest venc. The estimator uses all pairwise phase differences and the inherent correlations among them to minimize the estimation error. Correlations are directly estimated from multi-coil data without requiring knowledge of coil sensitivity maps, dephasing factors, or the actual per-voxel signal-to-noise ratio. Derivation of the estimator yields explicit probabilities of unwrapping errors and the probability distribution for the velocity estimate; this, in turn, allows for optimized design of the phase-encoded acquisition. These probabilities are also incorporated into spatial post-processing to further mitigate wrapping errors. Simulation, phantom, and in vivo results for three-point PC-MRI acquisitions validate the benefits of reduced estimation error, increased recovered velocity range, optimized acquisition, and fast computation. 

![Figure_INVIVO](https://user-images.githubusercontent.com/62859186/173437311-65f93e02-4876-4c09-a0cc-238e0e46268a.png)

## Implementation
* To implement and generate the figure of in vivo processing results, run Figure_In_Vivo.m

## Extension to General Linear Mixed Integer Least Squares
https://github.com/tuckerda/Mixed-Integer-Linear-Models-MLE

## References
1. Zhao, Shen, Rizwan Ahmad, and Lee C. Potter. "Venc Design and Velocity Estimation for Phase Contrast MRI." arXiv preprint arXiv:2109.12481 (2021). https://arxiv.org/abs/2109.12481
2. Zhao, Shen, Rizwan Ahmad, and Lee C. Potter. "Maximizing Unambiguous Velocity Range in Phase-Contrast MRI with Multipoint Encoding." 2022 IEEE 19th International Symposium on Biomedical Imaging (ISBI). IEEE, 2022. https://ieeexplore.ieee.org/document/9761589

