# xbohos-A-Novel-adapative-matching-pursuit-algorithm-for-moving-force-identification-
Ref: B. H. Xu and L. Yu, A novel regularized adaptive matching pursuit for moving force identification using multiple criteria and prior knowledge, Int. J. Struct. Stab. Dyn. 23(10) (2023) 2350117.
If you use this code, please cite 'Xu B H, Yu L. A novel regularized adaptive matching pursuit for moving force identification using multiple criteria and prior knowledge[J]. International Journal of Structural Stability and Dynamics
, 2023, 23(10): 2350117.' If you have any questions, you can contact with me through xbohos@163.com.

Note : huisuRAMP.m and selection1.m are the main function for NRAMP algorithm. They are contained in single impluse force (main0515.m) and two moving forces (main0514_1.m) simulations. There's no need for further enhancement of the
function to simulate additional axles. Users can improve or modify it as needed since this algorithm is easy to understand! 😊

ydhzdouble_force.m : the function to obtain the force vecoter.

CS_OMP.m : the function of orthogonal matching pursuit algorithm (OMP).

CS_ROMP.m : the function of regularized OMP algorithm. 'Regularize.m' is its sub-function. If you use this code, please cite 'D. Needell and R. Vershynin, Uniform uncertainty principle and signal recovery via regularized orthogonal matching pursuit, Found. Comput. Math. 9(3) (2009) 317–334.'

main0515.m: an example for single impluse force under different noise levels for the proposed method. It should be noted that the 'error', 'f1' and 'f2' obtained during the calculation process contain the corresponding values of the OMP and ROMP algorithms.

main0514_1.m : an example for two moving forces under different response combinations in the case of 10% Noise level. The calculation results are saved in 'error_220512.mat', 'f_220512.mat' and 't_220512.mat'. User can use 'figure_plot_f.m' function to further obtain the figure of comparative studies of OMP, ROMP and NRAMP methods.
