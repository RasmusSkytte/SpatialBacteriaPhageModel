% Scripts which generates all the figures

% Figure 2
!python3 functions/fig2.py
if ~exist('Fig_2','dir')
    mkdir('Fig_2')
end
!convert +append ../cpp/data/3D_Example_Full_Model/Visualize3D_py/Cells_density/007.png ../cpp/data/3D_Example_Full_Model/Visualize3D_py/Phages/007.png Fig_2/Fig2.tif

% Figure 3
run('functions/fig3')

!convert +append Fig_3/Dynamics_Model_1.tif Fig_3/Dynamics_Model_2.tif Fig_3/upper.tif
!convert +append Fig_3/Dynamics_Model_3.tif Fig_3/Dynamics_Model_4.tif Fig_3/lower.tif
!convert -append Fig_3/upper.tif Fig_3/lower.tif Fig_3/Fig3.tif
!rm Fig_3/Dynamics_Model_1.tif Fig_3/Dynamics_Model_2.tif Fig_3/Dynamics_Model_3.tif Fig_3/Dynamics_Model_4.tif
!rm Fig_3/upper.tif Fig_3/lower.tif

% Figure 4
run('functions/fig4')

% Figure 5
run('functions/fig5b')
run('functions/fig5c_and_figS10a')
warning('Figure 5 must be finished manually')

% Figure S1
run('functions/figS1')

% Figure S2 and S3
run('functions/figS2_and_figS3')

!convert +append Fig_S2/ShieldingFunction_S1_eta_1e+04.tif Fig_S2/ShieldingFunction_S2_eta_1e+04.tif Fig_S2/upper.tif
!convert +append Fig_S2/ShieldingFunction_S3_eta_1e+04.tif Fig_S2/ShieldingFunction_S4_eta_1e+04.tif Fig_S2/lower.tif
!convert -append Fig_S2/upper.tif Fig_S2/lower.tif Fig_S2/FigS2.tif
!rm Fig_S2/ShieldingFunction_S1_eta_1e+04.tif Fig_S2/ShieldingFunction_S2_eta_1e+04.tif Fig_S2/ShieldingFunction_S3_eta_1e+04.tif Fig_S2/ShieldingFunction_S4_eta_1e+04.tif
!rm Fig_S2/upper.tif Fig_S2/lower.tif

% Figure S4
run('functions/figS4')

% Figure S5
run('functions/figS5')

% Figure S6
run('functions/figS6')

% Figure S7
run('functions/figS7')

% Figure S8
run('functions/figS8')

% Figure S9
run('functions/figS9')

% Figure S10
run('functions/figS10b')

!convert -append Fig_S10/upper.tif Fig_S10/lower.tif Fig_S10/FigS10.tif
!rm Fig_S10/upper.tif Fig_S10/lower.tif

% Compress images
!find . -name '*.tif' -exec convert {} -compress lzw {} \;

% Cleanup
close all
