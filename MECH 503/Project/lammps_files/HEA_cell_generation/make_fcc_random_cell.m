%% This program generates a fcc random solid solution hea with Ne elements
% The code will generate 10 different random cells for an equiatomic HEA
% containing 4 elements. You can read these files with the script provided
% in ELASTIC folder and modify the cell generation command by read_data as
% shown below (line 26 - 34). Potential file usage is also illustrated below. 
% This code is part of a much larger project hosted in github: 
% https://github.com/Mponga/HEA_vacancy_mcml 

clear all, clc, close all 

[Pars] = input_cnn();

[xcoord,ycoord,zcoord,Pars.Na]=fcc_cluster(Pars);

% Make several mixtures to assess the chem environment
generations = make_random_generations(Pars);

for i=1:Pars.Ngen

    best_gen=i;

    % Print LAMMPS file to read with read_data command.
    % to read these file use:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ﻿# --------------------- ATOM DEFINITION ------------------------
    %
    % read_data hea-genR1_5x5x5_HEA.lmp
    %
    % # ------------------------ FORCE FIELDS -----------------------
    % pair_style      eam/alloy
    % pair_coeff  * * Fe_Ni_Cr_Co_Al.setfl  Ni Fe Cr Co
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    formatSpec = 'hea-genR%d_%dx%dx%d_HEA.lmp';
    filename = sprintf(formatSpec,best_gen,Pars.Nx,Pars.Ny,Pars.Nz);

    print_lammps_file([xcoord ycoord zcoord], Pars, best_gen, generations, filename);

    % plot the atoms to vizualize them only
    figure(i)
    scatter3(xcoord,ycoord,zcoord,200,generations(i,:),'filled')

end

