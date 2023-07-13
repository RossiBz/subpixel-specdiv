%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Project: Uncovering the hidden: Leveraging sub-pixel spectral diversity to estimate plant diversity from space

% Script Purpose: Estimate endmember diversity for simulated communities, as done in Rossi and Gholizadeh (2023) Experiment 1
%# using SISAL and MVC-NMF as endmember extraction techniques

% Date: 31.01.2023

% Author: Christian Rossi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.


% define path of simulated hyperspectral data
data_path = 'Q:/prjdata/botany/grassland/tallgrass_prairie_oklahoma/code/output/';

%define file name
FILENAME_HS ='simulated_communities_field_red100_soil2_3speciesup03_SNR60.csv';

comm_data = readmatrix([data_path FILENAME_HS]);


% last 17 rows are the endmembers  
comm_data_spec =comm_data(1:15300,6:end)';

[B N]=size(comm_data_spec);

comm_data_cube = reshape(comm_data_spec',[153 100 2001]); 

% estimate number of endmembers with the noise-whitened Harsanyi–Farrand–Chang method -----------------------

numEndmembers=countEndmembersHFC(comm_data_cube);

% extract endmembers-------------------------------------------------------------------------

extractMethod ='sisal' % chose method sisal or mvcnmf

for i = 1:10
    
    switch extractMethod
        
        case  'mvcnmf'
            
            x_unit=comm_data_spec;
            [UU, SS, VV] = svds(x_unit,numEndmembers);
            
            % PCA
            PrinComp= pca(x_unit');
            meanData = mean(x_unit');
            
            % test mvcnmf
            tol = 1e-6;
            maxiter = 150;
            T = 0.015;
            showflag = 1;
            
            % use vca to intalize
            A_vca = vca(x_unit,'Endmembers', numEndmembers,'verbose','on');
            
            % FCLS
            warning off;
            AA = [1e-5*A_vca;ones(1,length(A_vca(1,:)))];
            s_fcls = zeros(length(A_vca(1,:)),N);
            for j=1:N
                r = [1e-5*x_unit(:,j); 1];
                %   s_fcls(:,j) = nnls(AA,r);
                s_fcls(:,j) = lsqnonneg(AA,r);
            end
            
            Ainit = A_vca;
            sinit = s_fcls;
            
            % use conjugate gradient to find A can speed up the learning
            [endmembersHSI, sest] = mvcnmf(x_unit,Ainit,sinit,UU,PrinComp,meanData,T,tol,maxiter,showflag,2,1);
            
            endmembersHSI=endmembersHSI;
            
            
        case 'sisal'
            
            %   The application of this projection ensures that the data is in
            %   an affine set.
            %
            %   Up is an isometric matrix that spans the subspace where Y lives
            [Y,Up,my,sing_val] = dataProj(comm_data_spec,numEndmembers,'proj_type','affine');
            
            [endmembersHSI] = sisal(Y,numEndmembers, 'spherize', 'yes','MM_ITERS',100, 'TAU',10, 'verbose',2);
            
    end
    
    % estimate endmember abundance per community------------------------------------------------
	
    abundance_comm = estimateAbundanceLS(comm_data_cube,endmembersHSI,'Method','fcls');
    
    listOfValues = reshape(abundance_comm, 153*100, numEndmembers);
    
	% save endmemember abundance 
    writematrix(horzcat(comm_data(1:15300,1:5),listOfValues),[extractMethod '_simu_communities_SNR60_soil2_red100_field' num2str(i) '.csv']);
end