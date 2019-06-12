function [ ] = Run_36Cl_Facet()
% Function of the crep program that calculates exposure ages
    addpath(genpath('Functions'));
% Loading data
    % ice
    % load samples data
    Data = Load_data_36('INPUT/DATA_IN.xlsx');
    %load('Data_36_2.mat')
% Parameters given by user
    NumGMDB = 1; % 1: Mush; 2: GLOPIS; 3: LSD; 4: own user geomagnetic db
    Scheme = 3; % time dependant scaling model (1: LAL-STONE with cutoff rigidity, 2: LSD, 3: LAL-STONE 2000 no cutoff)
    Muon_model = 1; % 1: Exponential, 2: numeric integration of the flux (Balco 2008)
    Atm = 0; % 0: ERA40 (Uppala et al. 2005), 1: standard atmosphere equation (NOAA 1976)
    Data.GMDB = 1;
% !!!!!!!!!!!!!!!!! to be given by the user: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % Constants
     Data.lambda36 = 2.303e-6 ;
     Data.lambda36_uncert = Data.lambda36 .* 0.0066;

    % Attenuation length
     Data.Lambda_f_e = 160; % effective fast neutron attenuation coefficient (g.cm-2)
     Data.Lambda_f_t = 208; % TRUE fast neutron attenuation coefficient (g.cm-2)
     Data.Lambda_mu = 1510 ; % slow muon attenuation length (g.cm-2)
     Data.Psi_mu_0 = 190 ; % slow negative muon stopping rate at land surface (muon/g/an), Heisinger et al. (2002)

    % Unscaled sample specific 36Cl production rate by spallation of target elements
     Data.Psi_Cl36_Ca_0 = 42.2 ;% spallation production rate for Ca, SLHL (at of Cl36 /g of Ca per yr)
     Data.Psi_Cl36_K_0 = 148.1 ;% Spallation production rate at surface of 39K (at of Cl36 /g of Ca per yr)
     Data.Psi_Cl36_Ti_0 = 13 ; % Spallation production rate at surface of Ti (at of Cl36 /g of Ca per yr)
     Data.Psi_Cl36_Fe_0 = 1.9 ; % Spallation production rate at surface of Fe (at of Cl36 /g of Ca per yr)
    %uncertainties
     Data.Psi_Cl36_Ca_0_uncert = 4.8 ;% spallation production rate for Ca, SLHL (at of Cl36 /g of Ca per yr)
     Data.Psi_Cl36_K_0_uncert = 7.8 ;% Spallation production rate at surface of 39K (at of Cl36 /g of Ca per yr)
     Data.Psi_Cl36_Ti_0_uncert = 3 ; % Spallation production rate at surface of Ti (at of Cl36 /g of Ca per yr)
     Data.Psi_Cl36_Fe_0_uncert = 0.2 ; % Spallation production rate at surface of Fe (at of Cl36 /g of Ca per yr) 
% !!!!!!!!!!!!!!!!! to be given by the user: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    % Open paths and load cst
        addpath Functions
        addpath Constants
        load('Constants/GMDB.mat');
        load('Constants/OtherCst.mat');
 
        Data.Atm = Atm;% atmospheric model
        
    % Geomagnetic data base
    if length(Data.GMDB)==1;

        %NumGMDB = 2; % 1: Mush; 2: GLOPIS; 3: LSD; 4: own user geomagnetic db
        if NumGMDB==1;
            SelGMDB=GMDB.Musch;
        elseif NumGMDB==2;
            SelGMDB=GMDB.GLOPIS;
        else %  (NumGMDB=3)
            SelGMDB=GMDB.LSD;
        end
    else
        NumGMDB=4;
        SelGMDB=Data.GMDB;
    end

    % Get number of samples
    VecLat=Data.Lat';
    [NbSpl,~]=size(VecLat);

    % Formatting Outputs
    ExitMat=zeros(NbSpl,4);
    CellPDF=cell(1,2*NbSpl);
    ErrCol=zeros(NbSpl,1);
    StatCell=cell(NbSpl,1);

% Get 36Cl parameters 
    addpath(genpath('Functions/36Cl_Functions'));
    
            % Muon
            if(Muon_model == 1)
                flag.muon = 'exp';  % Muon attenuation length approximated by exponential (Schimmelfenning 2009).  
           elseif(Muon_model == 2)
                flag.muon = 'num'; % Muon attenuation length calculated following Heisinger (2002a,b) and Balco et al. (2008)
            else
                Mess=sprintf('Invlid muon model');
            end
            % Scaling model
            if(Scheme == 1)
                flag.scaling_model = 'st'; % LAL-STONE scheme cutoff rigidy
            elseif(Scheme == 2)
                flag.scaling_model = 'sa'; % LSD scheme
            elseif(Scheme == 3)
                flag.scaling_model = 'st2000'; % LAL-STONE scheme 2000 (no cuttoff)
            else
                Mess=sprintf('Invlid scaling model');
            end  
    % Constants initialization
        [Const_cosmo] = Init_const(Data,flag);   
        
    % Variable and Data Initialization 
        % rock
        [Data_formatted,Param_site] = Init_var(Data);

    % Scaling factors initialization
    %w = 0.2; % water content for Sato & Niita (2006)
    w = -1; % water content =  default value
            
    % geomagnetic database
    flag.NumGMDB = NumGMDB;
    Sf = Func_Sf(Param_site,Data,Atm,w,SelGMDB,flag);
    
    % Production rates and constants initialization
    Param_cosmo = clrock(Data_formatted,Param_site,Const_cosmo,Sf);
    
    %% Shielding factors

    Shf = Shielding_facet(Param_site,Const_cosmo);
    
    %% Lets calculate 36Cl for a given model

        % Search parameters
        flag.min_bounds = 0.0; % minimum bound for the search
        flag.max_bounds = 400000; % maximum bound for the search
        flag.plot = 0;
        flag.plotP = 0;
        
        % test a scenario
        Age_max = 300000; % maximum age (yr) to compute 36Cl produced in sample
        Age_PG = 32000; % Post-glacial age (yr) -> the sample is direclty at surface from this time to today.
        SR = 2.8; % rate (mm/yr) of sample uprising toward the surface (along a direction // to the fault plane)
        Z_samples = Data.Z(:); % depth of the samples (in cm)

        N_36 = Model_direct_36Facet(Age_max,Age_PG,SR,...
                                  Const_cosmo,Param_cosmo,Param_site,...
                                  Sf,Shf,Z_samples,flag)

        figure(3)
        plot(N_36,Data.Alt-Data.alt_scarp_top,'o'); hold on
        plot(Data.NuclCon,Data.Alt-Data.alt_scarp_top,'o')
        
     %% MCMc Inversion of the 36Cl concentration
     % erosion rate
     param_mcmc.search_bd(1,1)= 0.0 ; % minimum search bound (mm/yr)
     param_mcmc.search_bd(1,2)= 5.0 ; % maximum search bound (mm/yr)
     param_mcmc.search_std(1) = 1.0 ; % standart deviation for the proposal function (mm/yr)
     % Post-glacial age
     param_mcmc.search_bd(2,1)= 0.0 ; % minimum search bound (yr)
     param_mcmc.search_bd(2,2)= 30000 ; % maximum search bound (yr)
     param_mcmc.search_std(2) = 2000 ; % standart deviation for the proposal function (yr)
     % length of the chain
     param_mcmc.n_it = 10000;
     % Prepare dataset
     dataset(1,:) = Data.NuclCon(:);
     dataset(2,:) = Data.NuclErr(:);
     dataset(3,:) = Data.Z(:);
     Z_samples = Data.Z(:); % depth of the samples (in cm)
     tic
     [Age_mean,Age_std,SR_mean,SR_std] = inv_MCMC(Const_cosmo, Param_cosmo, Param_site, Sf, Shf, dataset, Z_samples, Age_max, param_mcmc, flag)
     toc