function [ Data,Param ] = Load_data_36(file_name)
%Function in order to get data from the xsl file

% import data

    % import data
    [~,tmp1,raw] = xlsread(file_name,1);
    % number of Samples
    [ns,~]= size(tmp1);
    ns= ns-1;
    nrowmax = ns+1;

% samples names
    Data.Samples = raw(2:nrowmax,1);

% Location
    Data.Lat = cell2mat(raw(2:nrowmax,2));
    Data.Lon = cell2mat(raw(2:nrowmax,3));
    Data.Alt = cell2mat(raw(2:nrowmax,4));

% Depth. Must be = .0 if taken at surface
    Data.Z =  cell2mat(raw(2:nrowmax,5));

% [36Cl]
    Data.NuclCon = cell2mat(raw(2:nrowmax,6));
    Data.NuclErr = cell2mat(raw(2:nrowmax,7));

% Shielding
    Data.Shield = 1;%cell2mat(raw(2:nrowmax,8));

% density
    Data.Dens = cell2mat(raw(2:nrowmax,8));

% Thickness
    Data.Thick = cell2mat(raw(2:nrowmax,9));

% Erosion rate (cm/yr)
    Data.Eros = 0.0;%cell2mat(raw(2:nrowmax,11));


% Age of formation
    Data.Age_form = cell2mat(raw(2:nrowmax,10));
    Data.Age_form_er = cell2mat(raw(2:nrowmax,11));

%% Target Chemistry
    Data.Cl_targ = cell2mat(raw(2:nrowmax,12)); % (ppm)
    Data.Cl_targ_er = cell2mat(raw(2:nrowmax,13)); % (ppm)
    Data.SiO2_targ = cell2mat(raw(2:nrowmax,14)); % (in %)
    Data.Al2O3_targ = cell2mat(raw(2:nrowmax,15)); % (in %)
    Data.Fe2O3_targ = cell2mat(raw(2:nrowmax,16)); % (in %)
    Data.Fe2O3_targ_er = cell2mat(raw(2:nrowmax,17)); % (in %)
    Data.MnO_targ = cell2mat(raw(2:nrowmax,18)); % (in %)
    Data.MgO_targ = cell2mat(raw(2:nrowmax,19)); % (in %)
    Data.CaO_targ = cell2mat(raw(2:nrowmax,20)); % (in %)
    Data.CaO_targ_er = cell2mat(raw(2:nrowmax,21)); % (in %)
    Data.Na2O_targ = cell2mat(raw(2:nrowmax,22)); % (in %)
    Data.K2O_targ = cell2mat(raw(2:nrowmax,23)); % (in %)
    Data.K2O_targ_er = cell2mat(raw(2:nrowmax,24)); % (in %)
    Data.TiO2_targ = cell2mat(raw(2:nrowmax,25)); % (in %)
    Data.TiO2_targ_er = cell2mat(raw(2:nrowmax,26)); % (in %)
    Data.P2O5_targ = cell2mat(raw(2:nrowmax,27)); % (in %)
    Data.LOI_targ = cell2mat(raw(2:nrowmax,28)); % (in %)
    Data.H2O_targ = cell2mat(raw(2:nrowmax,29)); % (in %)
    Data.CO2_targ = cell2mat(raw(2:nrowmax,30)); % (in %)


%% Bulk Chemistry
    Data.SiO2_bulk = cell2mat(raw(2:nrowmax,31)); % (in %)
    Data.Al2O3_bulk = cell2mat(raw(2:nrowmax,32)); % (in %)
    Data.Fe2O3_bulk = cell2mat(raw(2:nrowmax,33)); % (in %)
    Data.MnO_bulk = cell2mat(raw(2:nrowmax,34)); % (in %)
    Data.MgO_bulk = cell2mat(raw(2:nrowmax,35)); % (in %)
    Data.CaO_bulk = cell2mat(raw(2:nrowmax,36)); % (in %)
    Data.NaO_bulk = cell2mat(raw(2:nrowmax,37)); % (in %)
    Data.K2O_bulk = cell2mat(raw(2:nrowmax,38)); % (in %)
    Data.TiO2_bulk = cell2mat(raw(2:nrowmax,39)); % (in %)
    Data.P2O5_bulk = cell2mat(raw(2:nrowmax,40)); % (in %)
    Data.H2O_bulk = cell2mat(raw(2:nrowmax,41)); % (in %)
    Data.CO2_bulk = cell2mat(raw(2:nrowmax,42)); % (in %)
    Data.S_bulk = cell2mat(raw(2:nrowmax,43)); % (in %)
    Data.LOI_bulk = cell2mat(raw(2:nrowmax,44)); % (in %)
    Data.Li_bulk = cell2mat(raw(2:nrowmax,45)); % (in %)
    Data.B_bulk = cell2mat(raw(2:nrowmax,46)); % (in %)
    Data.Cl_bulk = cell2mat(raw(2:nrowmax,47)); % (in %)
    %Data.Cr_bulk = cell2mat(raw(2:nrowmax,50)); % (in %)
    %Data.Co_bulk = cell2mat(raw(2:nrowmax,51)); % (in %)
    Data.Sm_bulk = cell2mat(raw(2:nrowmax,48)); % (in %)
    Data.Gd_bulk = cell2mat(raw(2:nrowmax,49)); % (in %)
    Data.Th_bulk = cell2mat(raw(2:nrowmax,50)); % (in %)
    %Data.Th_bulk_er = cell2mat(raw(2:nrowmax,55)); % (in %)
    Data.U_bulk = cell2mat(raw(2:nrowmax,51)); % (in %)
    %Data.U_bulk_er = cell2mat(raw(2:nrowmax,57)); % (in %)

%% Site parameters
    % import data
    [~,~,raw] = xlsread(file_name,2);

    Data.alt_scarp_top = cell2mat(raw(2,2)); 
    Data.alpha = cell2mat(raw(3,2)); 
    Data.beta = cell2mat(raw(4,2)); 
    Data.gamma = cell2mat(raw(5,2)); 
    Data.rho_coll = cell2mat(raw(6,2)); 

%% Inversion parameters
    % 36Cl
    Data.lambda36 = 2.303e-6 ;
    Data.lambda36_uncert = Data.lambda36 .* 0.0066;
    
       % Attenuation length
    Data.Lambda_f_e = cell2mat(raw(9,2));  % effective fast neutron attenuation coefficient (g.cm-2)
    Data.Lambda_f_t = cell2mat(raw(10,2));  % TRUE fast neutron attenuation coefficient (g.cm-2)
    Data.Lambda_mu = cell2mat(raw(11,2));  % slow muon attenuation length (g.cm-2)
    Data.Psi_mu_0 = cell2mat(raw(12,2));  % slow negative muon stopping rate at land surface (muon/g/an), Heisinger et al. (2002)

    % Unscaled sample specific 36Cl production rate by spallation of target elements
    Data.Psi_Cl36_Ca_0 = cell2mat(raw(13,2));  % spallation production rate for Ca, SLHL (at of Cl36 /g of Ca per yr)
    Data.Psi_Cl36_K_0 = cell2mat(raw(14,2)); % Spallation production rate at surface of 39K (at of Cl36 /g of Ca per yr)
    Data.Psi_Cl36_Ti_0 = cell2mat(raw(15,2));  % Spallation production rate at surface of Ti (at of Cl36 /g of Ca per yr)
    Data.Psi_Cl36_Fe_0 = cell2mat(raw(16,2));  % Spallation production rate at surface of Fe (at of Cl36 /g of Ca per yr)
    %uncertainties
    Data.Psi_Cl36_Ca_0_uncert = cell2mat(raw(17,2)); % spallation production rate for Ca, SLHL (at of Cl36 /g of Ca per yr)
    Data.Psi_Cl36_K_0_uncert = cell2mat(raw(18,2)); % Spallation production rate at surface of 39K (at of Cl36 /g of Ca per yr)
    Data.Psi_Cl36_Ti_0_uncert = cell2mat(raw(19,2));  % Spallation production rate at surface of Ti (at of Cl36 /g of Ca per yr)
    Data.Psi_Cl36_Fe_0_uncert = cell2mat(raw(20,2));  % Spallation production rate at surface of Fe (at of Cl36 /g of Ca per yr) 

    
    Param.NumGMDB = cell2mat(raw(21,2)); % 1: Mush; 2: GLOPIS; 3: LSD; 4: own user geomagnetic db
    Param.Scheme = cell2mat(raw(22,2)); % time dependant scaling model (1: LAL-STONE with cutoff rigidity, 2: LSD, 3: LAL-STONE 2000 no cutoff)
    Param.Muon_model = 1; % 1: Exponential, 2: numeric integration of the flux (Balco 2008)
    Param.Atm = cell2mat(raw(23,2)); % 0: ERA40 (Uppala et al. 2005), 1: standard atmosphere equation (NOAA 1976)
    Data.GMDB = 1;
    
    Param.Age_max = cell2mat(raw(26,2)); % maximum age (yr) to compute 36Cl produced in sample
    
    Param.PG_age_0 = cell2mat(raw(27,2)); % initial guess for post-glacial age (yr)
    Param.SR_0 = cell2mat(raw(28,2)); % initial guess for the fault slip-rate (mm/yr)
    
    Param.SRmin = cell2mat(raw(29,2)); % Inversion minimum slip-rate bound (mm/yr)
    Param.SRmax = cell2mat(raw(30,2)); % Inversion maximum slip-rate bound (mm/yr)
    Param.SR_std = cell2mat(raw(31,2));
    
    Param.Tmin = cell2mat(raw(32,2)); % Inversion minimum post-glacial age bound (yr)
    Param.Tmax = cell2mat(raw(33,2)); % Inversion maximim post-glacial age bound (yr)
    Param.T_std = cell2mat(raw(34,2));
    
    Param.n_walker = cell2mat(raw(35,2)); % number of chain
    Param.n_models_inversion = cell2mat(raw(36,2)); % number of models generated during the inversion
    Param.parallel_computing = logical(cell2mat(raw(37,2))); % inversion is operated using parallel computing (true, false)
    
    Param.BurnIn = cell2mat(raw(38,2)); % ( >=0 and <1) Proportion of the chain removed, must be >= 0 and <100
    Param.N_stat = cell2mat(raw(39,2)); % number of samples randomly picked to draw 36Cl concentrations from the inversion chain
    
    Param.test_forward = cell2mat(raw(42,2)); % test a forward model ?
    % input parameters to test a forward model
    Param.sr = cell2mat(raw(43,2)); % mm/yr
    Param.tpg = cell2mat(raw(44,2)); % yr

end

