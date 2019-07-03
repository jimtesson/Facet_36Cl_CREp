function [ResultStat] = Inversion_36Cl_Facet()
% Function of the crep program that calculates exposure ages
    disp('Importing data')
%% Initialize path and functions 
    addpath(genpath('Functions'));
    addpath Functions
    addpath Constants
    addpath(genpath('Functions/36Cl_Functions'));
        
%% Loading data
    [Data,ParamUser] = Load_data_36('Input/DATA_IN.xlsx');
    
%% Set Atmospheric, Geomagnetic and scaling parameters
    ParamUser = SetUserChoices(Data,ParamUser);
    
%% Constants initialization
    Const_cosmo = Init_const(Data,ParamUser.flag);   
        
%% Formatting data
    [Data_formatted,Param_site] = Init_var(Data);

%% Compute scaling factors
    Sf = Func_Sf(Param_site,Data,ParamUser.Atm,ParamUser.w,ParamUser.SelGMDB,ParamUser.flag);
    
%% Compute production rates
    Param_cosmo = clrock(Data_formatted,Param_site,Const_cosmo,Sf);
    
%% Compute Shielding factors
    Shf = Shielding_facet(Param_site,Const_cosmo);
    
%% Prepare dataset
    dataset(1,:) = Data.NuclCon(:);
    dataset(2,:) = Data.NuclErr(:);
    dataset(3,:) = Data.Z(:);
    Z_samples = Data.Z(:); % depth of the samples (in cm)
    
%% Set plot parameters
    ParamUser.flag.plot = 0;
    ParamUser.flag.plotP = 0;
    
%% Packing variables
% Packing all required variables in a single cell for simplicity
    data_mc.Const_cosmo = Const_cosmo;
    data_mc.Param_cosmo = Param_cosmo;
    data_mc.Param_site = Param_site;
    data_mc.Sf = Sf;
    data_mc.Shf = Shf;
    data_mc.dataset = dataset;
    data_mc.Z_samples = Z_samples;
    data_mc.Age_max = ParamUser.Age_max;
    data_mc.flag = ParamUser.flag;
    data_mc.Cl36 = dataset(1,:);
    data_mc.Cl36_uncer = data_mc.dataset(2,:);
    data_mc.flag = ParamUser.flag;
    
%% Parameters of the MCMc Inversion
 
    % function providing the misfit
    ssfun = @(x0,data_mc) Get_misfit(x0,data_mc);
    
    % log-Likelihood
    logLike=@(m)log((2*pi)^-0.5.*sum(1./data_mc.Cl36_uncer)* exp(-ssfun(m,data_mc)/2));

    % Make an initial guess for the model parameters.
    m0=[ParamUser.SR_0 ParamUser.PG_age_0]';

    % Prior information
    logprior = @(m)(m(1)>ParamUser.SRmin)&(m(1)<ParamUser.SRmax)&(m(2)>ParamUser.Tmin)&(m(2)<ParamUser.Tmax);
    
%% test_forward: compute theoretical 36Cl for a given model

    if(ParamUser.test_forward == 1)
        fprintf( 1, '\t -> Test a model:\n \t Slip-rate = %f\n\tPost-glacial duration = %f\n',ParamUser.sr,ParamUser.tpg);
        % model 36Cl
        N_36 = Model_direct_36Facet([ParamUser.sr ParamUser.tpg],data_mc);
        % Figures: 36Cl profiles
        plot_36Cl_profile(N_36,Data,ParamUser);
        % print results
        fprintf( 1, '\t -> RMSw=%f \t Likelog = %f\t\n',ssfun(m0,data_mc),logLike(m0));
        save('Results/results_gwmcmc.mat')
    else
        
%% Initial sampling of models for inversion
    Ini_models=GenerateInitialSamples(ParamUser,logprior);

%% GWMCMc Inversion 

    disp('Starting Inversion')

    m=gwmcmc(Ini_models,{logprior logLike},ParamUser.n_models_inversion,'burnin',0,'Parallel',ParamUser.parallel_computing);
    
    % save all
    save('Results/results_gwmcmc.mat')
    
    %% plot
    [ResultStat] = Plot_results_inversion();
    
    end
end

function plot_36Cl_profile(N_36_site,Data,ParamUser)
% plot the 36Cl profile of a model

        figure
        % plot modelled 36Cl concentrations
        plot(N_36_site,Data.Alt,'o'); hold on
        % plot measured 36Cl concentrations
            % error bar coordinates on measurments
            X = [Data.NuclCon-Data.NuclErr Data.NuclCon+Data.NuclErr]' ;
            Y = [Data.Alt(:) Data.Alt(:)]' ;
        plot(Data.NuclCon,Data.Alt,'ko','MarkerFaceColor','black')
        plot(X,Y,'-k')
        ax_x = xlim; xlim([0 ax_x(2)])
        title(sprintf('Slip-rate = %2.1f mm/yr, Post-glacial duration = %2.0f kyr',ParamUser.sr,ParamUser.tpg/1000))
        xlabel('36Cl concentration (at/gr)') 
        ylabel('Altitude of the sample') 
        legend('model','data','Location','southeast')
        saveas(gcf,['Results/36Cl_forwardmodel_' num2str(get(gcf,'Number')) '.fig'])

end


function ball=GenerateInitialSamples(ParamUser,logprior)
% Function to generate the initial samples of models for inversion
    ball=zeros(2,ParamUser.n_walker); 
    tmp=zeros(2,1);
    disp('Initializing Inversion')
    disp('First we initialize the ensemble of walkers')
    disp(['-> number of walkers: ' num2str(ParamUser.n_walker)])
    disp('-> guess models: ')
    i=0;
    while i<ParamUser.n_walker
        tmp(1)=ParamUser.SRmin+rand(1,1).*(ParamUser.SRmax-ParamUser.SRmin);
        tmp(2)=ParamUser.Tmin+rand(1,1).*(ParamUser.Tmax-ParamUser.Tmin);
        if(logprior(tmp)==1)
            i=i+1;
            ball(1,i)=tmp(1);
            ball(2,i)=tmp(2); 
            fprintf ( 1, 'model:%i -> Slip-rate = %4.2f, Post-glacial duration = %6.0f\n',i,ball(1,i),ball(2,i));

        end
    end
end


function ParamUser_out = SetUserChoices(Data,ParamUser_in)

    disp('Computing parameters')
    
    % Load constants
    load('Constants/GMDB.mat');
    load('Constants/OtherCst.mat');
    
    % Geomagnetic data base
    ParamUser_in.flag.NumGMDB = ParamUser_in.NumGMDB;
    
    if length(Data.GMDB)==1
        %ParamUser.NumGMDB = 2; % 1: Mush; 2: GLOPIS; 3: LSD; 4: own user geomagnetic db
        if ParamUser_in.NumGMDB==1
            ParamUser_in.SelGMDB=GMDB.Musch;
        elseif ParamUser_in.NumGMDB==2
            ParamUser_in.SelGMDB=GMDB.GLOPIS;
        elseif ParamUser_in.NumGMDB==3
            ParamUser_in.SelGMDB=GMDB.LSD;
        else
            load('Constants/GMDB_own.mat');
            ParamUser_in.SelGMDB=GMDB_own.GMDB;
        end
    end

%% Set the muon model
    if(ParamUser_in.Muon_model == 1)
        ParamUser_in.flag.muon = 'exp';  % Muon attenuation length approximated by exponential (Schimmelfenning 2009).  
   elseif(ParamUser_in.Muon_model == 2)
        ParamUser_in.flag.muon = 'num'; % Muon attenuation length calculated following Heisinger (2002a,b) and Balco et al. (2008)
    else
        fprintf('Invlid muon model');
    end
%% Set the Scaling model
    if(ParamUser_in.Scheme == 1)
        ParamUser_in.flag.scaling_model = 'st'; % LAL-STONE ParamUser.Scheme cutoff rigidy
    elseif(ParamUser_in.Scheme == 2)
        ParamUser_in.flag.scaling_model = 'sa'; % LSD ParamUser.Scheme
    elseif(ParamUser_in.Scheme == 3)
        ParamUser_in.flag.scaling_model = 'st2000'; % LAL-STONE ParamUser.Scheme 2000 (no cuttoff)
    else
        fprintf('Invlid scaling model');
    end  

%% Parameters for LSD scaling model    
    %w = 0.2; % water content for Sato & Niita (2006)
    ParamUser_in.w = -1; % water content =  default value
    
    ParamUser_out = ParamUser_in;
end