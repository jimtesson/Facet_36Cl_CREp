function [ResultStat] = Inversion_36Cl_Facet()
% Function of the crep program that calculates exposure ages


    disp('Importing data')
%% Initialize path and functions 
    addpath(genpath('Functions'));
    addpath Functions
    addpath Constants
    addpath(genpath('Functions/36Cl_Functions'));
    
%% Load constants
    load('Constants/GMDB.mat');
    load('Constants/OtherCst.mat');
    
%% Loading data
    [Data,ParamUser] = Load_data_36('Input/DATA_IN.xlsx'); %

    % number of site 
    n_site = Data{1}.n_site;
    
%% Set ParamUser.Atmospheric and Geomagnetic parameters

    disp('Computing parameters')
        
    % Geomagnetic data base
    if length(Data{1}.GMDB)==1
        %ParamUser.NumGMDB = 2; % 1: Mush; 2: GLOPIS; 3: LSD; 4: own user geomagnetic db
        if ParamUser{1}.NumGMDB==1
            SelGMDB=GMDB.Musch;
        elseif ParamUser{1}.NumGMDB==2
            SelGMDB=GMDB.GLOPIS;
        else %  (ParamUser.NumGMDB=3)
            SelGMDB=GMDB.LSD;
        end
    else
        ParamUser{1}.NumGMDB=4;
        SelGMDB=Data{1}.GMDB;
    end

%% Set the muon model
    if(ParamUser{1}.Muon_model == 1)
        flag.muon = 'exp';  % Muon attenuation length approximated by exponential (Schimmelfenning 2009).  
   elseif(ParamUser{1}.Muon_model == 2)
        flag.muon = 'num'; % Muon attenuation length calculated following Heisinger (2002a,b) and Balco et al. (2008)
    else
        Mess=sprintf('Invlid muon model');
    end
%% Set the Scaling model
    if(ParamUser{1}.Scheme == 1)
        flag.scaling_model = 'st'; % LAL-STONE ParamUser.Scheme cutoff rigidy
    elseif(ParamUser{1}.Scheme == 2)
        flag.scaling_model = 'sa'; % LSD ParamUser.Scheme
    elseif(ParamUser{1}.Scheme == 3)
        flag.scaling_model = 'st2000'; % LAL-STONE ParamUser.Scheme 2000 (no cuttoff)
    else
        Mess=sprintf('Invlid scaling model');
    end  
    
    Const_cosmo = cell(1,n_site);
    Data_formatted = cell(1,n_site);
    Param_site = cell(1,n_site);
    Sf = cell(1,n_site);
    Param_cosmo = cell(1,n_site);
    Shf = cell(1,n_site);
    dataset = cell(1,n_site);
    Z_samples = cell(1,n_site);
    data_mc = cell(1,n_site);
    
   for i = 1:n_site
       
    
        %% Constants initialization
        %eval(sprintf('[Const_cosmo_%d] = Init_const(Data_%d,flag);',i,i));
        Const_cosmo{i} = Init_const(Data{i},flag);
        
        %% Formatting data
        %eval(sprintf('[Data_formatted_%d,Param_site_%d] = Init_var(Data_%d);',i,i,i));
        [Data_formatted{i},Param_site{i}] = Init_var(Data{i});
        
        %% Compute scaling factors
        %w = 0.2; % water content for Sato & Niita (2006)
        w = -1; % water content =  default value
        flag.NumGMDB = ParamUser{1}.NumGMDB;
        %eval(sprintf('Sf_%d = Func_Sf(Param_site_%d,Data_%d,ParamUser_%d.Atm,w,SelGMDB,flag);',i,i,i,i));
        Sf{i} = Func_Sf(Param_site{i},Data{i},ParamUser{i}.Atm,w,SelGMDB,flag);
        
        %% Compute production rates
        %eval(sprintf('Param_cosmo_%d = clrock(Data_formatted_%d,Param_site_%d,Const_cosmo_%d,Sf_%d);',i,i,i,i,i));
        Param_cosmo{i} = clrock(Data_formatted{i},Param_site{i},Const_cosmo{i},Sf{i});
        
        %% Compute Shielding factors
        %eval(sprintf('Shf_%d = Shielding_facet(Param_site_%d,Const_cosmo_%d);',i,i,i));
        Shf{i} = Shielding_facet(Param_site{i},Const_cosmo{i});
        
        %% Prepare dataset
%         eval(sprintf('dataset_%d(1,:) = Data_%d.NuclCon(:);',i,i));
%         eval(sprintf('dataset_%d(2,:) = Data_%d.NuclErr(:);',i,i));
%         eval(sprintf('dataset_%d(3,:) = Data_%d.Z(:);',i,i));
%         eval(sprintf('Z_samples_%d = Data_%d.Z(:);',i,i)); % depth of the samples (in cm)
        
        n_sample=length(Data{i}.NuclCon);
        
        dataset{i}=zeros(3,n_sample);
        dataset{i}(1,:) = Data{i}.NuclCon(:);
        dataset{i}(2,:) = Data{i}.NuclErr(:);
        dataset{i}(3,:) = Data{i}.Z(:);
        
        Z_samples{i} = Data{i}.Z(:);
        
        %% Set plot parameters
        flag.plot = 0;
        flag.plotP = 0;
        
        %% Packing variables
        % Packing all required variables in a single cell for simplicity
        
        data_mc{i}.Const_cosmo = Const_cosmo{i};
        data_mc{i}.Param_cosmo = Param_cosmo{i};
        data_mc{i}.Param_site = Param_site{i};
        data_mc{i}.Sf = Sf{i};
        data_mc{i}.Shf = Shf{i};
        data_mc{i}.dataset = dataset{i};
        data_mc{i}.Z_samples = Z_samples{i};
        data_mc{i}.Age_max = ParamUser{i}.Age_max;
        data_mc{i}.flag = flag;
        data_mc{i}.Cl36 = dataset{i}(1,:);
        data_mc{i}.Cl36_uncer = dataset{i}(2,:);


    end
    
    
%% Parameters of the MCMc Inversion

    % function providing the misfit for a single site
    ssfun = @(x0,thedata_mc) Get_misfit(x0,thedata_mc);
    
    % get the vector of all 36Cl uncertainties
    std_36Cl_er = [];
    for i = 1:n_site
        std_36Cl_er = [std_36Cl_er data_mc{i}.Cl36_uncer];
    end
    
    str_like = sprintf('log((2*pi)^-0.5.*sum(1./std_36Cl_er)* exp(-(ssfun([m(1) m(%d)],data_mc{1})',n_site+1);
    for i = 2:n_site
       str_like = [str_like sprintf(' + ssfun([m(%d) m(%d)],data_mc{%d})',i,n_site+1,i)];
    end   
    str_like = [str_like ' )/2));'];
    logLike = eval(sprintf('@(m) %s',str_like)); % the log likelihood
    
    % Make an initial guess for the model parameters.
    m0=[];
    for i = 1:n_site
        m0=[m0 ParamUser{i}.SR_0];
    end
    m0=[m0 ParamUser{1}.PG_age_0]';
    
    % Prior information
    str_test = ['(m(1)>ParamUser{1}.SRmin) & (m(1)<ParamUser{1}.SRmax)'];
    for i = 2:n_site
       str_test = [str_test sprintf(' & (m(%d)>ParamUser{%d}.SRmin) & (m(%d)<ParamUser{%d}.SRmax)',i,i,i,i)];
    end
    str_test = [str_test sprintf(' & (m(%d)>ParamUser{1}.Tmin) & (m(%d)<ParamUser{1}.Tmax);',n_site+1,n_site+1)];
    
    logprior = eval(sprintf('@(m) %s',str_test)); % the log prior
    
    
%% test_forward: compute theoretical 36Cl for a given model

    if(ParamUser{1}.test_forward == 1)
        fprintf( 1, '\t -> Test a model:\n')
        for i = 1:n_site
            fprintf( 1, '\t \t SR site %i = %f mm/yr\n',i,ParamUser{i}.sr);
        end
        fprintf( 1, '\t \t Post-glacial duration = %f yr\n',ParamUser{1}.tpg);
        
        RMSw_all = .0;
        
        for i = 1:n_site
            % get 36Cl concentration
            N_36_site = Model_direct_36Facet([ParamUser{i}.sr  ParamUser{i}.tpg],data_mc{i});
            % get the RMSw for each site
            RMSw_site = ssfun(m0,data_mc{i});
            % Get the total RMSw of all site
            RMSw_all = RMSw_all + RMSw_site;
            fprintf( 1, '\t -> RMSw site %i=%f \n',i,RMSw_site);
            % Plot 36Cl profile
            plot_36Cl_profile(N_36_site,Data{i},ParamUser{i});
        end    
        
        fprintf( 1, '\t -> RMSw=%f \t Likelog = %f\t\n',RMSw_all,logLike(m0));
        
    else % If inversion
    
    %% Initial sampling
    ball=zeros(n_site+1,ParamUser{1}.n_walker); % initial guess models
    tmp=zeros(n_site+1,1);
    disp('Initializing Inversion')
    disp('First we initialize the ensemble of walkers')
    disp(['-> number of walkers: ' num2str(ParamUser{1}.n_walker)])
    disp('-> guess models: ')
    i=0;
    while i<ParamUser{1}.n_walker
        for k=1:n_site
            tmp(k)=ParamUser{k}.SRmin+rand(1,1).*(ParamUser{k}.SRmax-ParamUser{k}.SRmin); % slip-rate of each site 
        end
            tmp(n_site+1)=ParamUser{1}.Tmin+rand(1,1).*(ParamUser{1}.Tmax-ParamUser{1}.Tmin); % common post-glacial duration
            
        if(logprior(tmp)==1) % if the model respects the prior
            i=i+1;
            ball(:,i)=tmp(:);
            fprintf ( 1, ' model:%i -> ',i);
            for k=1:n_site
                fprintf ( 1, '\t SR site %d = %4.2f',k,ball(k,i));
            end
                fprintf ( 1, '\t Tpg = %5.0f\n',ball(n_site+1,i));

        end

    end

    
    %% Inversion
    fprintf ( 1, ' \n\n Starting Inversion \n  ');
    m=gwmcmc(ball,{logprior logLike},ParamUser{1}.n_models_inversion,'burnin',0,'Parallel',ParamUser{1}.parallel_computing);
    
    %% save all
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
