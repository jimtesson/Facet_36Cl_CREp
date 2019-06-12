function [ ] = Inversion_36Cl_Facet()
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
    Data = Load_data_36('INPUT/DATA_IN.xlsx');
    
%% Parameters given by user
    NumGMDB = 1; % 1: Mush; 2: GLOPIS; 3: LSD; 4: own user geomagnetic db
    Scheme = 3; % time dependant scaling model (1: LAL-STONE with cutoff rigidity, 2: LSD, 3: LAL-STONE 2000 no cutoff)
    Muon_model = 1; % 1: Exponential, 2: numeric integration of the flux (Balco 2008)
    Atm = 0; % 0: ERA40 (Uppala et al. 2005), 1: standard atmosphere equation (NOAA 1976)
    Data.GMDB = 1;
    
    Age_max = 300000; % maximum age (yr) to compute 36Cl produced in sample
    
    PG_age_0 = 20000; % initial guess for post-glacial age (yr)
    Denud_0 = 1.0; % initial guess for the fault slip-rate (mm/yr)
    
    SRmin = 0; % Inversion minimum slip-rate bound (mm/yr)
    SRmax = 5; % Inversion maximum slip-rate bound (mm/yr)
    SR_std = 1;
    
    Tmin = 5000; % Inversion minimum post-glacial age bound (yr)
    Tmax = 30000; % Inversion maximim post-glacial age bound (yr)
    T_std = 2000;
    
    n_walker=10; % number of chain
    n_models_inversion = 20000; % number of models generated during the inversion
    parallel_computing = true; % inversion is operated using parallel computing (true, false)
    
    test = 0; % test a forward model ?
    % input parameters to test a forward model
    sr = 2.8; % mm/yr
    tpg = 32000; % yr
    
%% Some constants
    % 36Cl
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

%% Set Atmospheric and Geomagnetic parameters

    disp('Computing parameters')
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

%% Set the muon model
    if(Muon_model == 1)
        flag.muon = 'exp';  % Muon attenuation length approximated by exponential (Schimmelfenning 2009).  
   elseif(Muon_model == 2)
        flag.muon = 'num'; % Muon attenuation length calculated following Heisinger (2002a,b) and Balco et al. (2008)
    else
        Mess=sprintf('Invlid muon model');
    end
%% Set the Scaling model
    if(Scheme == 1)
        flag.scaling_model = 'st'; % LAL-STONE scheme cutoff rigidy
    elseif(Scheme == 2)
        flag.scaling_model = 'sa'; % LSD scheme
    elseif(Scheme == 3)
        flag.scaling_model = 'st2000'; % LAL-STONE scheme 2000 (no cuttoff)
    else
        Mess=sprintf('Invlid scaling model');
    end  
%% Constants initialization
    [Const_cosmo] = Init_const(Data,flag);   
        
%% Formatting data
    [Data_formatted,Param_site] = Init_var(Data);

%% Compute scaling factors
    %w = 0.2; % water content for Sato & Niita (2006)
    w = -1; % water content =  default value
    flag.NumGMDB = NumGMDB;
    Sf = Func_Sf(Param_site,Data,Atm,w,SelGMDB,flag);
    
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
    flag.plot = 0;
    flag.plotP = 0;
    
%% Packing variables
% Packing all required variables in a single cell for simplicity
    data_mc.Const_cosmo = Const_cosmo;
    data_mc.Param_cosmo = Param_cosmo;
    data_mc.Param_site = Param_site;
    data_mc.Sf = Sf;
    data_mc.Shf = Shf;
    data_mc.dataset = dataset;
    data_mc.Z_samples = Z_samples;
    data_mc.Age_max = Age_max;
    data_mc.flag = flag;
    data_mc.Cl36 = dataset(1,:);
    data_mc.flag = flag;
    
%% Parameters of the MCMc Inversion

    % Anonymous functions used by the inversion procedure:
    % function of the model -> provide modeled concentrations
    forwardmodel = @(x0,data_mc) Model_direct_36Facet(x0,data_mc); 

    % function providing the misfit
    ssfun = @(x0,data_mc) Get_misfit(x0,data_mc);
    
    % Likelihood: We assume the data are normally distributed around the forward model.
    % First we define a helper function equivalent to calling log(normpdf(x,mu,sigma))
    % but has higher precision because it avoids truncation errors associated with calling 
    % log(exp(xxx)).
    lognormpdf=@(x,mu,sigma)-0.5*((x-mu)./sigma).^2  -log(sqrt(2*pi).*sigma);
    logLike=@(m)sum(lognormpdf(data_mc.Cl36,forwardmodel(m,data_mc),exp(m(3))));

    % Make an initial guess for the model parameters.
    m0=[Denud_0 PG_age_0]';
    % 36Cl error is also a parameter
    sigma=std(data_mc.Cl36-forwardmodel(m0,data_mc));
    m0=[m0 ; log(sigma)];

    % Prior information
    logprior = @(m)(m(1)>SRmin)&(m(1)<SRmax)&(m(2)>Tmin)&(m(2)<Tmax);
    
%% Test: compute theoretical 36Cl for a given model
    disp('Test forward model')
    if(test == 1)
        % test ssfun
        ssfun(m0,data_mc)
        % test logLike
        logLike(m0)
        % model 36Cl
        N_36 = Model_direct_36Facet([sr tpg],data_mc);
        
        figure(3)
        % plot modelled 36Cl concentrations
        plot(N_36,Data.Alt-Data.alt_scarp_top,'o'); hold on
        % plot measured 36Cl concentrations
        plot(Data.NuclCon,Data.Alt-Data.alt_scarp_top,'o')
        return
    end
%% Initial sampling
    m_var(1)=SR_std;
    m_var(2)=T_std;
    m_var(3)=log(sigma*0.1); 
    
    ball=zeros(3,n_walker); 
    tmp=zeros(3,1);
    disp('Initializing Inversion')
    disp('First we initialize the ensemble of walkers in a small gaussian ball')
    disp(['number of walkers: ' num2str(n_walker)])

    i=0;
    while i<n_walker
        tmp(1)=m0(1)+randn(1,1).*m_var(1);
        tmp(2)=m0(2)+randn(1,1).*m_var(2);
        tmp(3)=m0(3)+randn(1,1).*m_var(3); 
        if(logprior(tmp)==1)
            i=i+1;
            ball(1,i)=tmp(1);
            ball(2,i)=tmp(2);
            ball(3,i)=tmp(3); 
            disp([num2str(i) ' ' num2str(ball(1,i)) ' ' num2str(ball(2,i)) ' ' num2str(ball(3,i))])
        end

    end



%% Apply the hammer:
%
% Draw samples from the posterior. 
%
    disp(['Starting Inversion'])

    m=gwmcmc(ball,{logprior logLike},n_models_inversion,'burnin',0,'Parallel',parallel_computing);

    save('Results/results_gwmcmc.mat')
    
% Check Autocorrelation 
    figure
    [C,lags,ESS]=eacorr(m);
    plot(lags,C,'.-',lags([1 end]),[0 0],'k');
    grid on
    xlabel('lags')
    ylabel('autocorrelation');
    text(lags(end),0,sprintf('Effective Sample Size (ESS): %.0f_ ',ceil(mean(ESS))),'verticalalignment','bottom','horizontalalignment','right')
    title('Markov Chain Auto Correlation')

% Remove BurnIn period
    BurnIn = 0.2; % Percent of the chain removed, must be >= 0 and <100
    L_chain = length(m(1,1,:)); % length of the chains
    crop=ceil(L_chain*BurnIn); % number of models removed
    m2 = m; % copy the results
    m2(:,:,1:crop)=[]; % removed models
    

% Corner plot of parameters
    figure
    ecornerplot(m2,'ks',true,'color',[0 0 0],'names',{'SR' 'T_PG' 'log(\sigma)'})

%% Statistics
    fprintf ( 1, 'Statistics from the inversion:\n' );
    % SR
        SR_mean = mean(reshape(m2(1,:,:), 1, [])); % mm/yr
        SR_std = std(reshape(m2(1,:,:), 1, [])); % mm/yr

    % PG Age
        T_mean = mean(reshape(m2(2,:,:), 1, [])); % yr
        T_std = std(reshape(m2(2,:,:), 1, [])); % yr

    % [36Cl] Error
        Error_36_mean = exp(mean(reshape(m2(3,:,:), 1, []))); % at/gr
        Error_36_std = exp(std(reshape(m2(3,:,:), 1, []))); % at/gr
        
    % Print results    
        fprintf ( 1, '\t -> Fault slip-rate:\n' );
        fprintf ( 1, '\t\t mean : %2.3f mm/yr\n',SR_mean);
        fprintf ( 1, '\t\t st. dev. : %2.3f mm/yr\n',SR_std);
        
        fprintf ( 1, '\t -> Post-glacial duration:\n' );
        fprintf ( 1, '\t\t mean : %6.0f yr\n',T_mean);
        fprintf ( 1, '\t\t st. dev. : %6.0f yr\n',T_std);
        
 %% Plot 36Cl Profile
    fprintf ( 1, 'Generating plot \n');
        % initialize some variables
        Nb_samples = length(Z_samples); % number of cosmogenic samples
        N_stat = 1000; % number of samples randomly picked to draw 36Cl concentrations
        N36=zeros(1000,Nb_samples); % number of cosmogenic samples
        
        % flatten m
        m3 = m(:,:)';
        
        % Progress Bar anonymous function
        progress=@textprogress2;

        % Get 36Cl concentration from 1000 random models generated by the inversion
        for kk=1:N_stat
            r=ceil(rand*size(m3,1));
            model=forwardmodel(m3(r,:),data_mc);
            N36(kk,:)=model;
            progress(kk/N_stat)
        end
        
        % Let's plot the models
        figure;
        hold on;
        % error bar coordinates on measurments
        X = [Data.NuclCon-Data.NuclErr Data.NuclCon+Data.NuclErr]' ;
        Y = [Data.Alt(:) Data.Alt(:)]' ;
        plot(X,Y,'-k')
        for kk=1:Nb_samples
            % plot 36Cl concentration
            plot(Data.NuclCon(kk),Data.Alt(kk),'ko','MarkerFaceColor','black')
            % Compute pdf
            [F,X,~]=ksdensity(N36(kk,:));
            X=X([1,1:end,end]);F=[0,F,0];
            %Scale pdf to altitude plot
            dalt=5;
            F_sca = F .* dalt ./ max(F(:))+Data.Alt(kk);
            % plot pdf
            fill(X,F_sca,[0.5 0.5 0.5],'edgecolor','none')
        end
            % set xlim
            xl = xlim; xl(1)=.0; xlim(xl);
            % set transparency
            alpha(.5); 
            
            % save all
            save('Results/results_gwmcmc.mat')
end

function textprogress2(pct)
    persistent lastNchar lasttime starttime
    if isempty(lastNchar)||pct==0
        lasttime=cputime-10;starttime=cputime;lastNchar=0;
        pct=1e-16;
    end
    if pct==1
        fprintf('%s',repmat(char(8),1,lastNchar));lastNchar=0;
        return
    end
    if (cputime-lasttime>0.1)

        ETA=datestr((cputime-starttime)*(1-pct)/(pct*60*60*24),13);
        progressmsg=[183-uint8((1:40)<=(pct*40)).*(183-'*') ''];
        %progressmsg=['-'-uint8((1:40)<=(pct*40)).*('-'-'?') ''];
        %progressmsg=[uint8((1:40)<=(pct*40)).*'#' ''];
        %curmtxt=sprintf('% 9.3g\n',curm(1:min(end,20),1));
        %curmtxt=mat2str(curm);
        progressmsg=sprintf('\nPlot progress %5.1f%% [%s] %s\n',pct*100,progressmsg,ETA);

        fprintf('%s%s',repmat(char(8),1,lastNchar),progressmsg);
        drawnow;lasttime=cputime;
        lastNchar=length(progressmsg);
    end
end
