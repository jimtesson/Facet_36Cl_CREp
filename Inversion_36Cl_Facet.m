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
        test = 1;
    else
        
%% Initial sampling of models for inversion
    Ini_models=GenerateInitialSamples(ParamUser,logprior);

%% GWMCMc Inversion 

    disp(['Starting Inversion'])

    m=gwmcmc(Ini_models,{logprior logLike},ParamUser.n_models_inversion,'burnin',0,'Parallel',ParamUser.parallel_computing);

    save('Results/results_gwmcmc.mat')
    
% Check Autocorrelation 
    PlotAutocorrelation(m)
    
% Remove BurnIn period
    L_chain = length(m(1,1,:)); % length of the chains
    crop=ceil(L_chain*ParamUser.BurnIn); % number of models removed
    % crop chain
    m_crop = m; % copy the results
    m_crop(:,:,1:crop)=[]; % removed models
    % flatten m
    m_flat = m_crop(:,:)';
    
% Corner plot of parameters
    PlotCorner(m_crop)

%% Statistics
    n_site = length(m_crop(:,1,1))-1;
    ResultStat = GetStatistics(m_crop,n_site);
 
%% Plot results from the inversion
    Plot_inversion_results(1,m_flat,Data,ParamUser,data_mc)

%% save results
    save('Results/results_gwmcmc.mat','ResultStat','m','m_flat')
    
    end
end

function ResultStat = GetStatistics(m_crop,nsite)
    fprintf ( 1, 'Statistics from the inversion:\n' );
    % Slip-rate
    for i=1:nsite
        ResultStat.SR_mean(i) = mean(reshape(m_crop(i,:,:), 1, [])); % mm/yr
        ResultStat.SR_std(i) = std(reshape(m_crop(i,:,:), 1, [])); % mm/yr
        fprintf ( 1, 'Site %i Slip-rate: mean = %3.1f, std = %3.1f \n',i,ResultStat.SR_mean(i),ResultStat.SR_std(i));
    end
    
    % PG Age
      ResultStat.T_mean = mean(reshape(m_crop(end,:,:), 1, [])); % yr
      ResultStat.T_std = std(reshape(m_crop(end,:,:), 1, [])); % yr
      fprintf ( 1, 'Mean Tpg = %5.0f, std Tpg = %5.0f \n',ResultStat.T_mean,ResultStat.T_std);

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

function Plot_inversion_results(i_site,m,Data,ParamUser,data_mc)
% Plot the 36Cl concentrations modelled from the inversion

        % Anonymous functions used by the inversion procedure:
        % function of the model -> provide modeled concentrations
        forwardmodel = @(model,data) Model_direct_36Facet(model,data); 
        
        % Progress Bar anonymous function
        progress=@textprogress2;   
        
        % Get statistics
        SR_mean = mean(reshape(m(:,i_site), 1, [])); % mm/yr
        SR_std = std(reshape(m(:,i_site), 1, [])); % mm/yr
        
        T_mean = mean(reshape(m(:,end), 1, [])); % mm/yr
    
        % initialize some variables
        Nb_samples = length(Data.NuclCon); % number of cosmogenic samples
        if(ParamUser.N_stat>length(m)) 
            nstat = length(m); % number of samples randomly picked to draw 36Cl concentrations
        else
            nstat = ParamUser.N_stat;
        end
        
        N36=zeros(nstat,Nb_samples); % number of cosmogenic samples
        
        % Get 36Cl concentration from 1000 random models generated by the inversion
        for kk=1:ParamUser.N_stat
            r=ceil(rand*size(m,1));
            m_temp = m(r,:);
            model=forwardmodel([m_temp(i_site) m_temp(end)],data_mc);
            N36(kk,:)=model;
            progress(kk/ParamUser.N_stat)
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
            %plot(N36(:,kk),ones(length(N36(:,kk))).*Data.Alt(kk),'ro')
            % Compute pdf
            [F,X,~]=ksdensity(N36(:,kk));
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
            title(sprintf('Site %i: Mean slip-rate = %2.1f mm/yr +/- %2.1f, Mean Post-glacial duration = %2.0f kyr',i_site,SR_mean,SR_std,T_mean/1000))
            xlabel('36Cl concentration (at/gr)') 
            ylabel('Altitude of the sample') 
            saveas(gcf,['Results/36Cl_PDF_' num2str(i_site) '.fig'])
 
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

function PlotAutocorrelation(m)
    figure
    [C,lags,ESS]=eacorr(m);
    plot(lags,C,'.-',lags([1 end]),[0 0],'k');
    grid on
    xlabel('lags')
    ylabel('autocorrelation');
    text(lags(end),0,sprintf('Effective Sample Size (ESS): %.0f_ ',ceil(mean(ESS))),'verticalalignment','bottom','horizontalalignment','right')
    title('Markov Chain Auto Correlation')
    saveas(gcf,'Results/ACM.fig')
end

function PlotCorner(m)
    figure
    ecornerplot(m,'ks',true,'color',[0 0 0],'names',{'SR' 'T_PG' 'log(\sigma)'})
    saveas(gcf,'Results/PDF.fig')
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