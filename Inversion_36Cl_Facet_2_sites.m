function [ ] = Inversion_36Cl_Facet_2_sites()
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
    [Data_1,ParamUser_1] = Load_data_36('Input/DATA_IN_1.xlsx'); % MA3
    [Data_2,ParamUser_2] = Load_data_36('Input/DATA_IN_2.xlsx'); % MA1


    
%% Set ParamUser.Atmospheric and Geomagnetic parameters

    disp('Computing parameters')
    Data_1.Atm = ParamUser_1.Atm;% ParamUser.Atmospheric model
    Data_2.Atm = ParamUser_2.Atm;% ParamUser.Atmospheric model
        
    % Geomagnetic data base
    if length(Data_1.GMDB)==1
        %ParamUser.NumGMDB = 2; % 1: Mush; 2: GLOPIS; 3: LSD; 4: own user geomagnetic db
        if ParamUser_1.NumGMDB==1
            SelGMDB=GMDB.Musch;
        elseif ParamUser_1.NumGMDB==2
            SelGMDB=GMDB.GLOPIS;
        else %  (ParamUser.NumGMDB=3)
            SelGMDB=GMDB.LSD;
        end
    else
        ParamUser_1.NumGMDB=4;
        SelGMDB=Data_1.GMDB;
    end

%% Set the muon model
    if(ParamUser_1.Muon_model == 1)
        flag.muon = 'exp';  % Muon attenuation length approximated by exponential (Schimmelfenning 2009).  
   elseif(ParamUser_1.Muon_model == 2)
        flag.muon = 'num'; % Muon attenuation length calculated following Heisinger (2002a,b) and Balco et al. (2008)
    else
        Mess=sprintf('Invlid muon model');
    end
%% Set the Scaling model
    if(ParamUser_1.Scheme == 1)
        flag.scaling_model = 'st'; % LAL-STONE ParamUser.Scheme cutoff rigidy
    elseif(ParamUser_1.Scheme == 2)
        flag.scaling_model = 'sa'; % LSD ParamUser.Scheme
    elseif(ParamUser_1.Scheme == 3)
        flag.scaling_model = 'st2000'; % LAL-STONE ParamUser.Scheme 2000 (no cuttoff)
    else
        Mess=sprintf('Invlid scaling model');
    end  
%% Constants initialization
    [Const_cosmo_1] = Init_const(Data_1,flag);
    [Const_cosmo_2] = Init_const(Data_2,flag);
        
%% Formatting data
    [Data_formatted_1,Param_site_1] = Init_var(Data_1);
    [Data_formatted_2,Param_site_2] = Init_var(Data_2);

%% Compute scaling factors
    %w = 0.2; % water content for Sato & Niita (2006)
    w = -1; % water content =  default value
    flag.NumGMDB = ParamUser_1.NumGMDB;
    Sf_1 = Func_Sf(Param_site_1,Data_1,ParamUser_1.Atm,w,SelGMDB,flag);
    Sf_2 = Func_Sf(Param_site_2,Data_2,ParamUser_2.Atm,w,SelGMDB,flag);
    
%% Compute production rates
    Param_cosmo_1 = clrock(Data_formatted_1,Param_site_1,Const_cosmo_1,Sf_1);
    Param_cosmo_2 = clrock(Data_formatted_2,Param_site_2,Const_cosmo_2,Sf_2);
    
%% Compute Shielding factors
    Shf_1 = Shielding_facet(Param_site_1,Const_cosmo_1);
    Shf_2 = Shielding_facet(Param_site_2,Const_cosmo_2);
    
%% Prepare dataset
    % Site 1
    dataset_1(1,:) = Data_1.NuclCon(:);
    dataset_1(2,:) = Data_1.NuclErr(:);
    dataset_1(3,:) = Data_1.Z(:);
    Z_samples_1 = Data_1.Z(:); % depth of the samples (in cm)
    
    % Site 2
    dataset_2(1,:) = Data_2.NuclCon(:);
    dataset_2(2,:) = Data_2.NuclErr(:);
    dataset_2(3,:) = Data_2.Z(:);
    Z_samples_2 = Data_2.Z(:); % depth of the samples (in cm)
    
%% Set plot parameters
    flag.plot = 0;
    flag.plotP = 0;
    
%% Packing variables
% Packing all required variables in a single cell for simplicity
    % Site 1
    data_mc_1.Const_cosmo = Const_cosmo_1;
    data_mc_1.Param_cosmo = Param_cosmo_1;
    data_mc_1.Param_site = Param_site_1;
    data_mc_1.Sf = Sf_1;
    data_mc_1.Shf = Shf_1;
    data_mc_1.dataset = dataset_1;
    data_mc_1.Z_samples = Z_samples_1;
    data_mc_1.Age_max = ParamUser_1.Age_max;
    data_mc_1.flag = flag;
    data_mc_1.Cl36 = dataset_1(1,:);
    data_mc_1.Cl36_uncer = data_mc_1.dataset(2,:);
    data_mc_1.flag = flag;
    % Site 2
    data_mc_2.Const_cosmo = Const_cosmo_2;
    data_mc_2.Param_cosmo = Param_cosmo_2;
    data_mc_2.Param_site = Param_site_2;
    data_mc_2.Sf = Sf_2;
    data_mc_2.Shf = Shf_2;
    data_mc_2.dataset = dataset_2;
    data_mc_2.Z_samples = Z_samples_2;
    data_mc_2.Age_max = ParamUser_2.Age_max;
    data_mc_2.flag = flag;
    data_mc_2.Cl36 = dataset_2(1,:);
    data_mc_2.Cl36_uncer = data_mc_2.dataset(2,:);
    data_mc_2.flag = flag;
    
%% Parameters of the MCMc Inversion

    % Anonymous functions used by the inversion procedure:
    % function of the model -> provide modeled concentrations
    forwardmodel = @(x0,data_mc) Model_direct_36Facet(x0,data_mc); 

    % function providing the misfit for a single site
    ssfun = @(x0,data_mc) Get_misfit(x0,data_mc);
    
    % Log-Likelihood:
    logLike=@(m)log((2*pi)^-0.5.*sum(1./[data_mc_1.Cl36_uncer data_mc_2.Cl36_uncer])* exp(-(ssfun([m(1) m(3)],data_mc_1)+ssfun([m(2) m(3)],data_mc_2))/2));

    % Make an initial guess for the model parameters.
    m0=[ParamUser_1.SR_0 ParamUser_2.SR_0 ParamUser_1.PG_age_0]';

    % Prior information
    logprior = @(m)(m(1)>ParamUser_1.SRmin) & (m(1)<ParamUser_1.SRmax) ...
                 & (m(2)>ParamUser_2.SRmin) & (m(2)<ParamUser_2.SRmax)...
                 & (m(3)>ParamUser_1.Tmin) & (m(3)<ParamUser_1.Tmax);
    
%% test_forward: compute theoretical 36Cl for a given model

    if(ParamUser_1.test_forward == 1)
        fprintf( 1, '\t -> Test a model:\n \t SR site 1 = %f\n\t SR site 2 = %f\n\tPost-glacial duration = %f\n',ParamUser_1.sr,ParamUser_2.sr,ParamUser_1.tpg);
        % ParamUser.test_forward ssfun
        %ssfun(m0,data_mc);
        % ParamUser.test_forward logLike
        %logLike(m0);
        % model 36Cl
        N_36_site1 = Model_direct_36Facet([ParamUser_1.sr  ParamUser_1.tpg],data_mc_1);
        N_36_site2 = Model_direct_36Facet([ParamUser_2.sr  ParamUser_2.tpg],data_mc_2);
        
        N36 = [N_36_site1 N_36_site2];
        
        fprintf( 1, '\t -> RMSw=%f \t Likelog = %f\t\n',ssfun(m0,data_mc_1)+ssfun(m0,data_mc_2),logLike(m0));
        fprintf( 1, '\t -> RMSw site 1=%f \t RMSw site 1=%f\t\n',ssfun(m0,data_mc_1),ssfun(m0,data_mc_2));

        % site 1
        figure
        % plot modelled 36Cl concentrations
        plot(N_36_site1,Data_1.Alt,'o'); hold on
        % plot measured 36Cl concentrations
            % error bar coordinates on measurments
            X = [Data_1.NuclCon-Data_1.NuclErr Data_1.NuclCon+Data_1.NuclErr]' ;
            Y = [Data_1.Alt(:) Data_1.Alt(:)]' ;
        plot(Data_1.NuclCon,Data_1.Alt,'ko','MarkerFaceColor','black')
        plot(X,Y,'-k')
        ax_x = xlim; xlim([0 ax_x(2)])
        title(sprintf('Slip-rate = %2.1f mm/yr, Post-glacial duration = %2.0f kyr',ParamUser_1.sr,ParamUser_1.tpg/1000))
        xlabel('36Cl concentration (at/gr)') 
        ylabel('Altitude of the sample') 
        legend('model','data','Location','southeast')
        saveas(gcf,'Results/36Cl_forwardmodel.fig')
        
        % site 2
        figure
        % plot modelled 36Cl concentrations
        plot(N_36,Data_2.Alt,'o'); hold on
        % plot measured 36Cl concentrations
            % error bar coordinates on measurments
            X = [Data_2.NuclCon-Data_2.NuclErr Data_2.NuclCon+Data_2.NuclErr]' ;
            Y = [Data_2.Alt(:) Data_2.Alt(:)]' ;
        plot(Data_2.NuclCon,Data_2.Alt,'ko','MarkerFaceColor','black')
        plot(X,Y,'-k')
        ax_x = xlim; xlim([0 ax_x(2)])
        title(sprintf('Slip-rate = %2.1f mm/yr, Post-glacial duration = %2.0f kyr',ParamUser_2.sr,ParamUser_1.tpg/1000))
        xlabel('36Cl concentration (at/gr)') 
        ylabel('Altitude of the sample') 
        legend('model','data','Location','southeast')
        saveas(gcf,'Results/36Cl_forwardmodel.fig')
    else
%% Initial sampling
    
    ball=zeros(3,ParamUser_1.n_walker); 
    tmp=zeros(3,1);
    disp('Initializing Inversion')
    disp('First we initialize the ensemble of walkers')
    disp(['-> number of walkers: ' num2str(ParamUser_1.n_walker)])
    disp('-> guess models: ')
    i=0;
    while i<ParamUser_1.n_walker
        tmp(1)=ParamUser_1.SRmin+rand(1,1).*(ParamUser_1.SRmax-ParamUser_1.SRmin); % slip-rate site 1
        tmp(2)=ParamUser_2.SRmin+rand(1,1).*(ParamUser_2.SRmax-ParamUser_2.SRmin); % slip-rate site 2
        tmp(3)=ParamUser_1.Tmin+rand(1,1).*(ParamUser_1.Tmax-ParamUser_1.Tmin); % common post-glacial duration

        if(logprior(tmp)==1)
            i=i+1;
            ball(1,i)=tmp(1);
            ball(2,i)=tmp(2);
           ball(3,i)=tmp(3); 
            fprintf ( 1, 'model:%i -> SR 1 = %4.2f, SR 2 = %4.2f, Post-glacial duration = %6.0f\n',i,ball(1,i),ball(2,i),ball(3,i));

        end

    end



%% Apply the hammer:
%
% Draw samples from the posterior. 
%
    disp(['Starting Inversion'])

    m=gwmcmc(ball,{logprior logLike},ParamUser_1.n_models_inversion,'burnin',0,'Parallel',ParamUser_1.parallel_computing);

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
    saveas(gcf,'Results/ACM.fig')
    
% Remove BurnIn period
    ParamUser_1.BurnIn = 0.2; % Percent of the chain removed, must be >= 0 and <100
    L_chain = length(m(1,1,:)); % length of the chains
    crop=ceil(L_chain*ParamUser_1.BurnIn); % number of models removed
    m2 = m; % copy the results
    m2(:,:,1:crop)=[]; % removed models
    

% Corner plot of parameters
    figure
    ecornerplot(m2,'ks',true,'color',[0 0 0],'names',{'SR_1' 'SR_2' 'T_PG'})
    saveas(gcf,'Results/PDF.fig')
%% Statistics
    fprintf ( 1, 'Statistics from the inversion:\n' );
    % SR
        SR_1_mean = mean(reshape(m2(1,:,:), 1, [])); % mm/yr
        SR_1_std = std(reshape(m2(1,:,:), 1, [])); % mm/yr
     % SR
        SR_2_mean = mean(reshape(m2(2,:,:), 1, [])); % mm/yr
        SR_2_std = std(reshape(m2(2,:,:), 1, [])); % mm/yr

    % PG Age
        T_mean = mean(reshape(m2(3,:,:), 1, [])); % yr
        T_std = std(reshape(m2(3,:,:), 1, [])); % yr

    % [36Cl] Error
%         Error_36_mean = exp(mean(reshape(m2(3,:,:), 1, []))); % at/gr
%         Error_36_std = exp(std(reshape(m2(3,:,:), 1, []))); % at/gr
        
    % Print results    
        fprintf ( 1, '\t -> Fault slip-rate site 1:\n' );
        fprintf ( 1, '\t\t mean : %2.3f mm/yr\n',SR_1_mean);
        fprintf ( 1, '\t\t st. dev. : %2.3f mm/yr\n',SR_1_std);
        
        fprintf ( 1, '\t -> Fault slip-rate site 2:\n' );
        fprintf ( 1, '\t\t mean : %2.3f mm/yr\n',SR_2_mean);
        fprintf ( 1, '\t\t st. dev. : %2.3f mm/yr\n',SR_2_std);
        
        fprintf ( 1, '\t -> Post-glacial duration:\n' );
        fprintf ( 1, '\t\t mean : %6.0f yr\n',T_mean);
        fprintf ( 1, '\t\t st. dev. : %6.0f yr\n',T_std);
        
%         fprintf ( 1, '\t -> 36Cl error (x10^5):\n' );
%         fprintf ( 1, '\t\t mean : %6.0f at. 36Cl/gr\n',Error_36_mean*1E-5);
        
 %% Plot 36Cl Profile Site 1
    fprintf ( 1, 'Generating plot site 1 \n');
        
        % flatten m
        m3 = m(:,:)';
        
        % initialize some variables
        Nb_samples = length(Z_samples_1); % number of cosmogenic samples
        if(ParamUser_1.N_stat>length(m3)) 
            nstat = length(m3); % number of samples randomly picked to draw 36Cl concentrations
        else
            nstat = ParamUser_1.N_stat;
        end
        
        N36=zeros(nstat,Nb_samples); % number of cosmogenic samples
        
        % Progress Bar anonymous function
        progress=@textprogress2;

        % Get 36Cl concentration from 1000 random models generated by the inversion
        for kk=1:ParamUser_1.N_stat
            r=ceil(rand*size(m3,1));
            m_temp = m3(r,:);
            model=forwardmodel([m_temp(1) m_temp(3)],data_mc_1);
            N36(kk,:)=model;
            progress(kk/ParamUser_1.N_stat)
        end
        
        % Let's plot the models
        figure;
        hold on;
        % error bar coordinates on measurments
        X = [Data_1.NuclCon-Data_1.NuclErr Data_1.NuclCon+Data_1.NuclErr]' ;
        Y = [Data_1.Alt(:) Data_1.Alt(:)]' ;
        plot(X,Y,'-k')
        for kk=1:Nb_samples

            % plot 36Cl concentration
            plot(Data_1.NuclCon(kk),Data_1.Alt(kk),'ko','MarkerFaceColor','black')
            %plot(N36(:,kk),ones(length(N36(:,kk))).*Data.Alt(kk),'ro')
            % Compute pdf
            [F,X,~]=ksdensity(N36(:,kk));
            X=X([1,1:end,end]);F=[0,F,0];
            %Scale pdf to altitude plot
            dalt=5;
            F_sca = F .* dalt ./ max(F(:))+Data_1.Alt(kk);
            % plot pdf
            fill(X,F_sca,[0.5 0.5 0.5],'edgecolor','none')
        end
            % set xlim
            xl = xlim; xl(1)=.0; xlim(xl);
            % set transparency
            alpha(.5); 
            title(sprintf('Mean slip-rate = %2.1f mm/yr, Mean Post-glacial duration = %2.0f kyr',SR_1_mean,T_mean/1000))
            xlabel('36Cl concentration (at/gr)') 
            ylabel('Altitude of the sample') 
            saveas(gcf,'Results/36Cl_PDF.fig')

 %% Plot 36Cl Profile Site 2
    fprintf ( 1, 'Generating plot site 2 \n');
        
        % flatten m
        m3 = m(:,:)';
        
        % initialize some variables
        Nb_samples = length(Z_samples_2); % number of cosmogenic samples
        if(ParamUser_1.N_stat>length(m3)) 
            nstat = length(m3); % number of samples randomly picked to draw 36Cl concentrations
        else
            nstat = ParamUser_1.N_stat;
        end
        
        N36=zeros(nstat,Nb_samples); % number of cosmogenic samples
        
        % Progress Bar anonymous function
        progress=@textprogress2;

        % Get 36Cl concentration from 1000 random models generated by the inversion
        for kk=1:ParamUser_1.N_stat
            r=ceil(rand*size(m3,1));
            m_temp = m3(r,:);
            model=forwardmodel([m_temp(2) m_temp(3)],data_mc_2);
            N36(kk,:)=model;
            progress(kk/ParamUser_1.N_stat)
        end
        
        % Let's plot the models
        figure;
        hold on;
        % error bar coordinates on measurments
        X = [Data_2.NuclCon-Data_2.NuclErr Data_2.NuclCon+Data_2.NuclErr]' ;
        Y = [Data_2.Alt(:) Data_2.Alt(:)]' ;
        plot(X,Y,'-k')
        for kk=1:Nb_samples

            % plot 36Cl concentration
            plot(Data_2.NuclCon(kk),Data_2.Alt(kk),'ko','MarkerFaceColor','black')
            %plot(N36(:,kk),ones(length(N36(:,kk))).*Data.Alt(kk),'ro')
            % Compute pdf
            [F,X,~]=ksdensity(N36(:,kk));
            X=X([1,1:end,end]);F=[0,F,0];
            %Scale pdf to altitude plot
            dalt=5;
            F_sca = F .* dalt ./ max(F(:))+Data_2.Alt(kk);
            % plot pdf
            fill(X,F_sca,[0.5 0.5 0.5],'edgecolor','none')
        end
            % set xlim
            xl = xlim; xl(1)=.0; xlim(xl);
            % set transparency
            alpha(.5); 
            title(sprintf('Mean slip-rate = %2.1f mm/yr, Mean Post-glacial duration = %2.0f kyr',SR_2_mean,T_mean/1000))
            xlabel('36Cl concentration (at/gr)') 
            ylabel('Altitude of the sample') 
            saveas(gcf,'Results/36Cl_PDF.fig')
            
            % save all
            save('Results/results_gwmcmc.mat')
    end
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
