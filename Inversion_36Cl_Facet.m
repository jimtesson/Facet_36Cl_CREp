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
    [Data,ParamUser] = Load_data_36('Input/DATA_IN.xlsx');
    
%% Set ParamUser.Atmospheric and Geomagnetic parameters

    disp('Computing parameters')
    Data.Atm = ParamUser.Atm;% ParamUser.Atmospheric model
        
    % Geomagnetic data base
    if length(Data.GMDB)==1;
        %ParamUser.NumGMDB = 2; % 1: Mush; 2: GLOPIS; 3: LSD; 4: own user geomagnetic db
        if ParamUser.NumGMDB==1;
            SelGMDB=GMDB.Musch;
        elseif ParamUser.NumGMDB==2;
            SelGMDB=GMDB.GLOPIS;
        else %  (ParamUser.NumGMDB=3)
            SelGMDB=GMDB.LSD;
        end
    else
        ParamUser.NumGMDB=4;
        SelGMDB=Data.GMDB;
    end

%% Set the muon model
    if(ParamUser.Muon_model == 1)
        flag.muon = 'exp';  % Muon attenuation length approximated by exponential (Schimmelfenning 2009).  
   elseif(ParamUser.Muon_model == 2)
        flag.muon = 'num'; % Muon attenuation length calculated following Heisinger (2002a,b) and Balco et al. (2008)
    else
        Mess=sprintf('Invlid muon model');
    end
%% Set the Scaling model
    if(ParamUser.Scheme == 1)
        flag.scaling_model = 'st'; % LAL-STONE ParamUser.Scheme cutoff rigidy
    elseif(ParamUser.Scheme == 2)
        flag.scaling_model = 'sa'; % LSD ParamUser.Scheme
    elseif(ParamUser.Scheme == 3)
        flag.scaling_model = 'st2000'; % LAL-STONE ParamUser.Scheme 2000 (no cuttoff)
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
    flag.NumGMDB = ParamUser.NumGMDB;
    Sf = Func_Sf(Param_site,Data,ParamUser.Atm,w,SelGMDB,flag);
    
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
    data_mc.Age_max = ParamUser.Age_max;
    data_mc.flag = flag;
    data_mc.Cl36 = dataset(1,:);
    data_mc.Cl36_uncer = data_mc.dataset(2,:);
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
    

    %lognormpdf=@(x,mu,sigma)-0.5*((x-mu)./sigma).^2  -log(sqrt(2*pi).*sigma);
    %logLike=@(m)sum(lognormpdf(data_mc.Cl36,forwardmodel(m,data_mc),exp(m(3))));
    
    logLike=@(m)log((2*pi)^-0.5.*sum(1./data_mc.Cl36_uncer)* exp(-ssfun(m,data_mc)/2));

    % Make an initial guess for the model parameters.
    m0=[ParamUser.SR_0 ParamUser.PG_age_0]';
    % 36Cl error is also a parameter
    %sigma=std(data_mc.Cl36-forwardmodel(m0,data_mc));
    %m0=[m0 ; log(sigma)];

    % Prior information
    logprior = @(m)(m(1)>ParamUser.SRmin)&(m(1)<ParamUser.SRmax)&(m(2)>ParamUser.Tmin)&(m(2)<ParamUser.Tmax);
    
%% test_forward: compute theoretical 36Cl for a given model

    if(ParamUser.test_forward == 1)
        fprintf( 1, '\t -> Test a model:\n \t Slip-rate = %f\n\tPost-glacial duration = %f\n',ParamUser.sr,ParamUser.tpg);
        % ParamUser.test_forward ssfun
        %ssfun(m0,data_mc);
        % ParamUser.test_forward logLike
        %logLike(m0);
        % model 36Cl
        N_36 = Model_direct_36Facet([ParamUser.sr ParamUser.tpg],data_mc);
        
        fprintf( 1, '\t -> RMSw=%f \t Likelog = %f\t\n',ssfun(m0,data_mc),logLike(m0));

        figure
        % plot modelled 36Cl concentrations
        plot(N_36,Data.Alt,'o'); hold on
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
        saveas(gcf,'Results/36Cl_forwardmodel.fig')
    else
%% Initial sampling
    
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
        %tmp(3)=log(exp(m0(3))+(rand(1,1)-0.5).*sigma);
        if(logprior(tmp)==1)
            i=i+1;
            ball(1,i)=tmp(1);
            ball(2,i)=tmp(2);
            %ball(3,i)=tmp(3); 
            fprintf ( 1, 'model:%i -> Slip-rate = %4.2f, Post-glacial duration = %6.0f\n',i,ball(1,i),ball(2,i));

        end

    end



%% Apply the hammer:
%
% Draw samples from the posterior. 
%
    disp(['Starting Inversion'])

    m=gwmcmc(ball,{logprior logLike},ParamUser.n_models_inversion,'burnin',0,'Parallel',ParamUser.parallel_computing);

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
    ParamUser.BurnIn = 0.2; % Percent of the chain removed, must be >= 0 and <100
    L_chain = length(m(1,1,:)); % length of the chains
    crop=ceil(L_chain*ParamUser.BurnIn); % number of models removed
    m2 = m; % copy the results
    m2(:,:,1:crop)=[]; % removed models
    

% Corner plot of parameters
    figure
    ecornerplot(m2,'ks',true,'color',[0 0 0],'names',{'SR' 'T_PG' 'log(\sigma)'})
    saveas(gcf,'Results/PDF.fig')
%% Statistics
    fprintf ( 1, 'Statistics from the inversion:\n' );
    % SR
        SR_mean = mean(reshape(m2(1,:,:), 1, [])); % mm/yr
        ParamUser.SR_std = std(reshape(m2(1,:,:), 1, [])); % mm/yr

    % PG Age
        T_mean = mean(reshape(m2(2,:,:), 1, [])); % yr
        ParamUser.T_std = std(reshape(m2(2,:,:), 1, [])); % yr

    % [36Cl] Error
%         Error_36_mean = exp(mean(reshape(m2(3,:,:), 1, []))); % at/gr
%         Error_36_std = exp(std(reshape(m2(3,:,:), 1, []))); % at/gr
        
    % Print results    
        fprintf ( 1, '\t -> Fault slip-rate:\n' );
        fprintf ( 1, '\t\t mean : %2.3f mm/yr\n',SR_mean);
        fprintf ( 1, '\t\t st. dev. : %2.3f mm/yr\n',SR_std);
        
        fprintf ( 1, '\t -> Post-glacial duration:\n' );
        fprintf ( 1, '\t\t mean : %6.0f yr\n',T_mean);
        fprintf ( 1, '\t\t st. dev. : %6.0f yr\n',T_std);
        
%         fprintf ( 1, '\t -> 36Cl error (x10^5):\n' );
%         fprintf ( 1, '\t\t mean : %6.0f at. 36Cl/gr\n',Error_36_mean*1E-5);
        
 %% Plot 36Cl Profile
    fprintf ( 1, 'Generating plot \n');
        
        % flatten m
        m3 = m(:,:)';
        
        % initialize some variables
        Nb_samples = length(Z_samples); % number of cosmogenic samples
        if(ParamUser.N_stat>length(m3)) 
            ParamUser.N_stat = length(m3); % number of samples randomly picked to draw 36Cl concentrations
        end
        N36=zeros(1000,Nb_samples); % number of cosmogenic samples
        
        % Progress Bar anonymous function
        progress=@textprogress2;

        % Get 36Cl concentration from 1000 random models generated by the inversion
        for kk=1:ParamUser.N_stat
            r=ceil(rand*size(m3,1));
            model=forwardmodel(m3(r,:),data_mc);
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
            title(sprintf('Mean slip-rate = %2.1f mm/yr, Mean Post-glacial duration = %2.0f kyr',SR_mean,T_mean/1000))
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
