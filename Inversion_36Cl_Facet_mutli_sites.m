function [ ] = Inversion_36Cl_Facet_mutli_sites()
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
    
    % Log-Likelihood:
%     logLike=@(m)log((2*pi)^-0.5.*sum(1./[data_mc_1.Cl36_uncer data_mc_2.Cl36_uncer])* exp(-(ssfun([m(1) m(5)],data_mc_1) ...
%               +ssfun([m(2) m(5)],data_mc_2) ...
%               +ssfun([m(3) m(5)],data_mc_3) ... 
%               +ssfun([m(4) m(5)],data_mc_4) )/2));
%   

    str_like = sprintf('log((2*pi)^-0.5.*sum(1./std_36Cl_er)* exp(-(ssfun([m(1) m(%d)],data_mc{1})',n_site+1);
    for i = 2:n_site
       str_like = [str_like sprintf(' + ssfun([m(%d) m(%d)],data_mc{%d})',i,n_site+1,i)];
    end   
    str_like = [str_like ' )/2));'];
    logLike = eval(sprintf('@(m) %s',str_like)); % the log likelihood
    
    % Make an initial guess for the model parameters.
    %    m0=[ParamUser{1}.SR_0 ParamUser_2.SR_0 ParamUser_3.SR_0 ParamUser_4.SR_0 ParamUser{1}.PG_age_0]';
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
        % ParamUser.test_forward ssfun
        %ssfun(m0,data_mc);
        % ParamUser.test_forward logLike
        %logLike(m0);
        % model 36Cl
        
        RMSw_all = .0;
        
        for i = 1:n_site
        N_36_site = Model_direct_36Facet([ParamUser{i}.sr  ParamUser{i}.tpg],data_mc{i});
        RMSw_site = ssfun(m0,data_mc{i});
        RMSw_all = RMSw_all + RMSw_site;
        fprintf( 1, '\t -> RMSw site %i=%f \n',i,RMSw_site);
        % Figures: 36Cl profiles
        plot_36Cl_profile(N_36_site,Data{i},ParamUser{i});
        end
%         N_36_site_2 = Model_direct_36Facet([ParamUser_2.sr  ParamUser_2.tpg],data_mc_2);
%         N_36_site_3 = Model_direct_36Facet([ParamUser_3.sr  ParamUser_3.tpg],data_mc_3);
%         N_36_site_4 = Model_direct_36Facet([ParamUser_4.sr  ParamUser_4.tpg],data_mc_4);
%         
        
        fprintf( 1, '\t -> RMSw=%f \t Likelog = %f\t\n',RMSw_all,logLike(m0));

        
%         plot_36Cl_profile(N_36_site_1,Data{1},ParamUser_1);
%         plot_36Cl_profile(N_36_site_2,Data_2,ParamUser_2);
%         plot_36Cl_profile(N_36_site_3,Data_3,ParamUser_3);
%         plot_36Cl_profile(N_36_site_4,Data_4,ParamUser_4);
        
    else
%% Initial sampling
    
    ball=zeros(n_site+1,ParamUser{1}.n_walker); 
    tmp=zeros(n_site+1,1);
    disp('Initializing Inversion')
    disp('First we initialize the ensemble of walkers')
    disp(['-> number of walkers: ' num2str(ParamUser{1}.n_walker)])
    disp('-> guess models: \n')
    i=0;
    while i<ParamUser{1}.n_walker
        for k=1:n_site
            tmp(k)=ParamUser{k}.SRmin+rand(1,1).*(ParamUser{k}.SRmax-ParamUser{k}.SRmin); % slip-rate of each site 
        end
            tmp(n_site+1)=ParamUser{1}.Tmin+rand(1,1).*(ParamUser{1}.Tmax-ParamUser{1}.Tmin); % common post-glacial duration
            
        if(logprior(tmp)==1)
            i=i+1;
            ball(:,i)=tmp(:);
            fprintf ( 1, ' model:%i -> ',i);
            for k=1:n_site
                fprintf ( 1, '\t SR site %d = %4.2f',k,ball(k,i));
            end
                fprintf ( 1, '\t Tpg = %5.0f\n',ball(n_site+1,i));

        end

    end



%% Apply the hammer:
%
% Draw samples from the posterior. 
%
    disp(['Starting Inversion'])

    m=gwmcmc(ball,{logprior logLike},ParamUser{1}.n_models_inversion,'burnin',0,'Parallel',ParamUser{1}.parallel_computing);
    
    % save all
    save('Results/results_gwmcmc.mat')
    
    %% plot
    [ResultStat] = Plot_results_inversion();
    
    end
end

% function logLike = logLike(m,data_mc)
%     % function providing the misfit for a single site
%     ssfun = @(x0,data_mc) Get_misfit(x0,data_mc);
%     
%     n_site = length(data_mc);
%     logLike = .0;
%     
%     for i=1:n_site
%         logLike = logLike + ssfun([m(i) m(n_site)],data_mc{i});
%     end
%     
%     logLike= log((2*pi)^-0.5.*sum(1./[data_mc_1.Cl36_uncer data_mc_2.Cl36_uncer])* exp(-(logLike)/2));
%     
% end

function ResultStat = GetStatistics(m_crop,nsite)
    % Slip-rate
    for i=1:nsite
        ResultStat.SR_mean(i) = mean(reshape(m_crop(i,:,:), 1, [])); % mm/yr
        ResultStat.SR_std(i) = std(reshape(m_crop(i,:,:), 1, [])); % mm/yr
        fprintf ( 1, 'Site %i Slip-rate: mean = %3.1f, std = %3.1f \n',i,ResultStat.SR_mean(i),ResultStat.SR_std(i));
    end
    
    % PG Age
      ResultStat.T_mean = mean(reshape(m_crop(nsite+1,:,:), 1, [])); % yr
      ResultStat.T_std = std(reshape(m_crop(nsite+1,:,:), 1, [])); % yr
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
