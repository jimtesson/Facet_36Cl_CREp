function [ResultStat] = Plot_results_inversion()
% Function of the crep program that calculates exposure ages
    disp('Importing data')
%% Initialize path and functions 
    addpath(genpath('Functions'));
    addpath Functions
    addpath Constants
    addpath(genpath('Functions/36Cl_Functions'));
 
%% Import data from the inversion
    load('Results/results_gwmcmc.mat')
    
%% Loading data
    [Data,ParamUser] = Load_data_36('Input/DATA_IN.xlsx');
    % number of site 
    n_site = Data{1}.n_site;
    % string name of the site
    nameSite = cell(1,n_site+1);
    for i=1:n_site
      nameSite{i} = Data{i}.site_name;
    end
    nameSite{n_site+1} = 'Tpg';
    
% Check Autocorrelation 
    PlotAutocorrelation(m)
    
% Remove BurnIn period
    L_chain = length(m(1,1,:)); % length of the chains
    crop=ceil(L_chain*ParamUser{1}.BurnIn); % number of models removed
    % crop chain
    m_crop = m; % copy the results
    m_crop(:,:,1:crop)=[]; % removed models
    % flatten m
    m_flat = m_crop(:,:)';
    
% Corner plot of parameters
    PlotCorner(m_crop,nameSite)

%% Statistics

    ResultStat = GetStatistics(m_crop,n_site);
 
%% Plot results from the inversion
    for i =1: n_site
    Plot_inversion_results(i,m_flat,Data{i},ParamUser{i},data_mc{i})
    end
end

function ResultStat = GetStatistics(m_crop,nsite)
    fprintf ( 1, 'Statistics from the inversion:\n' );
    % Slip-rate
    for i=1:nsite
        ResultStat.SR_mean(i) = mean(reshape(m_crop(i,:,:), 1, [])); % mm/yr
        ResultStat.SR_std(i) = std(reshape(m_crop(i,:,:), 1, [])); % mm/yr
        fprintf ( 1, '-> Site %i \t Slip-rate: mean = %3.1f mm/yr, std = %3.1f mm/yr \n',i,ResultStat.SR_mean(i),ResultStat.SR_std(i));
    end
    
    % PG Age
      ResultStat.T_mean = mean(reshape(m_crop(end,:,:), 1, [])); % yr
      ResultStat.T_std = std(reshape(m_crop(end,:,:), 1, [])); % yr
      fprintf ( 1, '\t \t Mean Tpg = %5.0f yr, std Tpg = %5.0f yr \n',ResultStat.T_mean,ResultStat.T_std);

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

function PlotCorner(m,nameSite)
    figure
    ecornerplot(m,'ks',true,'color',[0 0 0],'names',nameSite)
    saveas(gcf,'Results/PDF.fig')
end
