function [ N_calc_tot] = depth_profile_speed(...
                   t_vector,z_history,Const_cosmo,...
                   Param_cosmo_in,...
                   Param_site_in,...
                   Sf_in,Shf_in,...
                   z_vector,flag)
               
% This function computes the theoretical 36Cl for a given vector of samples
% , an fixed erosion rate and exposure duration. User can choose between 2 
%   models: 1) the Schimmelfenning 2009 model, 2) numerical approach using 
%   the LSD model.
% INPUT :   
%           Const_cosmo : Constants for cosmogenic's calculation.
%           Param_cosmo_in : Struct. containing the samples specific 
%                            variable previously calculated.                          
%           Param_site_in : Site specific parameters.
%           sf : Scaling factors, previously calculated as function of time
%                   and depth.
%           z_vector : depth vector of the samples
%           t_expo : exposure age of the samples
%           erosion : user-provided erosion-rate in m/Myr
%           flag : user-provided flag to specify the model used.
%
%                   flag.model = 'exp' for the Schimmelpfennig 2009 model
%                   (attenuation length approximated by exponential)
%
%                   flag.model = 'num' for the numerical approach (Marrero
%                   et al. 2018, attenuation length are calculated from the
%                   energy flux))
%
%                   flag.scaling_model = 'st' for the Stone 1999 scaling 
%                                               scheme
%
%                   flag.scaling_model = 'sa' for the Lifton-Sato scaling 
%                                               scheme
%
%
%           
% Output : 
%             N_calc_tot: 36Cl concentration for each sample
%             N_calc_tot_uncert: uncertainties on the 36Cl.
% version 01/08/2018, written by TESSON J.

N_samples = length(z_vector);

    % Time and position increment
    deltat = 100; % (yr)
    
    % Theoretical 36Cl vector
    N_calc_tot=zeros(1,N_samples);
   
    for i=1:length(z_vector)
            
            % depth-time vector
            depth_time_vector = z_history + z_vector(i); %t_vector .* erosion + z_vector(i); % evolution of the sample depth (cm) over the time
            
            % Rock parameters
            Param_cosmo = Param_cosmo_in{i};
            Param_site = Param_site_in{i};
            
            % Scaling factors
            Sf = Sf_in{i};
            % Shielding factors
            Shf = Shf_in{i};
            
            % get current scaling factors.
            currentsf=getcurrentsf(Sf,t_vector,flag); 
            Sf.currentsf = currentsf;
            
            % get Production rate within the sample at surface z=.0;
                % get production rate over the whole depth time vector, and
                % averaging over the whole sample.
         
                n_thick = 10;
                d_thick = Param_site.thick/n_thick.*(0:1:n_thick); 
                
                % Depth vectors integrating at the depth of the sample
                d_integ_samp = d_thick + z_vector(i); % z_vector = top of the sample
                depth_thick_matrix = zeros(n_thick,length(depth_time_vector)); % matrix for the 11 "depth-steps" of the sample

                for j=1:(n_thick+1)                                     
                    depth_thick_matrix(j,:) = depth_time_vector + d_integ_samp(j);
                end

                depth_thick_matrix(depth_thick_matrix>max(Shf.Z_sc(:)))=max(Shf.Z_sc(:));
                
                % get the shielding shielding factor for all depth of the
                % matrix 
                    L_sp = interp1(Shf.Z_sc,Shf.L_sp,depth_thick_matrix);
                    so_sp = Shf.so_sp;
                    L_mu = interp1(Shf.Z_sc,Shf.L_mu,depth_thick_matrix);
                    so_mu = Shf.so_mu;
                
                P36_tot = zeros(1,length(depth_time_vector));
                P36_s = zeros(1,length(depth_time_vector));
                P36_th = zeros(1,length(depth_time_vector));
                P36_eth = zeros(1,length(depth_time_vector));
                P36_mu = zeros(1,length(depth_time_vector));
                P36_rad = zeros(1,length(depth_time_vector));
%   tic              
%                 for j=1:length(d_integ_samp) % loop over the thickness
%                     
%                     % Extract the depth vector for the given depth from the matrix
%                     depth_vect = depth_thick_matrix(j,:);
%                     
%                     L_sp_j = L_sp(j,:);
%                     L_mu_j = L_mu(j,:);
%                     
%                     % get production rates
%                     [~,P36_s_tmp,P36_th_tmp,P36_eth_tmp,P36_mu_tmp,P36_rad_tmp] = prodz36_speed2(Const_cosmo,Param_cosmo,Sf,L_sp_j,so_sp,L_mu_j,so_mu,flag,depth_vect*Param_site.rho_rock);
%                     
%                     % total production rate within the sample at surface
%                     P36_tot_tmp = P36_s_tmp+P36_th_tmp+P36_eth_tmp+P36_mu_tmp+P36_rad_tmp;
%                     
%                     % Sum production rate over the sample thickness
%                     P36_tot = P36_tot + P36_tot_tmp;
%                     
%                     % details for each pathway
%                     P36_s = P36_s + P36_s_tmp;
%                     P36_th = P36_th + P36_th_tmp;
%                     P36_eth = P36_eth + P36_eth_tmp;
%                     P36_mu = P36_mu + P36_mu_tmp;
%                     P36_rad = P36_rad + P36_rad_tmp;
%                     
%                 end
% 
%                 % Average production rate over the sample thickness
%                 P36_tot = P36_tot./length(d_integ_samp); 
%                 P36_s = P36_s./length(d_integ_samp);
%                 P36_th = P36_th./length(d_integ_samp);
%                 P36_eth = P36_eth./length(d_integ_samp);
%                 P36_rad = P36_rad./length(d_integ_samp);
%                 P36_mu = P36_mu./length(d_integ_samp);
% toc                

                [~,P36_s_tmp,P36_th_tmp,P36_eth_tmp,P36_mu_tmp,P36_rad_tmp] = prodz36_speed2(Const_cosmo,Param_cosmo,Sf,L_sp,so_sp,L_mu,so_mu,flag,depth_thick_matrix*Param_site.rho_rock);
                P36_tot_tmp = P36_s_tmp+P36_th_tmp+P36_eth_tmp+P36_mu_tmp+P36_rad_tmp;
                P36_tot=mean(P36_tot_tmp);
     
                if(flag.plotP)
                P36_s_percent = P36_s_tmp./P36_tot.*100;
                P36_th_percent = P36_th_tmp./P36_tot.*100;
                P36_eth_percent = P36_eth_tmp./P36_tot.*100;
                P36_mu_percent = P36_mu_tmp./P36_tot.*100;
                P36_rad_percent = P36_rad_tmp./P36_tot.*100;
                disp(['Sample :',num2str(i)])
                disp(['   Production rate at surface(at/gr) :',num2str(P36_tot(length(P36_tot)))])
                disp(['   Production pathway :',' Spal (%):' num2str(mean(P36_s_percent(:))) ' Th (%):' num2str(mean(P36_th_percent(:))) ' Eth (%):' num2str(mean(P36_eth_percent(:))) ...
                      ' Muon (%):' num2str(mean(P36_mu_percent(:))) ' Rad (%):' num2str(mean(P36_rad_percent(:)))])
                disp(['--'])
                end
                

         % Compute Total concentration
            N_inh=.0;         
            N_calc_tot(i)=.0;
            N36=zeros(1,length(t_vector));
           PP=zeros(1,length(t_vector));
           
           f1=exp(-Const_cosmo.lambda36*deltat);
           f2=(1.0-exp(-Const_cosmo.lambda36*deltat))./Const_cosmo.lambda36;
            P=P36_tot;
          for it=1:length(t_vector)
              N_calc_tot(i)=N_calc_tot(i).*f1 + P(it).*f2;
              N36(it)=N_calc_tot(i);
              PP(it)=P(it);
          end

% tic
%         V1 = exp(-Const_cosmo.lambda36.*fliplr(t_vector));
%         V2= deltat * cumsum(P ./ V1);
%         N36 = V1 .* (N_inh + V2);
%         N_calc_tot(i) = N36(length(t_vector));
% 
% toc          

        if(flag.plot)
            figure(1);  
      
            subplot(2,2,1); 
            plot(-t_vector,P36_tot); hold on;
            
            subplot(2,2,2); 
            plot(-t_vector,depth_time_vector./100); hold on;
            
            subplot(2,2,3);  
            plot(-t_vector,N36);hold on;
            
            subplot(2,2,4);
            plot(-t_vector,N36./N_calc_tot(i).*100);hold on;
            
        end
    end
      


