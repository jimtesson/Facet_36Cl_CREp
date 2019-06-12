function best_age = Find_age(v_thin, ice_thickness, t_vector,...
                            Const_cosmo, ...
                            Param_cosmo,Param_cosmo_ICE, ...
                            Param_site,Param_site_ICE,...
                            Sf, Sf_ICE, erosion, dataset, t_guess,flag)
% Function to obtain the best age of a single sample and a fixed erosion
% rate. It uses the fminsearchbnd function to find the age associated with
% the lowest misfit (here RMSw).
% INPUT :   
%           Const_cosmo : Constants for cosmogenic's calculation
%           Param_cosmo : Sample specific cosmogenic's variable previously
%                           calculated
%           Param_site : Site specific parameters.
%           Sf : Scaling factors.
%           dataset : vector containing the data [36Cl, 36Cl uncertainties, depth of the sample] .
%           Z : depth of the sample
%           erosion : user-provided erosion-rate in m/Myr.
%           t_guess : guess for the time =  starting point.
%           flag : user-provided flag to specify the model used.
%
%                   flag.model = 'exp' for the Schimmelpfennig ?2009 model
%                   (attenuation length approximated by exponential)
%
%                   flag.model = 'num' for the numerical approach (Marrero
%                   et al. 2018, attenuation length are calculated from the
%                   energy flux))
%
%                   flag.scaling_model = 'st' for the Stone 1999 scaling 
%                                               scheme.
%
%                   flag.scaling_model = 'sa' for the Lifton-Sato scaling 
%                                               scheme.
%                   flag.search = 'fminsearch' for the 'fminsearch' method. 
%                   flag.search = 'nmsmax' for the Nelder-Mead simplex
%                                   method.
%                                   Required: 
%                                       - flag.min_bounds: minimum bound
%                                       - flag.max_bounds: maximun bound
%                                 
%                   flag.search = 'fminbnd' for the 'fminsearch' method
%                                   using boundaries.
%                                               
%
%
%           
% Output : 
%            best_age: output age for the sample
%
% version 01/08/2018, written by TESSON J.

% Function to obtain the modeled 36Cl as function of the erosion rate and the exposure age
Zsamples = dataset(3,:);


f = @(t_expo,erosion) inv_misfit(dataset,t_vector,...
                                gener_z_vector(t_vector,t_expo,v_thin,ice_thickness),...
                                Const_cosmo,Param_cosmo,Param_cosmo_ICE,...
                                Param_site,Param_site_ICE,...
                                Sf,Sf_ICE,Zsamples,erosion,flag); 
guess = t_guess; % guess time exposure


% Search
if(strcmp(flag.search,'fminsearch')==1)
    best_age = fminsearch(@(x) f(x,erosion),guess);
elseif(strcmp(flag.search,'nmsmax')==1)
    precx = 10^-4;%'Simplex' accuracy
    [best_age,~,~] = nmsmax(@(x) -f(x,erosion),guess,[precx inf inf 0 0]);
elseif(strcmp(flag.search,'fminbnd')==1)
    min = flag.min_bounds; % minimum bound for the search
    max = flag.max_bounds; % maximum bound for the search
    best_age = fminbnd(@(x) f(x,erosion),min,max);
else
    stop
end

end

function z_history = gener_z_vector(t_vector,t_degla,v_thin,ice_thickness)
            % t_vector : vector of time
            % t_degla : time when the sample start to be completly
            % uncovered by ice
            % 
            % create vector of depth history
            z_history = (t_vector-t_degla) .* v_thin; % z is equal to time * thining rate
            i_gla = find(z_history>ice_thickness); % find all z > ice thickness
            z_history(i_gla)=ice_thickness; % all z > thickness are set equal to the thickness
            i_t_degla = find(t_vector<=t_degla);
            z_history(i_t_degla)=0.0;
end

function RMSw = inv_misfit(dataset,t_vector,z_history,...
                            Const_cosmo,Param_cosmo,Param_cosmo_ICE, ...
                            Param_site,Param_site_ICE,...
                            Sf,Sf_ICE,...
                            erosion,Z,flag)
                        
                        
                        
% Function to calculate the sample concentration and return the RMSw.

[N_depth_profile] = depth_profile_speed(t_vector,z_history,...
                            Const_cosmo,Param_cosmo,Param_cosmo_ICE, ...
                            Param_site,Param_site_ICE,...
                            Sf,Sf_ICE,...
                            erosion,Z,flag);

RMSw = sum((((dataset(1,:)-N_depth_profile).^2).^.5)./dataset(2,:));
end
