function [ N_36 ] = Model_direct_36Facet(x0,data_mc)
% Fonction to compute the sample pathway below the surface over the time
% for a given model (Post-glacial age, denudation rate). Then the pathway
% is used to compute 36Cl concentration
% 

%% Initialization
    % UnPacking all required variables in a single cell for simplicity
    Const_cosmo = data_mc.Const_cosmo;
    Param_cosmo = data_mc.Param_cosmo;
    Param_site = data_mc.Param_site;
    Sf = data_mc.Sf;
    Shf = data_mc.Shf;
    Z_samples = data_mc.Z_samples;
    Age_max = data_mc.Age_max;
    flag = data_mc.flag;

    % Attribute input parameter
    Age_PG = x0(2); 
    SR = x0(1);
                                    
    % Other variables
    dt = 100; % time step
        
%% Define the time-depth vectors, 
% describing the pathway of the sample 
% over the time from depth to surface

    Age_max = Age_max-Age_PG;
    T_expo_Glacial=fliplr((0:dt:Age_max-dt));
    SR = SR*1e-1; %(cm/yr)
    depth_vector_Gla = T_expo_Glacial .* SR; 
    T_expo_Glacial = T_expo_Glacial+Age_PG+dt;

    T_expo_PG=fliplr((0:dt:Age_PG));
    depth_vector_PG = T_expo_PG .* 0.0;
    T_expo = [T_expo_Glacial T_expo_PG];
    depth_vector = [depth_vector_Gla depth_vector_PG];

    t_vector = T_expo; % time vector
    z_history = depth_vector; % depth vector : position of the sample below the surface along the plane parallel to the fault-plane
    
%% Compute the theoretical 36Cl concentration

    N_36 = depth_profile_speed(t_vector,z_history,...
                            Const_cosmo,Param_cosmo, ...
                            Param_site,...
                            Sf,Shf,...
                            Z_samples,flag);

end

