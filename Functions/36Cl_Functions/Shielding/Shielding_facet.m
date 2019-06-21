function [ Shf ] = Shielding_facet( Param_site_in,Const_cosmo)
%SHIELDING_FACET Summary of this function goes here
%   Detailed explanation goes here

    % scaling factors

    zmax = 40000; % maximum depth for calculation of shielding factors (cm)
    dz_surf = 10; % calculated every 10 cm in sub-surface (up to 20 m deep);
    dz_deep = 100; % calculated every 100 cm in sub-surface (from 20 m deep to zmax)
    Z_sc = [(0:dz_surf:2000) (2100:dz_deep:zmax)]; % depth vector to calculate the shielding factors accounting for the topography (cm)
    [~,NbSpl] = size(Param_site_in); % number of samples
    
    Shf=cell(1,NbSpl); % initialization
    
    for i=1:NbSpl
        Param_site = Param_site_in{1,i};

    % Spallation 
        [L_tmp,so_tmp] = Scaling_fact(Param_site.rho_rock,Param_site.rho_coll,Const_cosmo.Lambda_f_t, ...
                                        Param_site.dist_from_fault, Param_site.alpha_rad, ...
                                        Param_site.beta_rad, Param_site.gamma_rad,Z_sc);
        Shf{i}.L_sp = L_tmp;
        Shf{i}.so_sp = so_tmp;
    % Muon
        [L_tmp,so_tmp] = Scaling_fact(Param_site.rho_rock,Param_site.rho_coll,Const_cosmo.Lambda_mu_t, ...
                                        Param_site.dist_from_fault, Param_site.alpha_rad, ...
                                        Param_site.beta_rad, Param_site.gamma_rad,Z_sc);
        Shf{i}.L_mu = L_tmp;
        Shf{i}.so_mu = so_tmp;
        
        Shf{i}.Z_sc = Z_sc; % store the depth vector for futur interpolation
        
        %figure(1)
        %plot(Shf{i}.so_sp.*exp(-Z_sc.*Param_site.rho_rock./Shf{i}.L_sp),-Z_sc./100,'--'); hold on
        %plot(Shf{i}.so_mu.*exp(-Z_sc.*Param_site.rho_rock./Shf{i}.L_mu),-Z_sc./100,'--')
    end

end

