function [Lambda_sc,So] = Scaling_fact(rho_rock,rho_coll,Lambda, A, alpha, beta, gamma,Z)


%% geometry of the fault and the facet

alpha = alpha; % slope of the colluvial wedge
beta = beta; % dip of the fault-plane
gamma = gamma; % slope of the upper peri-glacial surface

Z(Z==0) = 0.0001;
%% test scsurf_facet
    
    d_fault = A; %  Horizontal distance from the samples to the fault-plane (cm)
    
    Z_above_coll = A/(cos(gamma)/sin(gamma)-cos(beta)/cos(pi/2-beta))*100; % maximum Z above the colluvial wedge
    %Z_below_coll = zmax_1.*100 - Z_above_coll;

    S_1 = zeros(1,length(Z));
    
    for i = 1:length(Z)     % loop on z
        if(Z(i)<Z_above_coll) % The sample is only shielded by the surface of the facet
            a = scsurf_facet(Z(i),Lambda,beta.*180./pi,gamma.*180./pi,rho_rock);
            S_1(i) = a ;
        else % The sample is shielded by the colluvial wedge + the surface of the facet
            a = scdepth_facet(Z(i)-Z_above_coll,d_fault,Lambda, alpha.*180./pi, beta.*180./pi, gamma.*180./pi, rho_rock,rho_coll);
            S_1(i) = a ;
        end
    end
% %     
So = S_1(1);
Lambda_sc = -Z.*rho_rock./(log(S_1)-log(So));
Lambda_sc(1) = Lambda_sc(2); % avoid inf

% figure(1)
%      plot(Lambda_sc,-Z./100);hold on
% figure(2);
%     plot(S_1,-Z./100);hold on;
%     plot(S_1(1).*exp(-Z.*rho_rock./Lambda_sc),-Z./100,'--')

  
