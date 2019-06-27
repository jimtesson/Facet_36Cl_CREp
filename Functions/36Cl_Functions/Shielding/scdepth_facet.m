function sd = scdepth_facet(Z,A,Lambda,alpha,beta,gamma,rho_rock,rho_coll)

% sd = scdepth(Z,H,Lambda,alpha,beta,gamma,rho_rock,rho_coll)
%
%--------------------------------------------------------------------------
% Schlagenhauf A., Gaudemer Y., Benedetti L., Manighetti I., Palumbo L.,
% Schimmelpfennig I., Finkel R., Pou K.
% G.J.Int., 2010
%-------------------------- ? ---------------------------------------------
%
%--------------------------- scdepth.m ------------------------------------
%
% Calculates the scaling factor sd for the buried samples
% as a function of:
%   Z = depth (cm) measured on the scarp (0 at surface, < 0 underneath),
%   H = height of the scarp (cm),
%   Lambda = the true attenuation length (g.cm-2) (for ex. 208 for neutrons),
%   alpha = colluvium dip (degrees) ;
%   beta = scarp dip (degrees),
%   gamma = dip of upper eroded part of the scarp, above beta (degrees),
%   rho_rock = density (g.cm-3) of the rock,
%   rho_coll = density (g.cm-3) of the colluvium.
%--------------------------------------------------------------------------

Z(Z==0) = 1/1000 ; % prevents NaNs

m = 2.3 ; %  Lal exponent
alpha = alpha*pi/180 ;
beta = beta*pi/180 ;
gamma = gamma*pi/180 ;

[theta,phi] = meshgrid(0:90,0:180) ;
theta = theta*pi/180 ;
phi = phi*pi/180 ;
dphi = pi/180 ;
dtheta = pi/180 ;

% Portion of rock above the vertical of the colluvial wedge top (measured
% along the fault-plane)
    Z_above_coll = A/(cos(gamma)/sin(gamma)-cos(beta)/cos(pi/2-beta)); % maximum Z above the colluvial wedge
    Z_above_coll = Z_above_coll/tan(gamma)/cos(beta)-A/cos(beta);% Z projected along the fault-plane

% Angle from vetical to the top of the colluvial wedge.
beta_topcoll=acos((A-(cos(beta)*Z))/(((tan(beta)*(cos(beta)*Z))^2+(A-(cos(beta)*Z))^2)^0.5));

%% PART 1 : colluvial wedge portion + rock between the sample and the fault plane
%           theta = [0 beta_topcoll], phi = [pi 2*pi]

drock = f3(beta,theta,phi+pi,A); % distance from sample to the fault-plane (= thickness of rock in front of the sample)
coll_thick  = f4(alpha,beta,theta,phi+pi,A,Z); % colluvial wedge thickness view from the sample


if(beta_topcoll<pi/2) % case where the sample is not vertically located below the colluvial wedge, beta_topcoll < pi/2
    
    B = atan(tan(beta_topcoll).*sin(phi)) ; %  apparent dip of beta_topcoll in direction phi, in order to remove point with a theta > beta_topcoll
    
    % 1) Portion where cosmic ray travel across the colluvial wedge AND rock in front of the sample.
    %           phi = [pi 2*pi] , theta = [0 beta_topcoll]
    dcwr = exp( -(rho_rock .* drock + rho_coll .* coll_thick)./Lambda ) ; 
    dcwr = dcwr .* (sin(theta).^m) .*cos(theta) .* (theta < B) ;
    dcwr = dcwr * dphi * dtheta ;
    dcwr = sum(dcwr(:));  
    
    % 2) Portion where cosmic ray comes from the facet surface, and only travel across rock : 
    %           phi = [pi 2*pi] , theta = [beta_topcoll pi/2]
    hz = Z+Z_above_coll;
    hz(hz==0) = 1/10000; % to avoid NaNs
    
    dupcol = f5(gamma,beta,theta,phi+pi,hz) ;
    dupcol = exp(-rho_rock * dupcol / Lambda) ;
    dupcol = dupcol .* (sin(theta).^m).*cos(theta).* (theta >= B);
    dupcol = dupcol * dphi * dtheta ;
    dupcol = sum(dupcol(:)) ;
    
    D1 = dcwr + dupcol;
    
else % case where the sample is vertically located below the colluvial wedge, beta_topcoll > pi/2
    % phi = [pi 2*pi] , theta = [0 pi/2]
    
    B = pi/2; % max theta = pi/2
    
    % 1) Portion where cosmic ray travel across the colluvial wedge AND rock in front of the sample.
    %           phi = [pi 2*pi] , theta = [0 beta_topcoll]
    
    dcwr = exp( -(rho_rock .* drock + rho_coll .* coll_thick)./Lambda ) ; 
    dcwr = dcwr .* (sin(theta).^m) .*cos(theta) .* (theta < B) ;
    dcwr = dcwr * dphi * dtheta ;
    
    D1 = sum(dcwr(:)) ;
end
 

%% Part 2: upslope surface, from theta= pi/2 to the fault-plane
   % Downslope surface : phi = [0 pi] , theta = [beta pi/2]
    hz = Z+Z_above_coll;
    hz(hz==0) = 1/10000; % to avoid NaNs
    
    
    C = atan(tan(gamma).*sin(phi)) ; %  apparent dip of gamma in direction phi, in order to remove point with a theta < gamma
    
    if(beta_topcoll>pi/2) % case where the sample is vertically located below the colluvial wedge, beta_topcoll> pi/2
       
        B = atan(tan(pi-beta_topcoll).*sin(phi)) ; %  apparent dip of beta_topcoll in direction phi, in order to remove point with a theta > beta_topcoll
        
        % part below the colluvial wedge : phi = [0 pi] , theta = [beta_topcoll pi/2]
            drock = f3(beta,theta,phi,A); % distance from sample to the fault-plane (= thickness of rock in front of the sample)
            coll_thick  = f4(alpha,beta,theta,phi,A,Z); % colluvial wedge thickness view from the sample

            % Scaling including colluvial wedge + rock in front of the sample.
            d_upsurf1 = exp( -(rho_rock .* drock + rho_coll .* coll_thick)./Lambda ) ; 
            d_upsurf1 = d_upsurf1 .* (sin(theta).^m) .*cos(theta) .* (theta > B) ;
                    %figure(10);surf(d_upsurf1)
            d_upsurf1 = d_upsurf1 * dphi * dtheta ;
            d_upsurf1 = sum(d_upsurf1(:));
        
        % part below the upper surface : phi = [0 pi] , theta = [gamma beta_topcoll]
            d_upsurf2 = f5(gamma,beta,theta,phi,hz) ;
            d_upsurf2 = exp(-rho_rock * d_upsurf2 /Lambda) ;
            d_upsurf2 = d_upsurf2 .* (sin(theta).^m).*cos(theta).* (theta <= B).* (theta > C);
                    %figure(11);surf(d_upsurf2)
            d_upsurf2 = d_upsurf2*dphi*dtheta ;
            d_upsurf2 = sum(d_upsurf2(:));
        
            D2 = d_upsurf1 + d_upsurf2;
    
    else % case where the sample is not vertically located below the colluvial wedge, there is only the upper surface
        B = pi/2 ; % max of Theta
        
        % part below the upper surface: phi = [0 pi] , theta = [gamma pi/2]
        d_upsurf = f5(gamma,beta,theta,phi,hz) ;
        d_upsurf = exp(-rho_rock*d_upsurf/Lambda) ;
        d_upsurf = d_upsurf.*(sin(theta).^m).*cos(theta).* (theta <= B).* (theta > C);
                %figure(12);surf(d_upsurf)
        d_upsurf = d_upsurf*dphi*dtheta ;
        d_upsurf = sum(d_upsurf(:));
        
        D2 = d_upsurf;
    end
    
sd = (D1 + D2)*(m + 1)/(2*pi) ;

end

function d = f(alpha,beta,theta,phi)
num = sin(beta - alpha) ;
den = sin(theta).*cos(alpha) - sin(alpha).*cos(theta).*sin(phi) ;
d = abs(num./den) ;
end

function d = f2(alpha,beta,theta,phi,A,Z)
% Distance from sample to colluvial surface
% d_rockcoll = (Z.*sin(beta-alpha)+A.*sin(alpha))./(sin(theta).*cos(alpha)-cos(theta).*sin(alpha).*sin(phi));

num = (Z.*sin(beta-alpha)+A.*sin(alpha)) ;
den = (sin(theta).*cos(alpha)-cos(theta).*sin(alpha).*sin(phi)) ;
d = abs(num./den) ;
end

function d = f3(beta,theta,phi,A)
% Distance from sample to fault-plane
% d_rock = (A.*sin(beta))./(sin(teta).*cos(beta)-cos(teta).*sin(beta).*sin(phi));

num = (A.*sin(beta)) ;
den = (sin(theta).*cos(beta)-cos(theta).*sin(beta).*sin(phi)) ;
    if(num==0)
        d = den.*0;
    else
        d = abs(num./den) ;
    end
end

function d = f4(alpha,beta,theta,phi,A,Z)
% Colluvial wedge thickness view from the sample

num1 = (Z.*sin(beta-alpha)+A.*sin(alpha));
den1 = (sin(theta).*cos(alpha)-cos(theta).*sin(alpha).*sin(phi));
den1(den1==0)= 1E-10;

num2 = (A.*sin(beta));
den2 = (sin(theta).*cos(beta)-cos(theta).*sin(beta).*sin(phi));
den2(den2==0)= 1E-10;

d = abs(num1./den1-num2./den2) ;
end

function d = f5(alpha,beta,theta,phi,Z)
    if (beta-alpha==0); alpha=alpha-0.0001; end % to avoid sin(0)
    
    num = Z.*sin(beta - alpha);%+d_fault.*sin(alpha) ;
    den = sin(theta).*cos(alpha) - sin(alpha).*cos(theta).*sin(phi) ;
    d = abs(num./den) ;
end
