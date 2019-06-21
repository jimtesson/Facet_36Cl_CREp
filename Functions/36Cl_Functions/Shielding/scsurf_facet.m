function Ss = scsurf_facet(Z,Lambda,beta,gamma,rho_rock)

% Ss = scsurf(Z,H,Lambda,beta,gamma,rho_rock)
% 
%--------------------------------------------------------------------------
% Schlagenhauf A., Gaudemer Y., Benedetti L., Manighetti I., Palumbo L.,
% Schimmelpfennig I., Finkel R., Pou K.
% G.J.Int., 2010
%-------------------------- ? ---------------------------------------------
%
%--------------------------- scsurf.m -------------------------------------
%
% Calculates the scaling factor Ss for the exhumed samples
% as a function of:
%   Z = depth (cm) measured on the scarp (0 at surface, > 0 above),
%   H = height of the sarp (cm),
%   Lambda = the true attenuation length (g.cm-2) (for ex. 208 for neutrons), 
%   beta = scarp dip (degrees),
%   gamma = dip of upper eroded part of the scarp, above beta (degrees),
%   rho_rock = density (g.cm-3) of the rock.
%--------------------------------------------------------------------------

hz=Z; hz(hz==0) = 1/10000; % to avoid NaNs

m = 2.3 ; % Lal exponent
beta = beta*pi/180 ;
gamma = gamma*pi/180 ;

[theta,phi] = meshgrid(0:90,0:180) ;
theta = theta*pi/180 ;
phi = phi*pi/180 ;
dphi = pi/180 ;
dtheta = pi/180 ;

% Downslope surface : phi = [pi 2*pi] , theta = [0 pi/2]

da = f(gamma,beta,theta,phi+pi,hz) ;
da = exp(-rho_rock*da/Lambda) ;
da = da.*(sin(theta).^m).*cos(theta);
da = da*dphi*dtheta ;
da = sum(da(:)) ;

% Upslope surface : between vertical and slope of the surface
%       phi = [0 pi] , theta = [C(phi) pi/2]

C = atan(tan(gamma).*sin(phi)) ; % apparent dip of upper part of the scarp in direction phi

dr = f(gamma,beta,theta,phi,hz) ;
dr = exp(-rho_rock*dr/Lambda) ;
dr = dr.*(sin(theta).^m).*cos(theta).*(theta > C) ;
dr = dr*dphi*dtheta ;
dr = sum(dr(:)) ;

S_rock = (dr+da)*(m + 1)/(2*pi) ;
%S_rock = (dr)*(m + 1)/(2*pi) ;
Ss = S_rock ;

end

function d = f(alpha,beta,theta,phi,Z)
if (beta-alpha==0), alpha=alpha-0.0001; end % to avoid sin(0)
num = Z.*sin(beta - alpha);%+d_fault.*sin(alpha) ;
den = sin(theta).*cos(alpha) - sin(alpha).*cos(theta).*sin(phi) ;
d = abs(num./den) ;
end

