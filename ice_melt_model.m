
function [himax, IRR, dr] = ice_melt_model(p1,p2,hml,Qw,Qsw_mean,i_Qsw)

% main calculation routine for melt period model
% goal is to evaluate retreat date (dr) and removal rate (IRR) given varying initial thickness
% value (himax, array of size Nh)

dhi        = 0.05 ; % Thickness resolution (0.05)
ts_per_day = 400  ; % Number of time steps per day (400)

% Ice properties
himin = 0.5 ; himax = 2.5;
h(1,:) = himin:dhi:himax; N_h = numel(h); % Initial (max) ice thickness
A_ini = 0.99     ; A(1,1:N_h) = A_ini;           % Initial ice concentration

% Calendar
Nd  = 150.                              % Number of days
deltat = 86400. / ts_per_day;            % Time step
Nts = fix( Nd / ( deltat / 86400. ) ) ; % Number of time steps

% Physical "constants"
rhoi = 917.      ; % ice density
L = 335000.      ; % Latent heat
kappa_w = 1./23. ; % attenuation length in water
Rw      = 0.58   ; % fraction of radiation absorbed in the surface layer
alpha_w = 0.06   ; % surface albedo

% Solar flux calculation

if ( i_Qsw == 1 ); Qsw_inc(1:Nts) = Qsw_mean ; end % Constant heat flux
if ( i_Qsw == 2 ); % Linear increase in SW flux
   Qsw_min = 50. ; Qsw_max = 250.         ; 
   dQsw = ( Qsw_max - Qsw_min ) / (Nts-1) ; 
   Qsw_inc(1:Nts) = Qsw_min:dQsw:Qsw_max  ;
end 



% Main calculation loop
for jl = 1:N_h  ; % maximum thickness categories
    for i_time = 2:Nts % time steps

       % Net solar absorption; Lengaigne et al CD 2007, eq 1 
       zqfrac           = ( 1. - ( 1 - Rw ) * exp (- kappa_w * hml) );
       Qsw(i_time)      = Qsw_inc(i_time) * ( 1 - alpha_w ) * zqfrac ; 

       % Thickness change
       zfA              = p1 * ( 1. - A(i_time-1,jl) ) / A(i_time-1,jl) + ( 1 - p1 ) * ( 1 - A_ini ) / A_ini ; % SW repartition function 
       dh_dt(i_time,jl) = - ( Qw + Qsw(i_time) * zfA ) / ( rhoi * L ) ; 
       delta_h          = dh_dt(i_time,jl) * deltat ;
       h(i_time,jl)     = max( [ h(i_time-1,jl) + delta_h 0. ]);

       % Concentration change
       if ( h(i_time,jl) > 0. )
          dA_dt(i_time,jl)  = p2 * A(i_time-1,jl) / ( 2.*h(i_time,jl) ) * dh_dt(i_time,jl);
          delta_A           = dA_dt(i_time,jl) * deltat ;
          A(i_time,jl)      = max( [ A(i_time-1,jl) + delta_A  0. ])  ;
       else
          dA_dt(i_time,jl)  = 0. ; %- A(i_time-1,jl) / deltat;
          A(i_time,jl)      = 0.;
       end

    end
end


%% Diagnostics (IRR, himax, dr)

time_axis_days = (deltat:deltat:Nts*deltat)/86400. ;

himax = h(1,:);

for jl = 1:N_h
   % Mean removal rate
   zaddr = find(A(:,jl) > 0. ); zN =  numel(zaddr);
   IRR(jl) = - nanmean(dA_dt(zaddr(2:zN),jl)) * 86400. * 100.; % skip first time step
   
   % Retreat date
   %ztsr  = max( find(h(:,jl)>0) ) ;
   ztsr  = max( find(A(:,jl)>0.15) ) ;
   dr(jl) = time_axis_days(ztsr)  ;
end



end