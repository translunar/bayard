d2r = pi/180;
as2d = 1/3600;
new_data = 1;
%
% Bayard EM-3455-00-005 ST values
%
%st.b = 333e-6;  % bias in radians, 1-sigma, worst case axis
%st.nea = 333e-6; % NEA in radians, 1-sigma, worst case axis
%st.delta = 0.5;  % Star tracker update sampling period
%st.r = st.delta * st.nea^2;
%
% BCT Nano Star Tracker values
%
st_bct.b   = 60 * as2d * d2r;
st_bct.nea = 333e-6; %60 * as2d * d2r;
st_bct.delta = 0.2; % seconds
st_bct.r = st_bct.delta * st_bct.nea^2;

st_on_sm = st_bct;
st_on_sm.b = .1*pi/180;

%
% SDI500a values
%
gyro.sdi500a.bias_stability = 1.0 * as2d * d2r; % radians per second
gyro.sdi500a.random_walk    = 0.02 * 1/60 * d2r; % radians per root second
gyro.sdi500a.scale_factor   = 200e-6; % parts per million
gyro.sdi500a.nonortho       = (0.2 * 1/1000) / d2r; % radians

%
% SDI500b values
%
gyro.sdi500b.bias_stability = 3.0 * as2d * d2r; % radians/second
gyro.sdi500b.random_walk    = 0.02 * 1/60 * d2r; % radians per root second
gyro.sdi500b.scale_factor   = 300e-6; % parts per million
gyro.sdi500b.nonortho       = (0.2 * 1/1000) / d2r; % radians

%
% CEV MIMU values
%
gyro.cev_mimu.bias_stability = 0.05 * as2d * d2r; % radians per second
gyro.cev_mimu.random_walk    = 0.01 * 1/60 * d2r;   % radians per root second
gyro.cev_mimu.scale_factor   = 5/1e6;         % parts per million
gyro.cev_mimu.nonortho       = 5;             % arcseconds
gyro.cev_mimu.nonortho       = gyro.cev_mimu.nonortho/3600 *d2r;             % radians

%
% Bayard EM-3455-00-005 MIMU values
%
gyro.jpl_mimu.bias_stability = 0.05/3 * as2d * d2r; % radians per second
gyro.jpl_mimu.random_walk    = 0.025/3 * 1/60 * d2r;   % radians per root second



%
% Bayard EM-3455-00-005 LN100 values
%
gyro.jpl_ln100.bias_stability = 0.003 * as2d * d2r; % radians per second
gyro.jpl_ln100.random_walk    = 0.0007 * 1/60 * d2r;   % radians per root second

%
% MSL LN200 values
%
gyro.msl_ln200.bias_stability = 3 * as2d * d2r; % radians per second
gyro.msl_ln200.random_walk    = 0.15 * 1/60 * d2r;   % radians per root second

%
% MEMS GYRO values
%
gyro.mems.bias_stability = 10 * as2d * d2r; % radians per second
gyro.mems.random_walk    = 0.1 * 1/60 * d2r;   % radians per root second

if (new_data)
  dt_sep = 0; %20*60;
  time = [0:60:3600*100];
  ndat2 = length(time);
  
  gyro.cev_mimu.error  = zeros(1,ndat2);
  gyro.cev_mimu2.error = zeros(1,ndat2);
  gyro.jpl_mimu.error  = zeros(1,ndat2);
  gyro.jpl_ln100.error = zeros(1,ndat2);
  gyro.msl_ln200.error = zeros(1,ndat2);
  gyro.mems.error      = zeros(1,ndat2);
  gyro.sdi500a.error   = zeros(1,ndat2);
  gyro.sdi500b.error   = zeros(1,ndat2);
  
  
  %for k = 1:ndat
  % dt = dat.time(k);
  
  for k = 1:length(time)
    dt = time(k);
    
    gyro.cev_mimu.error(k)  = bayard_calc(gyro.cev_mimu, st_bct, dt);
    gyro.cev_mimu2.error(k)  = bayard_calc(gyro.cev_mimu, st_on_sm, dt);
    gyro.jpl_mimu.error(k)  = bayard_calc(gyro.jpl_mimu, st_bct, dt);
    gyro.jpl_ln100.error(k) = bayard_calc(gyro.jpl_ln100, st_bct, dt);
    gyro.msl_ln200.error(k) = bayard_calc(gyro.msl_ln200, st_bct, dt);
    gyro.mems.error(k)      = bayard_calc(gyro.mems, st_bct, dt);
    gyro.sdi500a.error(k)   = bayard_calc(gyro.sdi500a, st_bct, dt);
    gyro.sdi500b.error(k)   = bayard_calc(gyro.sdi500b, st_bct, dt);
    
  end
  
  time_aei = dt_sep;
  time_kosi = dt_sep + 3*3600 + 50*60 + 54;
  time_kosi_hold_end = time_kosi + 60;
  time_6dof = time_kosi_hold_end + 24*60 + 23;
  time_6dof_hold_end = time_6dof + 60*10;
  time_berthing = time_6dof_hold_end + 2*60 + 46;
  time_berthing_end = time_berthing + 60*5 % 5 minutes, rather than 10, in reqs
  time_ei_1 = dt_sep;
  time_apo = dt_sep + 710;
  time_ei_2 = dt_sep + 1320;
  time_chute = dt_sep + 3850;
end




figure(1)
  ylim = [1e-2 1e3];
  loglog(time/3600, [gyro.cev_mimu.error; gyro.jpl_mimu.error; ...
		   gyro.jpl_ln100.error; gyro.msl_ln200.error; ...
		   gyro.mems.error; gyro.sdi500a.error; gyro.sdi500b.error]*180/pi, ...
       [time_aei time_aei]/3600,ylim, 'b--',...
       [time_kosi_hold_end time_kosi_hold_end]/3600,ylim, 'b--',...
       [time_berthing_end time_berthing_end]/3600,ylim, 'b--');
axis([1e-2 1e2 ylim])
grid
legend('CEV MIMU','JPL MIMU','JPL LN100','JPL LN200','MEMS','SDI500A','SDI500B');
xlabel('Time from Last ST Update (Hrs)');
ylabel('Att Error (deg) per-axis 1\sigma');
title('INS Attitude Error vs. Time After Last ST Align')


text(1e-2, 1e2, 'Assumes Star Tracker Update at AI-20 minutes');

figure(2)
ylim = [1e-2 1];
loglog(time/3600, [gyro.cev_mimu.error; gyro.jpl_mimu.error; ...
		   gyro.jpl_ln100.error; gyro.msl_ln200.error; ...
		   gyro.mems.error; gyro.sdi500a.error; gyro.sdi500b.error]*180/pi, ...
       [time_aei time_aei]/3600,ylim, 'b--',...
       [time_kosi_hold_end time_kosi_hold_end]/3600,ylim, 'b--',...
       [time_berthing_end time_berthing_end]/3600,ylim, 'b--');
axis([1e-2 1e2 ylim])
grid
legend('CEV MIMU','JPL MIMU','JPL LN100','JPL LN200','MEMS','SDI500A','SDI500B');
xlabel('Time from Last ST Update (Hrs)');
ylabel('Att Error (deg) per-axis 1\sigma');
title('INS Attitude Error vs. Time After Last ST Align')


text(1e-2, .5, 'Assumes Star Tracker Update at AI-20 minutes');

figure(3)
ylim = [1e-2 1];
loglog(time/3600, [gyro.cev_mimu.error; gyro.jpl_mimu.error; ...
		   gyro.jpl_ln100.error; gyro.msl_ln200.error; ...
		   gyro.mems.error; gyro.sdi500a.error; gyro.sdi500b.error]*180/pi*3, ...
       [time_aei time_aei]/3600,ylim, 'b--',...
       [time_kosi_hold_end time_kosi_hold_end]/3600,ylim, 'b--',...
       [time_berthing_end time_berthing_end]/3600,ylim, 'b--');
axis([1e-2 1e2 ylim])
grid
legend('CEV MIMU','JPL MIMU','JPL LN100','JPL LN200','MEMS','SDI500A','SDI500B');
xlabel('Time from Last ST Update (Hrs)');
ylabel('Att Error (deg) per-axis 3\sigma');
title('INS Attitude Error (3\sigma) vs. Time After Last ST Align')


text(1e-2, .5, 'Assumes Star Tracker Update at AI-20 minutes');


figure(4)
xlim = [.1 5];
ylim = [1e-2 1];
loglog(time/3600, [gyro.cev_mimu.error; gyro.jpl_mimu.error; ...
		   gyro.cev_mimu2.error; gyro.sdi500a.error; gyro.sdi500b.error]*180/pi*3, ...
       [time_aei time_aei]/3600,ylim, 'b--',...
       [time_kosi_hold_end time_kosi_hold_end]/3600,ylim, 'b--',...
       [time_berthing_end time_berthing_end]/3600,ylim, 'b--');       
axis([xlim ylim])
grid
legend('CEV MIMU','JPL MIMU','CEV MIMU ST-ON SM','SDI500A','SDI500B');
xlabel('Time from Last ST Update (Hrs)');
ylabel('Att Error (deg) per-axis 3\sigma');
title('INS Attitude Error (3\sigma) vs. Time After Last ST Align')


text(1e-2, .5, 'Assumes Star Tracker Update at AI-20 minutes');


figure(5)
xlim = [0 5];
ylim = [0 3];
plot(time/3600, [gyro.cev_mimu.error; gyro.jpl_mimu.error; ...
		   gyro.cev_mimu2.error; gyro.sdi500a.error; gyro.sdi500b.error]*180/pi*3, ...
       [time_aei time_aei]/3600,ylim, 'b--',...
       [time_kosi_hold_end time_kosi_hold_end]/3600,ylim, 'b--',...
       [time_berthing_end time_berthing_end]/3600,ylim, 'b--');     
axis([xlim ylim])
grid
legend('CEV MIMU','JPL MIMU','CEV MIMU ST-ON SM','SDI500A','SDI500B');
xlabel('Time from Last ST Update (Hrs)');
ylabel('Att Error (deg) per-axis 3\sigma');
title('INS Attitude Error (3\sigma) vs. Time After Last ST Align')
pause


%text(1e-2, .52, 'Assumes Star Tracker Update at AI-20 minutes');
%text(.35,.95,'EI#1');
%text(.55, .85, 'Apo');
%text(.72, .75, 'EI#2');
%text(1.41, .65, 'Chute');

printflag = 1;
if (printflag)
  figure(1);
  print -dpng large_1sig
  figure(2);
  print -dpng zoom_1sig
  figure(3);
  print -dpng zoom_3sig
  figure(4);
  print -dpng st_on_sm_fig
  figure(5);
  print -dpng st_on_sm_linear
end
