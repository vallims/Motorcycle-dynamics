    %% =========================================================================
%  Stability Analysis of the Linearized Bicycle Model
%  Based on: Meijaard, Papadopoulos, Ruina & Schwab (2007)
%  "Linearized dynamics equations for the balance and steer of a bicycle:
%   a benchmark and review", Proc. R. Soc. A, 463, 1955-1982.
%
%  Equations of motion (linearised, lateral dynamics):
%    M*q'' + v*C1*q' + (g*K0 + v^2*K2)*q = f
%
%  Generalised coordinates: q = [phi; delta]
%    phi   = lean angle (roll)
%    delta = steer angle
%
%  Reference speed range for eigenvalue analysis: 0 to 10 m/s
%% =========================================================================

clear; clc; close all;

%% =========================================================================
%  1. BENCHMARK PARAMETERS  (Table 1 of Meijaard et al. 2007)
%% =========================================================================

% --- Geometry ---
w   = 1.02;          % [m]   wheelbase
c   = 0.08;          % [m]   trail
lam = pi/10;         % [rad] steer-axis tilt (head angle from vertical)
g   = 9.81;          % [m/s^2] gravitational acceleration

% --- Rear wheel (R) ---
rR   = 0.3;          % [m]   radius
mR   = 2.0;          % [kg]  mass
IRxx = 0.0603;       % [kg m^2] spin inertia (about x)
IRyy = 0.12;         % [kg m^2] roll inertia (about y)

% --- Rear frame + rider (B) ---
xB   =  0.3;         % [m]   CoM x-position
zB   = -0.9;         % [m]   CoM z-position (negative = above ground)
mB   = 85.0;         % [kg]
IBxx =  9.2;         % [kg m^2]
IBxz =  2.4;         % [kg m^2] product of inertia
IByy = 11.0;         % [kg m^2]
IBzz =  2.8;         % [kg m^2]

% --- Front fork + handlebar (H) ---
xH   =  0.9;         % [m]
zH   = -0.7;         % [m]
mH   =  4.0;         % [kg]
IHxx =  0.05892;     % [kg m^2]
IHxz = -0.00756;     % [kg m^2]
IHyy =  0.06;        % [kg m^2]
IHzz =  0.00708;     % [kg m^2]

% --- Front wheel (F) ---
rF   = 0.35;         % [m]
mF   =  3.0;         % [kg]
IFxx =  0.1405;      % [kg m^2]
IFyy =  0.28;        % [kg m^2]
%% =========================================================================
% Derived quantites
%% =========================================================================
ITxx = mR*rR^2 + IRxx ...          % rear wheel
     + mB*zB^2 + IBxx ...          % rear frame+rider
     + mH*zH^2 + IHxx ...          % front fork
     + mF*rF^2 + IFxx;   
ITxz = IBxz + IHxz ...                            % frame & fork contributions
     - mB*xB*zB - mH*xH*zH + mF*w*(rF);        % + mR*0*(-rR) = 0

IRzz = IRxx;
IFzz = IFxx;
ITzz = IRzz + IBzz + IHzz + IFzz + mB*xB^2 + mH*xH^2 + mF*w^2;

% Total mass and system CoM
mT = mR + mB + mH + mF;
xT = (xB*mB + xH*mH + w*mF) / mT;
zT = (-rR*mR + zB*mB + zH*mH - rF*mF) / mT;

% Front assembly (H + F)
mA = mH + mF;
xA = (xH*mH + w*mF) / mA;
zA = (zH*mH - rF*mF) / mA;

IAxx = IHxx + IFxx + mH*(zH - zA)^2 + mF*(-rF - zA)^2;
IAxz = IHxz         - mH*(xH - xA)*(zH - zA) - mF*(w - xA)*(-rF - zA);
IAzz = IHzz + IFxx  + mH*(xH - xA)^2 + mF*(w - xA)^2;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
sL  = sin(lam);
cL  = cos(lam);
tL  = tan(lam);
mu    = c/w * cL;
uA  = (xA - w - c)*cL - zA*sL;


IAlx = -mA*uA*zA + IAxx*sL + IAxz*cL;
IAll = mA*uA^2 + IAxx*sL^2 + 2*IAxz*sL*cL + IAzz*cL^2;
IAlz = mA * uA * xA + IAxz * sL + IAzz * cL;

% Gyroscopic coefficients
SR  = IRyy / rR;    % rear wheel spin angular momentum / v
SF  = IFyy / rF;    % front wheel spin angular momentum / v
ST  = SR + SF;      % total gyroscopic coefficient
SA = mA*uA + mu*mT*xT;

%% =========================================================================
%  M matrix
%% =========================================================================
Mphiphi = ITxx;
Mphidel = IAlx + mu * ITxz;
Mdelphi = Mphidel;
Mdeldel = IAll + 2*mu*IAlz + mu^2*ITzz;

M = [Mphiphi, Mphidel;
    Mdelphi, Mdeldel];

%% =========================================================================
% K0 Matrix
%% =========================================================================

K0phiphi = mT*zT;
K0phidel = -SA;
K0delphi = K0phidel;
K0deldel = -SA*sL;

K0 = [K0phiphi, K0phidel;
    K0delphi, K0deldel];

%% =========================================================================
%K2 Matrix
%% =========================================================================

K2phiphi = 0;
K2phidel = ((ST - mT*zT)/w)*cL;
K2delphi = 0;
K2deldel = ((SA + SF*sL)/w)*cL;

K2 = [K2phiphi, K2phidel;
      K2delphi, K2deldel];

%% =========================================================================
% C1 Matrix
%% =========================================================================

C1phiphi = 0;
C1phidel = mu*ST + SF*cL + (ITxz/w)*cL - mu*mT*zT;
C1delphi = -(mu*ST + SF*cL);
C1deldel = (IAlz/w)*cL + mu*(SA + (ITzz/w)*cL);

C1 = [C1phiphi, C1phidel;
    C1delphi, C1deldel];

%% =========================================================================

fprintf('=== Benchmark Matrices (Meijaard et al. 2007) ===\n');
disp('M  ='); disp(M);
disp('C1 ='); disp(C1);
disp('K0 ='); disp(K0);
disp('K2 ='); disp(K2);

%% =========================================================================
%  4. EIGENVALUE SWEEP OVER SPEED  v = 0 .. 10 m/s
%% =========================================================================

v_vec = 0:0.05:10;   % speed vector [m/s]
nv    = length(v_vec);
nq    = 2;           % DOFs -> 4 first-order states

real_parts = zeros(2*nq, nv);   % real parts of eigenvalues
imag_parts = zeros(2*nq, nv);   % imag parts

for k = 1:nv
    eigs_k = compute_eigs(M, C1, K0, K2, v_vec(k), g);
    % Sort by real part (descending) for consistent plotting
    [~, idx] = sort(real(eigs_k), 'descend');
    real_parts(:,k) = real(eigs_k(idx));
    imag_parts(:,k) = imag(eigs_k(idx));
end

%% =========================================================================
%  5. IDENTIFY STABILITY BOUNDARIES
%% =========================================================================

% Weave speed: lowest speed above which weave mode becomes stable
% Capsize speed: speed above which capsize eigenvalue crosses zero
% (Going unstable again)

% Find weave speed (first v where a complex mode has negative real part)
weave_v = NaN;
for k = 2:nv
    has_complex = any(abs(imag_parts(:,k)) > 0.01);
    all_neg_real = all(real_parts(abs(imag_parts(:,k)) > 0.01, k) < 0);
    if has_complex && all_neg_real
        weave_v = v_vec(k);
        break;
    end
end

% Find capsize speed (highest real eigenvalue crosses from neg to pos... 
% or from neg to zero going upward)
capsize_v = NaN;
for k = 2:nv
    % Capsize: a real eigenvalue that is crossing zero from below
    if max(real_parts(:,k)) > 0 && max(real_parts(:,k-1)) <= 0
        capsize_v = interp1([max(real_parts(:,k-1)), max(real_parts(:,k))], ...
                            [v_vec(k-1), v_vec(k)], 0);
    end
end

% Narrow stable speed range: weave_v to capsize_v
fprintf('\n=== Stability Boundaries ===\n');
fprintf('Weave-stable onset (approx):  v_weave  = %.2f m/s\n', weave_v);
fprintf('Capsize instability (approx): v_capsize= %.2f m/s\n', capsize_v);

%% =========================================================================
%  6. PLOTS
%% =========================================================================

fig_width  = 1200;
fig_height = 900;
figure('Name','Meijaard Bicycle Stability','Position',[100 100 fig_width fig_height],...
       'Color','w');

% Color scheme
col = lines(4);
stable_region_color = [0.85 0.95 0.85];

%% --- Subplot 1: Real parts of eigenvalues ---
ax1 = subplot(2,2,1);
hold on; box on; grid on;

% Shade stable region
if ~isnan(weave_v) && ~isnan(capsize_v)
    patch([weave_v capsize_v capsize_v weave_v], [-10 -10 10 10], ...
          stable_region_color, 'EdgeColor','none','FaceAlpha',0.6);
    xline(weave_v,  '--', 'Color',[0.2 0.6 0.2], 'LineWidth',1.2,...
          'Label','v_{weave}','LabelVerticalAlignment','bottom');
    xline(capsize_v,'--', 'Color',[0.8 0.3 0.1], 'LineWidth',1.2,...
          'Label','v_{capsize}','LabelVerticalAlignment','bottom');
end
yline(0, 'k-', 'LineWidth', 0.8);

for i = 1:2*nq
    plot(v_vec, real_parts(i,:), '-', 'LineWidth', 2, 'Color', col(i,:));
end

ylim([-10 10]);
xlabel('Speed  v  [m/s]', 'FontSize', 12);
ylabel('Re(\lambda)  [1/s]', 'FontSize', 12);
title('Real Parts of Eigenvalues', 'FontSize', 13, 'FontWeight','bold');
legend({'Mode 1','Mode 2','Mode 3','Mode 4'}, 'Location','northwest','FontSize',9);
set(ax1, 'FontSize', 10);

%% --- Subplot 2: Imaginary parts (oscillation frequencies) ---
ax2 = subplot(2,2,2);
hold on; box on; grid on;

for i = 1:2*nq
    plot(v_vec, abs(imag_parts(i,:)), '-', 'LineWidth', 2, 'Color', col(i,:));
end

xlabel('Speed  v  [m/s]', 'FontSize', 12);
ylabel('|Im(\lambda)|  [rad/s]', 'FontSize', 12);
title('Oscillation Frequencies', 'FontSize', 13, 'FontWeight','bold');
legend({'Mode 1','Mode 2','Mode 3','Mode 4'}, 'Location','northwest','FontSize',9);
set(ax2, 'FontSize', 10);

%% --- Subplot 3: Root locus in complex plane ---
ax3 = subplot(2,2,3);
hold on; box on; grid on;

xline(0,'k-','LineWidth',0.8);
yline(0,'k-','LineWidth',0.8);

for i = 1:2*nq
    % Color by speed (gradient)
    scatter(real_parts(i,:), imag_parts(i,:), 8, v_vec, 'filled');
end

cb = colorbar;
cb.Label.String = 'Speed [m/s]';
colormap(ax3, parula);
xlabel('Re(\lambda)  [1/s]', 'FontSize', 12);
ylabel('Im(\lambda)  [rad/s]', 'FontSize', 12);
title('Root Locus (colored by speed)', 'FontSize', 13, 'FontWeight','bold');
xlim([-10 10]);
set(ax3, 'FontSize', 10);

%% --- Subplot 4: Stability map & dominant mode ---
ax4 = subplot(2,2,4);
hold on; box on; grid on;

% Compute max real part (stability indicator)
max_real = max(real_parts, [], 1);
stable_mask = max_real <= 0;

area(v_vec, max_real .* stable_mask,   'FaceColor', [0.6 0.9 0.6], 'EdgeColor','none','FaceAlpha',0.5);
area(v_vec, max_real .* ~stable_mask,  'FaceColor', [1.0 0.7 0.7], 'EdgeColor','none','FaceAlpha',0.5);

plot(v_vec, max_real, 'k-', 'LineWidth', 2.5, 'DisplayName', 'max Re(\lambda)');
yline(0,'k--','LineWidth',1.2);

if ~isnan(weave_v)
    xline(weave_v,  '--', 'Color',[0.1 0.5 0.1], 'LineWidth',1.5,...
          'Label','v_{weave}','LabelVerticalAlignment','bottom');
end
if ~isnan(capsize_v)
    xline(capsize_v,'--', 'Color',[0.8 0.2 0.0], 'LineWidth',1.5,...
          'Label','v_{capsize}','LabelVerticalAlignment','bottom');
end

text(mean([weave_v, capsize_v]), -5, 'STABLE', ...
     'HorizontalAlignment','center','Color',[0.1 0.5 0.1],...
     'FontSize',13,'FontWeight','bold');

ylim([-10 5]);
xlabel('Speed  v  [m/s]', 'FontSize', 12);
ylabel('max Re(\lambda)  [1/s]', 'FontSize', 12);
title('Stability Envelope', 'FontSize', 13, 'FontWeight','bold');
legend({'Stable','Unstable','max Re(\lambda)'},'Location','northeast','FontSize',9);
set(ax4, 'FontSize', 10);

sgtitle('Linearized Bicycle Stability  —  Meijaard et al. (2007) Benchmark', ...
        'FontSize', 15, 'FontWeight','bold');

%% =========================================================================
%  7. SEPARATE FIGURE: TIME-DOMAIN SIMULATION  (linearised, at v = 4.5 m/s)
%% =========================================================================

figure('Name','Time-Domain Response','Position',[150 150 900 500],'Color','w');

v_sim = 4.5;    % m/s — inside the stable corridor
tspan = [0 8];
% Initial condition: 3-degree lean, small steer
phi0   = 3*pi/180;
delta0 = 0;
dphi0  = 0;
ddelta0= 0;
x0 = [phi0; delta0; dphi0; ddelta0];

K  = g*K0 + v_sim^2*K2;
C  = v_sim*C1;
Minv = M \ eye(2);

% Build 1st-order state space: A*x, x = [q; qdot]
A_sys = [zeros(2)    eye(2);
        -Minv*K     -Minv*C];

[t_out, x_out] = ode45(@(t,x) A_sys*x, linspace(tspan(1),tspan(2),500), x0);

subplot(2,1,1);
plot(t_out, x_out(:,1)*180/pi, 'b-','LineWidth',2); hold on;
plot(t_out, x_out(:,2)*180/pi, 'r-','LineWidth',2);
xlabel('Time [s]','FontSize',12);
ylabel('Angle [deg]','FontSize',12);
title(sprintf('Time Response — v = %.1f m/s (stable regime)', v_sim),'FontSize',13,'FontWeight','bold');
legend('\phi (lean)','\delta (steer)','Location','best','FontSize',11);
grid on;

subplot(2,1,2);
plot(t_out, x_out(:,3)*180/pi, 'b-','LineWidth',2); hold on;
plot(t_out, x_out(:,4)*180/pi, 'r-','LineWidth',2);
xlabel('Time [s]','FontSize',12);
ylabel('Angular Rate [deg/s]','FontSize',12);
title('Angular Velocities','FontSize',13,'FontWeight','bold');
legend('\dot{\phi}','\dot{\delta}','Location','best','FontSize',11);
grid on;

sgtitle('Linearised Bicycle — Time-Domain Simulation (ODE45)','FontSize',14,'FontWeight','bold');

%% =========================================================================
%  8. PRINT SUMMARY TABLE
%% =========================================================================

fprintf('\n=== Speed-Dependent Eigenvalue Summary ===\n');
fprintf('%-8s  %-14s  %-14s  %-14s  %-14s\n', 'v [m/s]', 'eig1', 'eig2', 'eig3', 'eig4');
fprintf('%s\n', repmat('-',1,72));
v_print = [0, 1, 2, 3, 4, 4.5, 5, 6, 7, 8, 9, 10];
for vp = v_print
    ev = compute_eigs(M, C1, K0, K2, vp, g);
    [~,idx] = sort(real(ev),'descend');
    ev = ev(idx);
    row = sprintf('%-8.2f', vp);
    for ei = 1:4
        if abs(imag(ev(ei))) < 1e-6
            row = [row, sprintf('  %+7.4f      ', real(ev(ei)))]; %#ok
        else
            row = [row, sprintf('  %+6.3f%+5.3fi', real(ev(ei)), imag(ev(ei)))]; %#ok
        end
    end
    fprintf('%s\n', row);
end

fprintf('\n=== Done. ===\n');

%% =========================================================================
%  LOCAL FUNCTION: compute eigenvalues at given speed
%% =========================================================================
function eigs_out = compute_eigs(M, C1, K0, K2, v, g)
%COMPUTE_EIGS  First-order eigenvalues of linearised bicycle at speed v.
%
%  EOM:  M*q'' + v*C1*q' + (g*K0 + v^2*K2)*q = 0
%  State: x = [q; q']
%  A = [0, I; -M\(g*K0+v^2*K2), -M\(v*C1)]

    n = size(M,1);
    K = g*K0 + v^2*K2;
    C = v*C1;
    Mi = M \ eye(n);

    A = [ zeros(n,n),  eye(n);
         -Mi*K,       -Mi*C ];

    eigs_out = eig(A);
end