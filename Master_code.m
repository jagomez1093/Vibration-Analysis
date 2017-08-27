% MAE 598 Project
% Jimmy Gomez

clc;
close all; 
clear all; 
format compact
load('project19.mat')

%% Section #1
% Determine the natural frequencies and mode shapes of the structure you were assigned when fixed to ground and in the frequency range of interest.
K_matrix = OUTPUT4_rd3('kgg.dat');
M_matrix = OUTPUT4_rd3('mgg.dat');

% The number of DOFs to be eleminated was 54 by 54. This is due to having nine ground nodes and six degrees of freedom
K_matrix(1:54,:) = [];  
K_matrix(:,1:54) = [];
M_matrix(1:54,:) = [];  
M_matrix(:,1:54) = [];
[Modes,W2] = eig(K_matrix, M_matrix);

for a = 1:length(W2)
    W(a) = W2(a,a);
end
[W_a,I] = sort(W,'ascend');
W_2 = W_a; 
W_2(835:1620) = []; 
W_2(1:22) = [];

for b = 1:length(I)   
	Modes(:,b) = Modes(:,I(b));
end

naturalFreq = sqrt(W_2);
frequency = naturalFreq ./ (2 * pi);

k = 1;
for counterFrequency = 1:length(frequency)
  f1 = frequency(counterFrequency);
  if f1 >= 0
    if f1 <= free_range(2,1)
      naturalF(k) = f1;
      fprintf('Frequency = %f Hz\n', naturalF(k))
      k = k + 1;
    end
  end
end

titles = {'Mode shape 1', 'Mode shape 2', 'Mode shape 3', 'Mode shape 4', 'Mode shape 5', 'Mode shape 6',...
          'Mode shape 7', 'Mode shape 8', 'Mode shape 9', 'Mode shape 10', 'Mode shape 11', 'Mode shape 12',...
          'Mode shape 13', 'Mode shape 14', 'Mode shape 15', 'Mode shape 16'};

c = 6;
for counterD = 23:38;
  if counterD == 36
    p = 2;
  else
    p = 3;
  end 
  for counterE = 1:270
    pval(counterE,counterD-22) = Modes(p,counterD);
    p = c+p;
  end
  figure
  plot3(nodes_coord(10:279,2), nodes_coord(10:279,3), pval(:,counterD-22), '*')
  xlabel('x');
  ylabel('y'); 
  zlabel('mode shape');
  title(sprintf(titles{counterD-22}))
end


%% Section #2
%{ Determine the damped free response of a set of points of the structure to given initial conditions. 
   Discuss convergence with respect to the number of modes retained in the analysis. The time step used in determining the response must
   be appropriately chosen.
   2.1 zeta                  
%}
f = damp_ratio_free(:,1);
zeta = damp_ratio_free(:,2);
f_0 = frequency(1:50);
for e = 1:length(f_0)
  fo = f_0(e);
  zDamp(e) = interp1(f, zeta, fo);
  fprintf('Damping = %f\n', zDamp(e))
end

%{ 2.2 active DOF in which initial velocities are imposed
   2.3 active DOF in which displacementes are to be observed
   2.4 active DOF in which velocities are to be observed     
%}
names = {'Vo', 'V', 'X'};
values = [V0(2), nodes_freeX(1,2), nodes_freeV(1,2)];
for y = 1:3
  fprintf('Active DOF %s = %f', names{y}, values(y));
end

%{ 2.5 List the initial velocities of the first 5 modes
   Take x0 as 0 since its an empty vector
   Solve for Qsi and Qci
   Qci will be equal to zero since the initial diplacement vector is empty.
   This in turn makes the position and velocity completely dependent of Qsi.
%}

X0 = 0;
w = sqrt(W_2(1:length(f_0))) .* sqrt(1 - zDamp .^ 2);
w25 = sqrt(W_2(1:length(f_0)));
time1 = linspace(0, time_free, 1900);
velZ = 0;
node = (99-9) * 6-3;
position = zeros(1, length(time1));
velocity = zeros(1, length(time1));

for counterM = 1:5
  wd = w(counterM);
  wo = w25(counterM);
  zet = zDamp(counterM);
  Qci = Modes(:,counterM + 22)' * M_matrix * X0 / (Modes(:,counterM + 22)' * M_matrix * Modes(:,counterM + 22));
  Qsi = Modes(:,counterM + 22)' * M_matrix * V0(1,3) / (wd * Modes(:,counterM + 22)' * M_matrix * Modes(:,counterM + 22));
  for counterT = 1:length(time1)
    time = time1(counterT);
    position(counterT) = position(counterT) + exp(-zet * wo * time) * Modes(node,counterM + 22) * (Qci(node) * cos(wd * time)+Qsi(node)...
                         * sin(wd * time));
    velocity(counterT) = velocity(counterT) + exp(-zet * wo * time) * Modes(node,counterM + 22) * (-zet * wo * (Qci(node) * cos(wd * time)...
                         + Qsi(node) * sin(wd * time)) + wd * (-Qci(node) * sin(wd * time) + Qsi(node) * cos(wd * time)));
  end
  velZ = velocity(counterT) + velZ;
  V(counterM) = velZ;
end

fprintf('velocity of 5 modes = %f\n', V)
modes = 1:50;
nodesOfInterest = nodes_freeX(:,1);
titles = {'Displacement of node 99 degree of freedom 3', 'Displacement of node 91 degree of freedom 3',...
          'Displacement of node 271 degree of freedom 3', 'Displacement of node 279 degree of freedom 3'};
legendary = {'Number of modes=5', 'Number of modes=10', 'Number of modes=20'};

for counterNodes = 1:length(nodesOfInterest)
  nodesI = nodesOfInterest(counterNodes);
  position = zeros(1,length(time1));
  velocity = zeros(1,length(time1));
  figure
  for counterModes = 1:length(modes)
    m = modes(counterModes);
    wd = w(m);
    wo = w25(m);
    ze = zDamp(m);
    nodes = (nodesI-9) * 6-3;
    Qci1 = Modes(:,m + 22)' * M_matrix * X0/(Modes(:,m + 22)' * M_matrix * Modes(:,m + 22));
    Qsi1 = Modes(:,m + 22)' * M_matrix * V0(1,3)/(wd * Modes(:,m + 22)' * M_matrix * Modes(:,m + 22));
    for counterTime = 1:length(time1)
      time = time1(counterTime);
      position(counterTime) = position(counterTime) + exp(-ze * wo * time) * Modes(nodes,m + 22) * (Qci1(nodes) * cos(wd * time) + ...
                              Qsi1(nodes) * sin(wd * time));
      velocity(counterTime) = velocity(counterTime) + exp(-ze * wo * time) * Modes(nodes,m + 22) * (-ze * wo * (Qci1(nodes) * cos(wd * time)...
                             + Qsi1(nodes) * sin(wd * time)) + wd * (-Qci1(nodes) * sin(wd * time) + Qsi1(nodes) * cos(wd * time)));
    end
    velocityEl(:,counterModes,counterNodes) = velocity;
    hold on
    if m == 5
      plot(time1, position, 'k')
    elseif m == 10
      plot(time1, position, '--')
    elseif m == 20
      plot(time1, position, 'b-')
    end
    hold off
    legend(sprintf(legendary{1}), sprintf(legendary{2}), sprintf(legendary{3}))
    xlabel('time');
    ylabel('Displacement');
    title(sprintf(titles{counterNodes}));
  end
end

titles2 = {'Velocity of node 99 degree of freedom 3', 'Velocity of node 91 degree of freedom 3',...
           'Velocity of node 271 degree of freedom 3', 'Velocity of node 279 degree of freedom 3'};
for counterU = 1:4
  figure
  plot(time1, velocityEl(:,5,counterU), 'b-', time1, velocityEl(:,10,counterU), '--', time1, velocityEl(:,20,counterU), 'r')
  xlabel('time');
  ylabel('velocity');
  title(sprintf(titles2{counterU}));
  legend(sprintf(legendary{1}), sprintf(legendary{2}), sprintf(legendary{3}))
end

%% Step 3
names = {'Force cosine', 'Force sine', 'steady state'};
values = [Fc0(1,3), Fs0(1,3), nodes_harm(1,2)];
for y = 1:3
  fprintf('Active DOF %s = %f\n', names{y}, values(y));
end

nodes_h = nodes_harm(:,1);
omega = linspace(0.1, range_harm(2,1) * 2 * pi,3000);
nod = (nodes_h-9) * 6-3;
titles3 = {'Displacement of node 99 degree of freedom 3', 'Displacement of node 91 degree of freedom 3',...
           'Displacement of node 271 degree of freedom 3', 'Displacement of node 279 degree of freedom 3'};

modC = 30;
for counterQ = 1:4  
  nodes = nod(counterQ);
  S = 0;
  C = 0;
  for counterR = 1:modC
    m = Modes(:,counterR + 22)' * M_matrix * Modes(:,counterR + 22);
    k = Modes(:,counterR + 22)' * K_matrix * Modes(:,counterR + 22);
    c = alpha * m + beta * k;
    fc = Modes(91,counterR + 22) * Fc0(1,3);
    fs = Modes(99,counterR + 22) * Fs0(1,3);

    %Results from solving 
    for counterA = 1:length(omega)
      w = omega(counterA);
      q_sj = (fs * (k - m * w^2) + w * c * fc)/((w * c)^2 + (k - m * w^2)^2);  
      q_cj = (fc - w * c * q_sj)/(k - m * w^2)^2;
      C(counterR,counterA) = q_cj * Modes(nodes,counterR + 22);  
      S(counterR,counterA) = q_sj * Modes(nodes,counterR + 22);
    end
  end
  CD = sum(C);
  SD = sum(S);
  Amp = sqrt(CD.^2 + SD.^2);

  Phase = atan2(SD,CD);
  if Phase < 0
    Phase(w) = Phase + 2 * pi;
  end

  figure
  semilogy(omega/(2 * pi),Amp)
  xlabel('Frequency[Hz]');
  ylabel('Amplitude');
  title(sprintf(titles3{counterQ}))

  figure
  semilogy(omega/(2 * pi),Phase)
  xlabel('Frequency[Hz]');
  ylabel('Phase');
  title(sprintf(titles3{counterQ}))
end