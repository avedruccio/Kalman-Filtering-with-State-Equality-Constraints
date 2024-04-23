function   ConstrainedKalmanFiltering

% function ConstrainedKalmanFiltering
% This m-file simulates a vehicle tracking problem.
% The vehicle state is estimated with an extended Kalman filter.
% In addition, with the a priori knowledge that the vehicle is on
% a particular road, the vehicle state is estimated with a 
% constrained Kalman filter.
% The state consists of the north and east position, and the
% north and east velocity of the vehicle.
% The measurement consists of ranges to two transponders.


tf = 300; % final time (seconds)
T = 3; % time step (seconds)

Q = diag([4, 4, 1, 1]); % Process noise covariance (m, m, m/sec, m/sec)
Qsqrt = sqrt(Q);

R = diag([900, 900]); % Measurement noise covariance (m, m)
sigma_r = sqrt(R);

% Measurement noise covariance for perfect measurement formulation.
R1 = diag([900, 900, 0, 0]);
sigma_r1 = sqrt(R1);

theta = pi / 3; % heading angle (measured CCW from east)
tantheta = tan(theta);
sintheta = sin(theta);
costheta = cos(theta);

% Define the initial state x, initial unconstrained filter estimate x_hat,
% and initial constrained Kalman filter estimate x_tilde.
x = [0; 0; 170; 100] ;
x_hat = x; % Unconstrained Kalman filter
x_hat_1 = x; % Kalman filter with perfect measurements
x_tilde = x; % Constrained Kalman filter (W=I)
x_tilde_P = x; % Constrained Kalman filter (W=inv(P))

% Initial estimation error covariance
P = diag([R(1,1), R(2,2), Q(1,1), Q(2,2)]);
% Initial estimation error covariance for perfect measurement formulation
P1 = P;

% AccelDecelFlag is used to simulate the vehicle alternately accelerating and
% decelerating, as if in traffic.
AccelDecelFlag = 1;

% Transponder locations.  The first transponder is located at rn1 meters north
% and re1 meters east.  The second transponder is located at rn2 meters north
% and re2 meters east.
%1e5 = 10000
rn1 = 0;
re1 = 0;
rn2 = 1e5 * tantheta;
re2 = 1e5;

% System matrix.
A = [1 0 T 0
   0 1 0 T
   0 0 1 0
   0 0 0 1];

% State constraint matrix.
D = [1 -tantheta 0 0; 
   0 0 1 -tantheta];

% Initialize arrays for saving data for plotting.
%x = [0; 0; 170; 100] ;
x_array = x;

 % Unconstrained Kalman filter
x_hat_array = [];
Constr_array = [];

% Kalman filter with perfect measurements
x_hat_1_array = [];
Constr_1_array = [];

% Constrained Kalman filter (W=I)
x_tilde_array = [];
Constr_Tilde_array = [];

% Constrained Kalman filter (W=inv(P))
x_tilde_P_array = [];
Constr_Tilde_P_array = [];



% Begin the simulation.
for t = T : T : tf
   % Get the measurement.
   z(1, 1) = (x(1)-rn1)^2 + (x(2)-re1)^2; %vector d landmark1 
   z(2, 1) = (x(1)-rn2)^2 + (x(2)-re2)^2; %vector d landmark2
   %error measurements z 
   MeasErr = sigma_r*randn(size(z));
   z = z + MeasErr;
   % Get the measurement for the perfect measurement formulation.
   z1(1, 1) = z(1, 1);
   z1(2, 1) = z(2, 1);
   z1(3, 1) = 0;
   z1(4, 1) = 0;
   
   %z1 = [	z(1,1)
	%		z(2,2)
	%		0
	%		0	  ]
	
	
   % Set the known input.
   if AccelDecelFlag == 1
      if (x(3) > 30) | (x(4) > 30)
         AccelDecelFlag = -1;
      end
   else
      if (x(3) < 5) | (x(4) < 5)
         AccelDecelFlag = 1;
      end
   end
   u = 1 * AccelDecelFlag;
   % Estimate the heading on the basis of the state estimate.
   headinghat = atan2(x_hat(3), x_hat(4));  
   
   
   % Run the Kalman filter.
   H = [2*(x_hat(1)-rn1) 2*(x_hat(1)-rn2) %0
      2*(x_hat(2)-re1) 2*(x_hat(2)-re2) %0
      0 0 
      0 0]; 
   
   %w(k+1) = K
   K = P * H * inv(H' * P * H + R);
   
   % Run the filter for the perfect measurement formulation.
   H1 = [2*(x_hat_1(1)-rn1) 2*(x_hat_1(1)-rn2) 1 0
      2*(x_hat_1(2)-re1) 2*(x_hat_1(2)-re2) -tantheta 0
      0 0 0 1
      0 0 0 -tantheta];
   if (cond(H1' * P1 * H1 + R1) > 1/eps)
      disp('ill conditioning problem');
      return;
   end
   K1 = P1 * H1 * inv(H1' * P1 * H1 + R1);
   
   h(1) = (x_hat(1)-rn1)^2 + (x_hat(2)-re1)^2;
   h(2) = (x_hat(1)-rn2)^2 + (x_hat(2)-re2)^2;
   
   
   
   x_hat = x_hat + K * (z - h');
   
   x_hat_array = [x_hat_array x_hat];
   
   % Find the constrained Kalman filter estimates.
   
   x_tilde = x_hat - D' * inv(D*D') * D * x_hat;
   x_tilde_array = [x_tilde_array x_tilde];
   
   x_tilde_P = x_hat - P * D' * inv(D*P*D') * D * x_hat;
   x_tilde_P_array = [x_tilde_P_array x_tilde_P];
   
   
   h1(1) = (x_hat_1(1)-rn1)^2 + (x_hat_1(2)-re1)^2;
   h1(2) = (x_hat_1(1)-rn2)^2 + (x_hat_1(2)-re2)^2;
   h1(3) = x_hat_1(1) - tantheta * x_hat_1(2);
   h1(4) = x_hat_1(3) - tantheta * x_hat_1(4);
   
   x_hat_1 = x_hat_1 + K1 * (z1 - h1');
   x_hat_1_array = [x_hat_1_array x_hat_1];
   
	B = [0; 0; T*sin(headinghat); T*cos(headinghat)];
   
   x_hat = A*x_hat + B*u;
   Constr_Err = D * x_hat;
   Constr_array = [Constr_array Constr_Err];
   
	B1 = [0; 0; T*sintheta; T*costheta];
   
   x_hat_1 = A*x_hat_1 + B1*u;
   Constr_1_Err = D * x_hat_1;
   Constr_1_array = [Constr_1_array Constr_1_Err];
   
   x_tilde = A*x_tilde + B*u;
   x_tilde_P = A*x_tilde_P + B*u;
   Constr_Tilde = D * x_tilde;
   Constr_Tilde_array = [Constr_Tilde_array Constr_Tilde];
   Constr_Tilde_P = D * x_tilde_P;
   Constr_Tilde_P_array = [Constr_Tilde_P_array Constr_Tilde_P];
   % Update the state estimation error covariance.
	P = (eye(4) - K * H') * P;
   P = A * P * A' + Q;   
   
   P1 = (eye(4) - K1 * H1') * P1;
   P1 = A * P1 * A' + Q;
   
   % Simulate the system.
   B = [0; 0; T*sin(theta); T*cos(theta)];
   x = A*x + B*u + Qsqrt*randn(size(x));
   % Constrain the vehicle (i.e., the true state) to the straight road.
   if abs(x(1) - tantheta * x(2)) > 2
      x(2) = (x(2) + x(1) * tantheta) / (1 + tantheta^2);
      x(1) = x(2) * tantheta;
   end
   if abs(x(3) - tantheta * x(4)) > 0.2
      x(4) = (x(4) + x(3) * tantheta) / (1 + tantheta^2);
      x(3) = x(4) * tantheta;
   end
   x_array = [x_array x];
end

% Process one more measurement.
z(1, 1) = (x(1)-rn1)^2 + (x(2)-re1)^2;
z(2, 1) = (x(1)-rn2)^2 + (x(2)-re2)^2;

MeasErr = sigma_r*randn(size(z));
z = z + MeasErr;

H = [2*(x_hat(1)-rn1) 2*(x_hat(1)-rn2)
   2*(x_hat(2)-re1) 2*(x_hat(2)-re2) 
   0 0 
   0 0 ];

K = P * H * inv(H' * P * H + R);

h(1) = (x_hat(1)-rn1)^2 + (x_hat(2)-re1)^2;
h(2) = (x_hat(1)-rn2)^2 + (x_hat(2)-re2)^2;


x_hat = x_hat + K * (z - h');
x_hat_array = [x_hat_array x_hat];


x_tilde = x_hat - D' * inv(D*D') * D * x_hat;
x_tilde_array = [x_tilde_array x_tilde];

x_tilde_P = x_hat - P * D' * inv(D*P*D') * D * x_hat;
x_tilde_P_array = [x_tilde_P_array x_tilde_P];


headinghat = atan2(x_hat(3), x_hat(4));  

z1(1, 1) = z(1, 1);
z1(2, 1) = z(2, 1);
z1(3, 1) = 0;
z1(4, 1) = 0;

H1 = [2*(x_hat_1(1)-rn1) 2*(x_hat_1(1)-rn2) 1 0
   2*(x_hat_1(2)-re1) 2*(x_hat_1(2)-re2) -tantheta 0
   0 0 0 1
   0 0 0 -tantheta];

K1 = P1 * H1 * inv(H1' * P1 * H1 + R1);

h1(1) = (x_hat_1(1)-rn1)^2 + (x_hat_1(2)-re1)^2;
h1(2) = (x_hat_1(1)-rn2)^2 + (x_hat_1(2)-re2)^2;
h1(3) = x_hat_1(1) - tantheta * x_hat_1(2);
h1(4) = x_hat_1(3) - tantheta * x_hat_1(4);

x_hat_1 = x_hat_1 + K1 * (z1 - h1');
x_hat_1_array = [x_hat_1_array x_hat_1];

% Compute averages.

Constr = sqrt(Constr_array(1,:).^2 + Constr_array(2,:).^2);
Constr = mean(Constr);
disp(['Average Constraint Error (Unconstrained) = ', num2str(Constr)]);

Constr1 = sqrt(Constr_1_array(1,:).^2 + Constr_1_array(2,:).^2);
Constr1 = mean(Constr1);
disp(['Average Constraint Error (Perfect Meas) = ', num2str(Constr1)]);

Constr_Tilde = sqrt(Constr_Tilde_array(1,:).^2 + Constr_Tilde_array(2,:).^2);
Constr_Tilde = mean(Constr_Tilde);
disp(['Average Constraint Error (W=I) = ', num2str(Constr_Tilde)]);

Constr_Tilde_P = sqrt(Constr_Tilde_P_array(1,:).^2 + Constr_Tilde_P_array(2,:).^2);
Constr_Tilde_P = mean(Constr_Tilde_P);
disp(['Average Constraint Error (W=inv(P)) = ', num2str(Constr_Tilde_P)]);

% Plot data.
close all;
t = 0 : T : tf;
figure;
plot(t, x_array(1, :),'r', t, x_array(2, :),'g');hold on;
plot(t, x_hat_array(1, :), t, x_hat_array(2, :));
title('True Position', 'FontSize', 12);

xlabel('seconds'); ylabel('meters');
legend('North Position TRUE', 'East Position TRUE','North Estimation Unconstrained', 'East Estimation Unconstrained');

figure;
plot(t, x_array(1, :) - x_hat_array(1, :), ...
   t, x_array(2, :) - x_hat_array(2, :));
title('Position Estimation Error (Unconstrained)', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('seconds'); ylabel('meters');
legend('North Estimation Error', 'East Estimation Error');



figure;
plot(t, x_array(1, :) - x_hat_1_array(1, :), ...
   t, x_array(2, :) - x_hat_1_array(2, :));
title('Position Estimation Error (Perfect Measurements)', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('seconds'); ylabel('meters');
legend('North Estimation Error', 'East Estimation Error');


figure;
plot(t, x_array(1, :) - x_tilde_array(1, :), ...
   t, x_array(2, :) - x_tilde_array(2, :));
title('Position Estimation Error (Constrained, W=I)', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('seconds'); ylabel('meters');
legend('North Estimation Error', 'East Estimation Error');

figure;
plot(t, x_array(1, :) - x_tilde_P_array(1, :), ...
   t, x_array(2, :) - x_tilde_P_array(2, :));
title('Position Estimation Error (Constrained, W=inv(P))', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('seconds'); ylabel('meters');
legend('North Estimation Error', 'East Estimation Error');