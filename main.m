clc
clear all
close all

%T0e computation (desired transformation matrix)

%given configuration at the instance the picture is captured
DegQc = [ 78.74 18.46 0.35 -65.23 1.43 7.11 -1.19 ];
RadQc =  deg2rad( DegQc );       %converting initial angles to radians

%transformation matrix T07 at the instance the picture is captured
T0e = t_matrix(RadQc);

%TeC computation
%given position and orientation from end effector to camera
xC = -2.1 * 0.0254;
yC = 1.72 * 0.0254;
zC = -1.55 * 0.0254;
thetaYc = -68 * pi/180;
RotYc = myrotmat(thetaYc, 'y');
phiXc = pi/2;
RotXc = myrotmat(phiXc, 'x');
CamRot = RotYc * RotXc;
TCe = [CamRot; 0 0 0];
TCe = [TCe, [ xC; yC; zC; 1 ]];
TeC = inv(TCe);

%TCA computation
%given pose matrix of aruco marker wrt camera
x = 0.169432340063091;
y = -0.084399496991377;
z = 0.760080832922170;
phi = -177.1661574999052 * pi/180;
theta = -21.5174231910357 * pi/180;
psi = -12.2087112074803 * pi/180;

%transformation matrix TCA aruco marker wrt camera
eul = [phi theta psi];
TCA = eye(4,4);
TCA(1:3,1:3) = eul2rotm(eul,'XYZ');
TCA(1:3,4) = [x; y; z];

%TAt computation
xA = 0.375 * 0.0254;
yA = -3.125 * 0.0254;
zA = 0;
psiZA = -pi/2;
RotZA = myrotmat(psiZA, 'z');
thetaYA = -pi;
RotYA = myrotmat(thetaYA, 'x');
ArRot = RotZA * RotYA;
TAt = [ArRot; 0 0 0];
TAt = [TAt, [ xA; yA; zA; 1 ]];

%T0t at Qc
T0t = T0e * TeC * TCA * TAt;

%Tet at Qf (final orientation of robot at target t)
xt = 0.85 * 0.0254;
yt = -1.8 * 0.0254;
zt = -1.403 * 0.0254;
psiZt = -22 * pi/180;
RotZt = myrotmat(psiZt, 'z');
Tte = [RotZt; 0 0 0];
Tte = [Tte, [ xt; yt; zt; 1 ]];
Tet = inv(Tte);

%Final transformation matrix or desired matrix
T0eF = T0t * inv(Tet);

%Trajectory Generation

%desired position pd
pd = T0eF(1:3,4);

%desired orientation oriD
phiD = atan2(T0eF(2,3),T0eF(1,3));
thetaD = atan2(sqrt((T0eF(1,3))^2 + (T0eF(2,3))^2),(T0eF(3,3)));
psiD = atan2(T0eF(3,2),-T0eF(3,1));
oriD = [phiD; thetaD; psiD];

steps = 1000000;

%Storing DH parameters
alpha = [-pi/2, pi/2, pi/2, -pi/2, -pi/2, pi/2, 0];
a = zeros(1,7);
d = [0.340, 0, 0.400, 0, 0.400, 0, 0.126];
q = zeros(7, steps);

%given initial configuration for trajectory generation
q_initial = [ 58.2686; 75.3224; 11.7968; 45.9029; -22.1081; -31.2831; -42.3712 ];
q(:,1) = deg2rad(q_initial);    %storing as radians in 1st column of joint angle matrix

%error matrix initialization
e = zeros(6, steps);

%finding Tphi and Ta for analytical jacobian
Tphi = [0 -sin(phiD) cos(phiD)*sin(thetaD);
        0 cos(phiD) sin(phiD)*sin(thetaD);
        1 0 cos(thetaD)];
I = eye(3,3);
O = zeros(3,3);
Ta = [I O; O Tphi];

%assigning gain value
K = 10*eye(6,6);

for i=1:steps
    xe = pose_matrix(q(:,i));               %storing pose matrix at q
    Ja = jacobian(q(:,i), pd, Ta);          %storing analytical jacobain at q
    e(:,i) = [pd; oriD] - xe;               %computing error between desired and actual pose matrix
    pseudo_inv = pinv(Ja);                  %finding pseudo inverse of Jacobian matrix
    qdot = pseudo_inv*K*e(:,i);             %finding qdot with respest to gain K
    q(:,i+1) = q(:,i) + qdot*0.01;          %incrementing q with qdot*0.01
        if  (max(abs(e(:,i))) < 0.00001)    %breaking loop when error < 0.00001  
         final_config = q(:,i);             %final configuration is last set of q
    break;
        end
     for j = 1:7
        for k = 1:steps
            if q(j,k) > pi
             	q(j,k) = q(j,k) - pi;       %constraining the q values above 180 degrees
                else if q(j,k) < -pi       
                q(j,k) = q(j,k) + pi;       %constraining the q values below 180 degrees
                end
            end            
        end
    end
end

%Testing with forward kinematics giving final configuration
test = pose_matrix(final_config)';
disp("Pose Matrix of Final Configuration (x, y, z, phi, theta, psi)");
disp(test);

%Printing outputs
q_final = final_config';
disp("Computed Configuration (q values)");
disp(q_final);

disp("Desired Transformation Matrix (Computed)");
disp(T0eF);

disp("Actual Transformation Matrix");
disp(t_matrix(final_config));

%Trajectory Design
%setting final time t as 5 seconds
t = 0:0.005:5;

%initializing space for joint angles 1 to 7
q1 = zeros(length(t),1);
q2 = zeros(length(t),1);
q3 = zeros(length(t),1);
q4 = zeros(length(t),1);
q5 = zeros(length(t),1);
q6 = zeros(length(t),1);
q7 = zeros(length(t),1);

%initializing space for joint velocities 1 to 7
q1dot = zeros(length(t),1);
q2dot = zeros(length(t),1);
q3dot = zeros(length(t),1);
q4dot = zeros(length(t),1);
q5dot = zeros(length(t),1);
q6dot = zeros(length(t),1);
q7dot = zeros(length(t),1);

%getting a3 and a2 values from a_values function
a = a_values(q_final);

%solving cubic polynomial for joint angles
q1(:,1) = a(1,1)*t.^3 + a(1,2)*t.^2 + 1.0710;
q2(:,1) = a(2,1)*t.^3 + a(2,2)*t.^2 + 1.3146;
q3(:,1) = a(3,1)*t.^3 + a(3,2)*t.^2 + 0.2059;
q4(:,1) = a(4,1)*t.^3 + a(4,2)*t.^2 + 0.8012;
q5(:,1) = a(5,1)*t.^3 + a(5,2)*t.^2 - 0.3859;
q6(:,1) = a(6,1)*t.^3 + a(6,2)*t.^2 - 0.5460;
q7(:,1) = a(7,1)*t.^3 + a(7,2)*t.^2 - 0.7395;

%storing joint angles in degrees
qtra = [q1 q2 q3 q4 q5 q6 q7];
qtra01 = qtra';

%storing joint angles in degrees in text file
fileID = fopen('Topno_SujalAmrit.txt','w');
fprintf(fileID,'%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n',qtra01);
fclose(fileID);

%solving derivative of cubic polynomial for joint velocities
q1dot(:,1) = 3*a(1,1)*t.^2 + 2*a(1,2)*t;
q2dot(:,1) = 3*a(2,1)*t.^2 + 2*a(2,2)*t;
q3dot(:,1) = 3*a(3,1)*t.^2 + 2*a(3,2)*t;
q4dot(:,1) = 3*a(4,1)*t.^2 + 2*a(4,2)*t;
q5dot(:,1) = 3*a(5,1)*t.^2 + 2*a(5,2)*t;
q6dot(:,1) = 3*a(6,1)*t.^2 + 2*a(6,2)*t;
q7dot(:,1) = 3*a(7,1)*t.^2 + 2*a(7,2)*t;

%storing joint velocities in degrees
qdottra = [q1dot q2dot q3dot q4dot q5dot q6dot q7dot];

%checking for joint limit or velocity limit violations
ct = 0;
for m = 1:length(qtra)
    for n = 1:7
        if n==1
            if qtra(m,n) >=170 || qtra(m,n) <=-170 ||...    %joint limit for joint 1
                 qdottra(m,n) >= 98                         %velocity limit for joint 1
                ct = ct + 1;
            end
        elseif n==2
            if qtra(m,n) >=120 || qtra(m,n) <=-120 ||...
                    qdottra(m,n) >= 98 
                ct = ct + 1;
            end
        elseif n==3
            if qtra(m,n) >=170 || qtra(m,n) <=-170 ||...
                    qdottra(m,n) >= 100
                ct = ct + 1;
            end
        elseif n==4
            if qtra(m,n) >=120 || qtra(m,n) <=-120 ||...
                    qdottra(m,n) >= 130
                ct = ct + 1;
            end
        elseif n==5
            if qtra(m,n) >=170 || qtra(m,n) <=-170 ||...
                    qdottra(m,n) >= 140
                ct = ct + 1;
            end
        elseif n==6
            if qtra(m,n) >=120 || qtra(m,n) <=-120 ||...
                    qdottra(m,n) >= 180
                ct = ct + 1;
            end
        elseif n==7
            if qtra(m,n) >=175 || qtra(m,n) <=-175 ||...
                    qdottra(m,n) >= 180
                ct = ct + 1;
            end
        end
    end
end

if ct > 0
    disp('Joint Limit Exceeded : YES');
    else
    disp('Joint Limit Exceeded : NO');
end

%Plotting Angles, Velocities and Accelerations
subplot(2,1,1);
plot(t,qtra)
title('Joint Angles (Rad)')

subplot(2,1,2); 
plot(t,qdottra)
title('Joint Velocities (Rad/5ms)')