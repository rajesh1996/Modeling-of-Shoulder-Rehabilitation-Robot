% ---------ENPM662 Modeling project----------
% -------REACH Robot-Human Interaction-------
% -----SRIKUMAR MURALIDHARAN (116950572)-----
% ---------RAJESHWAR N S (116921237)---------
clc
clear
% Matlab code part to find out the initial 
% angle conditions for shoulder and elbow
l1 = 0.31806; %Shoulder-Elbow portion's link length
l2 = 0.24966; %Elbow-Wrist portion's link length

% D is the distance between the Human midriff and
% the REACH robot's lower boundary.
D = 0.35; %in m

% Distance between midriff and shoulder joint 
C = 0.22059; %in m

% Expanding on the cosine rule mentioned in the report to find
% the angles given the lengths of the triangle, constituted
% by the links L1, L2 and hypotenuse of the C and D components. 
A = sqrt(C^2 + D^2); % hypotenuse

% angle opposite hypotenuse =>
a = acosd((l1^2 + l2^2 - A^2)/(2*l1*l2)); 

% angle opposite shoulder-elbow part =>
b = acosd((A^2 + l2^2 - l1^2)/(2*A*l2));

% angle opposite elbow-wrist part =>
c = acosd((A^2 + l1^2 - l2^2)/(2*A*l1));

% computing the initial joint angles (Shoulder, Elbow, Wrist)  
% for the gripper position at lower boundary condition
theta_e = 180 - a; %elbow joint angle - initial
theta_s = 180 - (atan(D/C)+c); %shoulder joint angle - initial
theta_w = atan2d(C,D)+b; %wrist joint angle - initial
% ------------------------------------------------
% To verify the joint angles at any other point.
% D1 is the distance between the human mid riff and the robot gripper 
% position. Also, 17cm is the total range of the REACH robot to 
% function along the Y axis. So, for D1 has the following range:
%   D <= D1 <= D+0.17 (D1 in m)
D1=D+0.17; %distance in m

% calculating new hypotenuse
A1 = sqrt(C^2 + D1^2);

% angle opposite hypotenuse =>
a1 = acosd((l1^2 + l2^2 - A1^2)/(2*l1*l2));

% angle opposite shoulder-elbow part =>
b1 = acosd((A1^2 + l2^2 - l1^2)/(2*A1*l2));

% angle opposite elbow-wrist part =>
c1 = acosd((A1^2 + l1^2 - l2^2)/(2*A1*l1));

theta_e1 = 180 - a1;%elbow joint angle for end point
theta_s1 = 180 - (atand(D1/C)+c1);%shoulder joint angle at end point
theta_w1 = atan2d(C,D1)+b1;%wrist joint angle at end point
% We have set the joint angle at the initial pose as the reference
% angles. So, to validate the values obtained from VRep, we have
% to compare the following angles as a difference between obtained
% and initial joint angles.
theta1_e=theta_e-theta_e1 %Elbow joint angle measured wrt initial pose 
theta1_s=theta_s-theta_s1 %Shoulder jt angle measured wrt initial pose
theta1_w=theta_w1-theta_w %Wrist jt angle measured wrt initial pose

%---------------------------------------------------------------------
% From the DH tables created, we can frame the Transformation matrices 
% using the rots function attached below, as follows:
A1_0=rots(theta_s1,0,l1,0); %First row of DH table
A2_0=rots(theta_e1,0,l2,0); %Second row of DH table
A3_0=rots(theta_w1+90,0,0,90); %Third row of DH table
A4_0=rots(0,D1-D,0,0); %Fourth row of DH table

% Corresponding transformation matrices
H1_0=A1_0; %For Shoulder (revolute) joint
H2_0=H1_0*A2_0; %For Elbow (revolute) joint
H3_0=H2_0*A3_0; %For Wrist (revolute) joint
H4_0=H3_0*A4_0; %For gripper (prismatic) actuator

% Extracting R and O matrices to obtain Jacobian matrix 
O0=[0 0 0]'; % Origin coordinates
R1=H1_0(1:3,1:3);
O1=H1_0(1:3,4);
R2=H2_0(1:3,1:3);
O2=H2_0(1:3,4);
R3=H3_0(1:3,1:3);
O3=H3_0(1:3,4);
R4=H4_0(1:3,1:3);
O4=H4_0(1:3,4);

% Calculating the values of Z for Jacobian matrix
Z0=eye(3)*[0 0 1]';
Z1=R1*[0 0 1]';
Z2=R2*[0 0 1]';
Z3=R3*[0 0 1]';

% Computing the Jacobian Matrix componenents and framing J
J11=cross(Z0,(O4-O0)); %Jv component of Shoulder jt
J12=cross(Z1,(O4-O1)); %Jv component of Elbow jt
J13=cross(Z2,(O4-O2)); %Jv component of Wrist jt
J14=Z3; %Jv component of Gripper prismatic jt
J21=Z0; %Jw component of Shoulder jt
J22=Z1; %Jw component of Elbow jt
J23=Z2; %Jw component of Wrist jt
J24=[0 0 0]'; %Jw component of Gripper prismatic jt
J=[J11 J12 J13 J14;J21 J22 J23 J24]; %Jacobian matrix created for the same.

%----------------------------------------------------
% Now, we compute the values for joint torques, using the
% joint angles and joint tortional stiffness values. The torsional 
% stiffness values for the same are given as follows:
K_s=100; %Torsional stiffness for shoulder = 100Nm/rad.
K_e=60; %Torsional stiffness for elbow = 60Nm/rad.
K_w=40; %Torsional stiffness for wrist = 40Nm/rad.
K_gripper=50; %Stiffness for gripper modeled as a spring = 50N/m
T_shoulder=K_s*deg2rad(theta_s1);
T_elbow=K_e*deg2rad(theta_e1);
T_wrist=K_w*deg2rad(theta_w1);
F_gripper=K_gripper*(D1-D);
% Creating the new joint Torque matrix:
T= [T_shoulder;T_elbow;T_wrist;F_gripper];
F=pinv(J')*T;

% Defining the syntax of the function used above.
% ---------------------------------------------------------------------
function T = rots(theta, d, a, alpha)
T=[cosd(theta) -sind(theta)*cosd(alpha) sind(theta)*sind(alpha) a*cosd(theta); sind(theta) cosd(theta)*cosd(alpha) -cosd(theta)*sind(alpha) a*sind(theta); 0 sind(alpha) cosd(alpha) d; 0 0 0 1];
end
% ---------------------------------------------------------------------