%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%Motion estimation using 4 points algorithm%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%                 vs                       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%      2 points algorithm with IMU         %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%        planar scene : homography         %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name : Tsagkatatis Ioannis

% test 2
% Example with noise, propose a test with different camera positions $(R1 =
% I, T1 = 0)$ with angles of rotation between $\rad{0}$ and $\rad{45}$ and translation of 0 to
% 100 AND white noise in image points of camera 2 between 0 to 1 pixel std (use
% RANSAC functions).
%

close all;
clear all;
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    data generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lowerbound=-300;
upperbound=300;
nbpoints = 50;
zposition = 500;
P=[lowerbound + (upperbound-lowerbound).*rand(nbpoints,1),...
    lowerbound + (upperbound-lowerbound).*rand(nbpoints,1),...
    zposition*ones(nbpoints,1),...
    ones(nbpoints,1)];

% the data belong on a plane  of equation N'X+d=0 d=-zposition N = [0,0,1]'

%camera parameter (camera is calibrated)
f=1; u0 = 0; v0 = 0;
K=[f 0 u0;0 f v0;0 0 1];
K1=[f 0 u0 0;0 f v0 0;0 0 1 0];

%%%%%%%%%%%%%%%%%%%%%%%    camera 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% camera position at time 1
tx1=0; ty1=0; tz1=0;
% rotation angles in degree
roll1=0; pitch1=0; yaw1=0;

% rotation of the camera 1

Rp1=[cos(deg2rad(pitch1)) 0 sin(deg2rad(pitch1));...
    0 1 0;...
    -sin(deg2rad(pitch1)) 0 cos(deg2rad(pitch1))];
Ry1=[cos(deg2rad(yaw1)) -sin(deg2rad(yaw1)) 0;...
    sin(deg2rad(yaw1)) cos(deg2rad(yaw1)) 0;...
    0 0 1];
Rr1=[1 0 0;...
    0 cos(deg2rad(roll1)) -sin(deg2rad(roll1));...
    0 sin(deg2rad(roll1)) cos(deg2rad(roll1))];
R1=Ry1*Rp1*Rr1;

T1 = [tx1,ty1,tz1]';

% camera image 1 :
for i =1 : nbpoints
    P1(i,:) = K1*[R1' , -R1'*T1; 0 0 0 1]*P(i,:)';
end

%%%%%%%%%%%%%%%%%%%%%%%    camera 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% camera position at time 2
tx2=10; ty2=-6; tz2=10;
% rotation angles in degree
roll2=3; pitch2=5; yaw2=5;


% rotation of the camera 2

Rp2=[cos(deg2rad(pitch2)) 0 sin(deg2rad(pitch2));...
    0 1 0;...
    -sin(deg2rad(pitch2)) 0 cos(deg2rad(pitch2))];
Ry2=[cos(deg2rad(yaw2)) -sin(deg2rad(yaw2)) 0;...
    sin(deg2rad(yaw2)) cos(deg2rad(yaw2)) 0;...
    0 0 1];
Rr2=[1 0 0;...
    0 cos(deg2rad(roll2)) -sin(deg2rad(roll2));...
    0 sin(deg2rad(roll2)) cos(deg2rad(roll2))];
R2 =Ry2*Rp2*Rr2;

T2 = [tx2,ty2,tz2]';

% camera image 2 :

for i =1 : nbpoints
    P2(i,:) = K1*[R2' , -R2'*T2; 0 0 0 1]*P(i,:)';
end

% Add white Noise
std1 = 0.5;
P2 = P2 + std1.*rand(size(P2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    Display data    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure
% hold on
% plot3(P(:,1),P(:,2),P(:,3),'*')
% plot3(tx1,ty1,tz1,'g*')
% Rx = R1*[100,0,0]';
% line([tx1,tx1+Rx(1)],[ty1,ty1+Rx(2)],[tz1,tz1+Rx(3)],'Color','b')
% Ry = R1*[0,100,0]';
% line([tx1,tx1+Ry(1)],[ty1,ty1+Ry(2)],[tz1,tz1+Ry(3)],'Color','g')
% Rz = R1*[0,0,100]';
% line([tx1,tx1+Rz(1)],[ty1,ty1+Rz(2)],[tz1,tz1+Rz(3)],'Color','r')
% plot3(tx2,ty2,tz2,'r*')
% Rx = R2*[100,0,0]';
% line([tx2,tx2+Rx(1)],[ty2,ty2+Rx(2)],[tz2,tz2+Rx(3)],'Color','b')
% Ry = R2*[0,100,0]';
% line([tx2,tx2+Ry(1)],[ty2,ty2+Ry(2)],[tz2,tz2+Ry(3)],'Color','g')
% Rz = R2*[0,0,100]';
% line([tx2,tx2+Rz(1)],[ty2,ty2+Rz(2)],[tz2,tz2+Rz(3)],'Color','r')
% axis equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    4 pts algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Theoretical Homography H :
N = [0,0,1]';
d = (N'*T1-zposition); 

% get the rotation and trsnlation between two camera and calculate
% theoreotical homography H
Rt = R2'*R1;
Tt = R2'*(T1-T2);
H = Rt-Tt*(R1'*N)'/d;
H=H/H(3,3);

% Homography estimation
% H4pt=homography2d(P1',P2');
[H4pt, inliers] = ransacfithomography(P1', P2', 0.005);
H4pt=H4pt/H4pt(3,3);


%
% Display results
%
disp('### TEST 2 ##');
disp('The Theoretical Homography');
H
disp('Estimated Homography (4 Point Algorithm)');
H4pt

%
% Error calculation
% 
%4 point error with noise
error = sum(sum((H - H4pt).^2));
disp(['Error(SSD) in 4 point estimation, with noise is : ', num2str(error)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    2 pts algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HYPOTHESIS roll and pitch angles are known :

% Virtual image from camera 2 (Z axis correpsonds to the vertical)
PV1 = (Rp1*Rr1*P1(:,1:3)')';
PV2 = (Rp2*Rr2*P2(:,1:3)')';

% 2 points algorithm
N =[0 0 1]';

% Theoretical Homogrphy calculation
HV2 = Ry2'*Ry1-Ry2'*(T1-T2)/(N'*T1-zposition)*(Ry1'*N)';

% True  homography 
TrueHomography = HV2/HV2(3,3);

% Estimated homography from 2 point emthod
% H = homography2d2Pt(PV1',PV2');
[H2pt, inliers] = ransacfithomography2pt(PV1', PV2', 0.005);
EstimatedH=H/H(3,3);

%
% Display results
%
disp('Theoretical Homography (2 angles known)');
HV2
disp('True Theoretical Homography (2 angles known)');
TrueHomography
disp('Estimated Homography (2+1 Point Algorithm)');
EstimatedH

%
% Error calculation
% 
%2 point error with noise
error = sum(sum((TrueHomography - EstimatedH).^2));
disp(['Error(SSD) in 2 point estimation, with noise is : ', num2str(error)]);
