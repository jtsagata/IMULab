%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%Motion estimation using 4 points algorithm%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%                 vs                       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%      2 points algorithm with IMU         %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%        planar scene : homography         %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name :
%
%
%
%

close all
clear all

%test 1 : example with a particular data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    data generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lowerbound=-300;
upperbound=300;
nbpoints = 50
zposition = 500
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
tx1=10; ty1=-2; tz1=20;
% rotation angles in degree
roll1=5; pitch1=5; yaw1=30;

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

%%%%%%%%%%%%%%%%%%%%%%%    camera 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%    display data     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% explain why d is expressed like this!



Rt = R2'*R1;
Tt = R2'*(T1-T2);
H = Rt-Tt*(R1'*N)'/d;
H=H/H(3,3)
HP1 = H*P1(1,:)';
HP1 = HP1/HP1(3)
GroundTruth = (P2(1,:)/P2(1,3))'


% Homography estimation
H4pt=homography2d(P1',P2');
H4pt=H4pt/H4pt(3,3)

% Homography decomposition
% solutions = invhomog(H4pt);
% 
% solutions(1).T(:,4)/solutions(1).T(3,4);
% solutions(1).T(1:3,1:3)
% 
% solutions(2).T(:,4)/solutions(2).T(3,4);
% solutions(2).T(1:3,1:3);



% HYPOTHESIS roll and pitch angles are known :

% Virtual image from camera 2 (Z axis correpsonds to the vertical)
PV1 = (Rp1*Rr1*P1(:,1:3)')';
PV2 = (Rp2*Rr2*P2(:,1:3)')';

% 2 points algorithm
N =[0 0 1]';

HV2 = Ry2'*Ry1-Ry2'*(T1-T2)/(N'*T1-zposition)*(Ry1'*N)'

HV2*PV1(2,:)'/norm(HV2*PV1(2,:)')
GroundTruth = (PV2(2,:)/norm(PV2(2,:)))'
TrueHomography = HV2/HV2(3,3)

H = homography2d2Pt(PV1',PV2');

EstimatedH=H/H(3,3)




%%%%%%%%%%%%%%%%Yaw Estimation (see exercice 1)

Yaw= atan2(H(2,1),H(1,1))
YawGroundTruth=deg2rad(yaw1-yaw2)

Azimuth=atan2(H(2,3),H(1,3));
SEalpha=-((sin(Yaw)/H(2,1))-1);
CEalpha=(-sin(Yaw)/H(2,1)*H(1,3))/cos(Azimuth);
Elevation=atan2(SEalpha,CEalpha);
T2E = [cos(Azimuth)*cos(Elevation);sin(Azimuth)*cos(Elevation);sin(Elevation)];


T2t=Ry2'*(T1-T2);
T2t/T2t(3)
T2E/T2E(3)




%test 2 : example with different datas
% propose a test with different camera position
% R1 = I, T1 = 0
% angles of rotation of camera 2 between 0� and 45�
% translation of 0 to 100



%test 3 : example with noise
% propose a test with different camera position
% R1 = I, T1 = 0
% angles of rotation of camera 2 between 0� and 45�
% translation of 0  to 100
% AND white noise in image points of camera 2 between 0 to 1 pixel std
%


%test 4 : example with noise on IMU inforamtion
% propose a test with different camera position
% R1 = I, T1 = 0
% angles of rotation of camera 2 between 0� and 45�
% translation of 0  to 100
% AND white noise in image points of camera 2 between 0 to 1 pixel std
% AND white noise in IMU 2 between 0 to 2�

