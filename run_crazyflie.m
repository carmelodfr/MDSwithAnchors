% Copyright (c) 2020
% Author: Carmelo Di Franco 
% Email: carmelodfr@gmail.com
% This code is licensed under MIT license (see LICENSE.txt for details)
 
clc
clear
close all
 
MDSparam.iter       = 400;
MDSparam.verbose    = 1;
MDSparam.xhistory   = 'off';
MDSparam.rtol       = 0.000001;
MDSparam.atol       = 0.001;



X_anchors = [ -2.03   -2.14   0.17; 
            -1.51   1.8     0.15;
            1.78	1.76 	2.41;
            1.30	-1.94	0.27;
            -1.83	-1.84	2.50;
            -1.56	1.81	2.38; 
            1.82	1.80 	0.15; 
            1.29	-1.83	2.53];

 
nodes = 9;
n_anchors = 8;
dimension = 3;

tracked_nodes = 1;

%%
 
load('log_crazyflie.mat')
twrdatarun = twrdatarun2;
groundTruthrun = groundTruthrun2;


 
Xt = [];
Et = [];
Xtgt = [];
      
    
Xtgts = groundTruthrun{:,2:4};
Tgts = groundTruthrun{:,1};
 


D = zeros(nodes,nodes); 
W = zeros(nodes,nodes);
D(2:end,2:end) = pdist2(X_anchors,X_anchors);
W(2:end,2:end) = ~eye(n_anchors,n_anchors);
D_time = zeros(nodes,nodes);
errd = [];

firstrun = 1;
t = 1;
tlen = size(twrdatarun,1);
tmds = 1;
time_0 = 0;


distances = twrdatarun{:,2};
dlmwrite('distances.csv', distances, 'precision', '%6.5f');  

%%
while    t <= tlen -n_anchors
    
   
        for i = 1 : n_anchors

            time = twrdatarun{t,1};
            distance_meas = twrdatarun{t,2};
            anchor_id = twrdatarun{t,3} +1;
            anchor_index = anchor_id + 1;
 
            [~,index] = min(abs(time - Tgts));  
            distance = norm(X_anchors(anchor_id,:) - Xtgts(index,:));
            errd(anchor_id, tmds) = distance_meas -  distance;
            dtrue(anchor_id,tmds) = distance;
            dmeas(anchor_id,tmds) = distance_meas;

            D(1,anchor_index)  = distance_meas +0.26 ;
            D(anchor_index,1)  = distance_meas +0.26 ;
 
            D_time(1,anchor_index)= time;
            D_time(anchor_index,1)= time;
            W(1,anchor_index) = 1;
            W(anchor_index,1) = 1;
            W(1,anchor_index) = 1;
            W(anchor_index,1) = 1;

            t = t+1;
        end
   
    if (firstrun == 1)
        minS = inf;
        for i = 1 : 10
            X0 = [ rand(tracked_nodes,dimension)*10 ; X_anchors(:,1:dimension) ];
            [X,hist] = smacofAnchors(D,X0,W,n_anchors,MDSparam,0);

            if hist.s(end) < minS
                minS = hist.s(end);
                minX = X;
            end  
        end
            X = minX;
            firstrun = 0;
    else
            [X   ,hist   ] = smacofAnchors(D,X0,W,n_anchors,MDSparam,0);
    end
    [~,index] = min(abs(time - Tgts));    
    X0 = X;
    Xt(:,:,tmds) = X;
    Xtgt(:,:,tmds) = Xtgts(index,:);
    
    
    dhatmeas(:,tmds) = D(2:end,1);
    timemds(:,tmds) = time;
     
    Et(:,:,tmds) = X - Xtgt(:,:,tmds); 
       
    tmds = tmds+1;
    time_0 = time;
end
%%
Err1 = squeeze(mean(sqrt(sum(Et(1:tracked_nodes,:,:).^2,2)),1));
 stringTXT = ['Position error of the crazyflie, mean = ' num2str(mean(Err1)) 'm , std = ' num2str(std(Err1)) ' m'];
 
figure()
plot(timemds/1000,Err1);
legend(stringTXT);
xlabel('seconds [s]');
ylabel('RMSE [m]');
box on 
set(gcf,'color','w');
 
figure()
plot3(squeeze(Xtgt(1,1,:)),squeeze(Xtgt(1,2,:)),squeeze(Xtgt(1,3,:)),'b');
hold on

plot3(X_anchors(:,1),X_anchors(:,2),X_anchors(:,3),'ko');
plot3(squeeze(Xt(1,1,:)),squeeze(Xt(1,2,:)),squeeze(Xt(1,3,:)),'rx');

plot3(squeeze(Xt(1,1,:)),squeeze(Xt(1,2,:)),squeeze(Xt(1,3,:)),'r-');

grid on 
legend('Ground truth','Anchors positions','estimated position of the crazyflie')
axis equal
box on 
set(gcf,'color','w');
xlabel('meter [m]');
ylabel('meter [m]');

