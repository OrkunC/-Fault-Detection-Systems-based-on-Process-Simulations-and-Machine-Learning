close all
clearvars
clc
tic
row = 0 ;
t.rn=8;
tn=4;
%tn = 2:4
spindle_speed = 6000;
for pN= spindle_speed
    for pa = 1:0.1:10.1 
        for pb = 0.1*t.rn:0.05*t.rn:t.rn
            for pf = 0.01:0.01:0.211
%% Initial Data
row = row+1
% Tool Parameters %

t.n=tn;% No. of Edges %
t.p=[90 90 90 90];           % Pitch Angels %

t.r=[t.rn t.rn t.rn t.rn];       % Tool Radius %
t.h=[30 30 30 30];           % Helix Angle %
phid=[0 90 180 270];
t.rake=5;                    % Rake Angle %
t.L=25;
t.runout=0;

%Process Parametes %

p.N = pN;                   % Spindle Speed %
p.V = pi*p.N*2*t.r*1e-3;      % Cutting Speed %
p.f = pf;                   % Feed per rev per edge %
p.a=pa;                        % Axial Depth of Cut %
p.b=pb;                      % Width of Cut %
p.mode=1;                     % Up Milling = 1 , Down Milling = 2 %


%Simulation Parameters%
s.nz=100;                     % # of Axial Elements  %
s.nr=1;                      % # of tool revolution
s.nphi=360;                   % Number of Angular Increments in each revolution %
s.phi=360/s.nphi;             % Angular Increment Deg %
dz=t.L/(s.nz);                   % Element Thickness
%% Failures

tb.E=0;
tr.E=1;

%% Tool Geometry
if tr.E==1
    for i=1:t.n
        t.r(i)=t.r(i)+t.runout*cos(phid(i));
    end
end
if tb.E==0
    phi=deg2rad(phid);               % Tool Edges Positions in the Tip of the Tool (Rad) %%
    tn=t.n;
    
    for k=1:t.n
        j=0;
        for i=0:dz:t.L
            j=j+1;
            loc.lag(j,k)=rad2deg( phi(k)-i*tand(t.h(k))/t.r(k));  %x cordinate of edge points before considering pitch angle in LP
            
        end
    end
    t.pu=t.p;
elseif tb.E==1
    tb.cbt=0;  %Count of failed teeth
    tb.nbt=[]; %position of failed teeth
    tn=t.n-tb.cbt;  %number of existing teeth
    phid(:,tb.nbt)=[];  %position of teeth of broken tool at tool tip
    phi=deg2rad(phid);               % position of teeth of broken tool at tool tip(Rad) %%
    
    
    for k=1:tn
        j=0;
        for i=0:dz:t.L
            j=j+1;
            loc.lag(j,k)=rad2deg( phi(k)-i*tand(t.h(k))/t.r(k));  %lag angle of the exisiting teeth
            
        end
    end
    
    for k=1:tn
        if k==tn
            t.pu(k)=phid(1)-phid(k)+360;      %updated pitch angles between exiiting teeth
        else
            t.pu(k)= phid(k+1)-phid(k);
        end
    end
end



%% Force Calculation

if p.mode==1    %% Entry and Exit angle calculation for Down Milling or Up Milling %%
    phi_st=0;
    phi_ex=acosd(1-p.b/t.rn);
else
    phi_st=180-acosd(1-p.b/t.rn);
    phi_ex=180;
end


s.nza=p.a/dz;

beta=zeros(s.nphi*s.nr,s.nz,t.n);
PHI=zeros(s.nphi*s.nr,s.nz,t.n);
Ktc=zeros(s.nphi*s.nr,s.nz,t.n);
Krc=zeros(s.nphi*s.nr,s.nz,t.n);
Kac=zeros(s.nphi*s.nr,s.nz,t.n);
h=zeros(s.nphi*s.nr,s.nz,t.n);
dFt=zeros(s.nphi*s.nr,s.nz,t.n);
dFr=zeros(s.nphi*s.nr,s.nz,t.n);
dFa=zeros(s.nphi*s.nr,s.nz,t.n);
dFx=zeros(s.nphi*s.nr,s.nz,t.n);
dFy=zeros(s.nphi*s.nr,s.nz,t.n);
dFz=zeros(s.nphi*s.nr,s.nz,t.n);
%AL7075
% tau=297.1+t.rake*1.1;         %shear stress
%Ti6ALV7
tau=613;

% Kte=23.4;    %Edge Force coefficent
% Kre=35.2;
% Kae=5;
Kte=24;    %Edge Force coefficent
Kre=43;
Kae=5;
for i=1:s.nza    %differential forces for each edge on esch elemnt at every rotation angle
    for n=1:tn
        for fi=1:s.nphi*s.nr
          
            im=mod(loc.lag(i,n)+fi,360);
            if im<phi_st || im>phi_ex
                h=0;
            else
                if n==1
                h=t.pu(n)/360*t.n*(t.r(n)-t.r(end)+p.f)*sind(im);  %Chip Thickness
                else
                h=t.pu(n)/360*t.n*(t.r(n)-t.r(n-1)+p.f)*sind(im);  %Chip Thickness
                end
            end
            
            %             beta(fi+1,i,n)=18.8+6.7*h+0.0076*p.V+0.26*t.rake; % Friction Angle for AL7075%
            %             PHI(fi+1,i,n)=24.2+36.7*h+0.005*p.V+0.3*t.rake;   % Shear Angle for AL7075
            beta=22.58;
            PHI=30;
            c=sqrt(cosd(PHI+beta-t.rake)^2+tand(t.h(n))^2*sind(beta)^2);
            Ktc=tau*(cosd(beta-t.rake)+tand(t.h(n))^2*sind(beta))/(sind(PHI)*c);
            
            Krc=tau*(sind(beta-t.rake))/(cosd(t.h(n))*sind(PHI)*c);
            Kac=tau*(cosd(beta-t.rake)*tand(t.h(n))-tand(t.h(n))*sind(beta))/(sind(PHI)*c);
            
            g=1;
            if h==0
                g=0;
            end
            
            dFt(fi,i,n)=g*(Ktc*h+Kte)*dz;
            dFr(fi,i,n)=g*(Krc*h+Kre)*dz;
            dFa(fi,i,n)=g*(Kac*h+Kae)*dz;
            dFx(fi,i,n)=-dFr(fi,i,n)*sind(im)-dFt(fi,i,n)*cosd(im);
            dFy(fi,i,n)=-dFr(fi,i,n)*cosd(im)+dFt(fi,i,n)*sind(im);
            dFz(fi,i,n)=-dFa(fi,i,n);
            
        end
    end
end
Fx=zeros(s.nphi*s.nr,1);
Fy=zeros(s.nphi*s.nr,1);
Fz=zeros(s.nphi*s.nr,1);


for fi=1:s.nphi*s.nr    % Total forces in x y z directions
    Fx(fi)=0;
    Fy(fi)=0;
    Fz(fi)=0;
    for i=1:s.nz
        for n=1:t.n
            Fx(fi)=Fx(fi)+dFx(fi,i,n);
            Fy(fi)=Fy(fi)+dFy(fi,i,n);
            Fz(fi)=Fz(fi)+dFz(fi,i,n);
        end
    end
end
%% Key indicators
Fx_max= max(Fx);
Fx_min= min(Fx);
Fx_maxabs= max (abs (Fx));
Fx_avg= mean (Fx);
Fx_var= max(Fx) - min(Fx);
Fx_max_to_avg= max (abs (Fx))/Fx_avg;
MAX_fxmax_fxmin = max (abs (Fx_max),abs (Fx_min));

Fy_max= max(Fy);
Fy_min= min(Fy);
Fy_maxabs= max (abs (Fy));
Fy_avg= mean (Fy);
Fy_var= max(Fy) - min(Fy);
Fy_max_to_avg= max (abs (Fy))/Fy_avg;
MAX_fymax_fymin = max (abs (Fy_max),abs (Fy_min));

Fz_max= max(Fz);
Fz_min= min(Fz);
Fz_maxabs= max (abs (Fz));
Fz_avg= mean (Fz);
Fz_var= max(Fz) - min(Fz);
Fz_max_to_avg= max (abs (Fz))/Fz_avg;
MAX_fzmax_fzmin = max (abs (Fz_max),abs (Fz_min));

xy_max_ratio= Fx_maxabs/Fy_maxabs;
xz_max_ratio= Fx_maxabs/Fz_maxabs;
yz_max_ratio= Fy_maxabs/Fz_maxabs;

xy_avg_ratio= mean (Fx)/ mean (Fy);
xz_avg_ratio= mean (Fx)/ mean (Fz);
yz_avg_ratio= mean (Fy)/ mean (Fz);

%% Write to excel

% results = [Fx_max Fx_min MAX_fxmax_fxmin Fx_maxabs  Fx_var Fx_max_to_avg Fy_max Fy_min MAX_fymax_fymin Fy_maxabs Fy_avg Fy_var Fy_max_to_avg Fz_max Fz_min MAX_fzmax_fzmin Fz_maxabs Fz_avg Fz_var Fz_max_to_avg xy_max_ratio xz_max_ratio yz_max_ratio xy_avg_ratio xz_avg_ratio yz_avg_ratio];
% parameters = [tn t.r(1) t.h(1) t.rake t.runout p.N p.f p.a p.b p.mode];
% output(row,:) = [parameters results];

results = [Fx_maxabs abs(Fx_avg) Fy_maxabs abs(Fy_avg) Fz_maxabs abs(Fz_avg)];
parameters = [p.N p.f p.a p.b ];
output(row,:) = [parameters results];

            end
        end
    end
end
 filename = 'Simdatawithoutspindlevar.xlsx';
 writematrix(output,filename,'Sheet',1);
