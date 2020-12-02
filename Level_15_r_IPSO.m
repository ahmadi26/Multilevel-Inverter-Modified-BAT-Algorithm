
% Fifteen level Inverter
%clc
 
%clear all
 
%close all
 
ma=0.01:0.01:1;
 
for kk=1:length(ma)
 
 
XX=linspace(0,pi/2,8);
 
xmin=XX(1,1:7); % Upper Bound of Variables
 
xmax=XX(1,2:8); % Lower Bound of Variables
 
n=length(xmax);   % Number of Decision Variables
 
N=10*n;              % Population Size (Swarm Size)
 
    maxit=5000;      % Maximum Number of Iterations
 
% Velocity Limits
 
VelMax=0.1*(xmax-xmin);
 
VelMin=-VelMax;
 
 
wmin=0.2;
wmax=0.9;
c1i=2.5;
c2i=0.2;
c1f=0.2;
c2f=2.5;
alfa =1.5;
epsilon=0.001;
Q=-1*(1/maxit)*log(epsilon);
 
for iter=1:maxit
    c1(iter)=((c1f-c1i).*(iter/maxit))+c1i;
    c2(iter)=((c2f-c2i).*(iter/maxit))+c2i;
    w1(iter)=wmax-(wmax-wmin).*(1+tanh(Q*(iter-(maxit/2))))/2;
    w(iter)=wmin+rand.*w1(iter);
end
 
cost_t=0;
 
for i=1:N
    



    particle(i,:)=rand(1,n).*(xmax-xmin)+xmin; % initial Particles
    
    x=particle(i,:);
    
    if x(1)<0.0175
       cost_t=1e10;
    end
    if x(1)>=x(2)
        cost_t=1e10;
    end
    if x(2)>=x(3)
        cost_t=1e10;
    end
    if x(3)>=x(4)
        cost_t=1e10;
    end
    if x(4)>=x(5)
        cost_t=1e10;
    end
    if x(5)>=x(6)
        cost_t=1e10;
    end
    if x(6)>=x(7)
        cost_t=1e10;
    end
    if x(7)>((pi/2)-0.0175)
        cost_t=1e10;
    end
    if x(2)-x(1)<0.0175
        cost_t=1e10;
    end
    if x(3)-x(2)<0.0175
        cost_t=1e10;
    end
    if x(4)-x(3)<0.0175
        cost_t=1e10;
    end
    if x(5)-x(4)<0.0175
        cost_t=1e10;
    end
    if x(6)-x(5)<0.0175
        cost_t=1e10;
    end
    if x(7)-x(6)<0.0175
        cost_t=1e10;
    end
    
    fitness(i,1)=cost_t+(100*((7*ma(kk))-(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7))))/(7*ma(kk))).^4 +...
                (1/5)*((50*(1/5)*(cos(5*x(1))+cos(5*x(2))+cos(5*x(3))+cos(5*x(4))+cos(5*x(5))+cos(5*x(6))+cos(5*x(7)))/(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7)))).^2)+...
                (1/7)*((50*(1/7)*(cos(7*x(1))+cos(7*x(2))+cos(7*x(3))+cos(7*x(4))+cos(7*x(5))+cos(7*x(6))+cos(7*x(7)))/(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7)))).^2)+...
                (1/11)*((50*(1/11)*(cos(11*x(1))+cos(11*x(2))+cos(11*x(3))+cos(11*x(4))+cos(11*x(5))+cos(11*x(6))+cos(11*x(7)))/(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7)))).^2)+...
                (1/13)*((50*(1/13)*(cos(13*x(1))+cos(13*x(2))+cos(13*x(3))+cos(13*x(4))+cos(13*x(5))+cos(13*x(6))+cos(13*x(7)))/(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7)))).^2)+...
                (1/17)*((50*(1/17)*(cos(17*x(1))+cos(17*x(2))+cos(17*x(3))+cos(17*x(4))+cos(17*x(5))+cos(17*x(6))+cos(17*x(7)))/(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7)))).^2)+...
                (1/19)*((50*(1/19)*(cos(19*x(1))+cos(19*x(2))+cos(19*x(3))+cos(19*x(4))+cos(19*x(5))+cos(19*x(6))+cos(19*x(7)))/(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7)))).^2);
    
    cost_t=0;
             
    vel(i,:)=rand(1,n).*(VelMax-VelMin)+VelMin; % initial Velocity
 
end
 
    %  Personal Best
 
    Pbest=particle;  
    
    fitness_Pbest=fitness;
    
    %  Global Best
    
    particle_mix=[particle,fitness];
 
    particle_sort=sortrows(particle_mix,n+1);
 
    Gbest=particle_sort(1,1:n);  
    
    fitness_Gbest=particle_sort(1,n+1:end);
    
    xmin=[  0   0    0     0    0    0     0 ];
    xmax=[pi/2 pi/2 pi/2 pi/2 pi/2 pi/2  pi/2];
 
for iter=1:maxit
    
    for i=1:N
        
        % Update Velocity
        vel(i,:)=w(iter)*vel(i,:)+...
                 c1(iter)*rand(1,n).*(Pbest(i,:)-particle(i,:))+...
                 c2(iter)*rand(1,n).*(Gbest(1,:)-particle(i,:));
                 
        % Apply Velocity Limits
        vel(i,:)=min([VelMax;vel(i,:)]);
        vel(i,:)=max([VelMin;vel(i,:)]);
        
        
        % Update Position
        particle_new(i,:)=particle(i,:)+vel(i,:);
        
        

        % Apply Position Limits
        particle_new(i,:)=min([xmax;particle_new(i,:)]);
        particle_new(i,:)=max([xmin;particle_new(i,:)]);
        
        % Evaluation
 
        x=particle_new(i,:);
 
        if x(1)<0.0175
           cost_t=1e10;
        end
        if x(1)>=x(2)
            cost_t=1e10;
        end
        if x(2)>=x(3)
            cost_t=1e10;
        end
        if x(3)>=x(4)
            cost_t=1e10;
        end
        if x(4)>=x(5)
            cost_t=1e10;
        end
        if x(5)>=x(6)
            cost_t=1e10;
        end
        if x(6)>=x(7)
            cost_t=1e10;
        end
        if x(7)>((pi/2)-0.0175)
            cost_t=1e10;
        end
        if x(2)-x(1)<0.0175
            cost_t=1e10;
        end
        if x(3)-x(2)<0.0175
            cost_t=1e10;
        end
        if x(4)-x(3)<0.0175
            cost_t=1e10;
        end
        if x(5)-x(4)<0.0175
            cost_t=1e10;
        end
        if x(6)-x(5)<0.0175
            cost_t=1e10;
        end
        if x(7)-x(6)<0.0175
            cost_t=1e10;
        end
       
        
        fitness_new(i,1)=cost_t+(100*((7*ma(kk))-(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7))))/(7*ma(kk))).^4 +...
                         (1/5)*((50*(1/5)*(cos(5*x(1))+cos(5*x(2))+cos(5*x(3))+cos(5*x(4))+cos(5*x(5))+cos(5*x(6))+cos(5*x(7)))/(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7)))).^2)+...
                         (1/7)*((50*(1/7)*(cos(7*x(1))+cos(7*x(2))+cos(7*x(3))+cos(7*x(4))+cos(7*x(5))+cos(7*x(6))+cos(7*x(7)))/(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7)))).^2)+...
                         (1/11)*((50*(1/11)*(cos(11*x(1))+cos(11*x(2))+cos(11*x(3))+cos(11*x(4))+cos(11*x(5))+cos(11*x(6))+cos(11*x(7)))/(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7)))).^2)+...
                         (1/13)*((50*(1/13)*(cos(13*x(1))+cos(13*x(2))+cos(13*x(3))+cos(13*x(4))+cos(13*x(5))+cos(13*x(6))+cos(13*x(7)))/(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7)))).^2)+...
                         (1/17)*((50*(1/17)*(cos(17*x(1))+cos(17*x(2))+cos(17*x(3))+cos(17*x(4))+cos(17*x(5))+cos(17*x(6))+cos(17*x(7)))/(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7)))).^2)+...
                         (1/19)*((50*(1/19)*(cos(19*x(1))+cos(19*x(2))+cos(19*x(3))+cos(19*x(4))+cos(19*x(5))+cos(19*x(6))+cos(19*x(7)))/(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7)))).^2);
    % Fitness Function

        cost_t=0;
        
        if fitness_new(i,1)<fitness(i,1)
            
            delta=zeros(1,n);
            
            for r=1:n
            
                if particle_new(i,r) > particle(i,r)
                      delta(1,r)=1;
                      
                elseif particle_new(i,r) < particle(i,r)
                      delta(1,r)=-1;
                end
            end
                
            particle(i,:)=particle(i,:)+alfa*delta.*abs(vel(i,:));
            
            
            % Apply Position Limits
            particle(i,:)=min([xmax;particle(i,:)]);
            particle(i,:)=max([xmin;particle(i,:)]);
            
            % Evaluation
            x=particle(i,:);
 
            if x(1)<0.0175
               cost_t=1e10;
            end
            if x(1)>=x(2)
                cost_t=1e10;
            end
            if x(2)>=x(3)
                cost_t=1e10;
            end
            if x(3)>=x(4)
                cost_t=1e10;
            end
            if x(4)>=x(5)
                cost_t=1e10;
            end
            if x(5)>=x(6)
                cost_t=1e10;
            end
            if x(6)>=x(7)
                cost_t=1e10;
            end
            if x(7)>((pi/2)-0.0175)
                cost_t=1e10;
            end
            if x(2)-x(1)<0.0175
                cost_t=1e10;
            end
            if x(3)-x(2)<0.0175
                cost_t=1e10;
            end
            if x(4)-x(3)<0.0175
                cost_t=1e10;
            end
            if x(5)-x(4)<0.0175
                cost_t=1e10;
            end
            if x(6)-x(5)<0.0175
                cost_t=1e10;
            end
            if x(7)-x(6)<0.0175
                cost_t=1e10;
            end
 
            fitness(i,1)=cost_t+(100*((7*ma(kk))-(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7))))/(7*ma(kk))).^4 +...
                         (1/5)*((50*(1/5)*(cos(5*x(1))+cos(5*x(2))+cos(5*x(3))+cos(5*x(4))+cos(5*x(5))+cos(5*x(6))+cos(5*x(7)))/(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7)))).^2)+...
                         (1/7)*((50*(1/7)*(cos(7*x(1))+cos(7*x(2))+cos(7*x(3))+cos(7*x(4))+cos(7*x(5))+cos(7*x(6))+cos(7*x(7)))/(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7)))).^2)+...
                         (1/11)*((50*(1/11)*(cos(11*x(1))+cos(11*x(2))+cos(11*x(3))+cos(11*x(4))+cos(11*x(5))+cos(11*x(6))+cos(11*x(7)))/(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7)))).^2)+...
                         (1/13)*((50*(1/13)*(cos(13*x(1))+cos(13*x(2))+cos(13*x(3))+cos(13*x(4))+cos(13*x(5))+cos(13*x(6))+cos(13*x(7)))/(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7)))).^2)+...
                         (1/17)*((50*(1/17)*(cos(17*x(1))+cos(17*x(2))+cos(17*x(3))+cos(17*x(4))+cos(17*x(5))+cos(17*x(6))+cos(17*x(7)))/(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7)))).^2)+...
                         (1/19)*((50*(1/19)*(cos(19*x(1))+cos(19*x(2))+cos(19*x(3))+cos(19*x(4))+cos(19*x(5))+cos(19*x(6))+cos(19*x(7)))/(cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4))+cos(x(5))+cos(x(6))+cos(x(7)))).^2);
    
            cost_t=0;
 
            else
            
            particle(i,:)=particle_new(i,:);
            
            fitness(i,1)=fitness_new(i,1);
        end
 
        % Update Personal Best
        if fitness(i,1)<fitness_Pbest(i,1)
           Pbest(i,:)=particle(i,:);
           fitness_Pbest(i,1)= fitness(i,1);
        end
        
        % Update Global Best
        if fitness(i,1)<fitness_Gbest(1,1)
           Gbest(1,:)=particle(i,:);
           fitness_Gbest(1,1)= fitness(i,1);
        end
 
    end 
    
    final(iter,:)=[Gbest fitness_Gbest];
 
end
 
%% Results
 
V1star=7*ma(kk);
V1=cos(Gbest(1))+cos(Gbest(2))+cos(Gbest(3))+cos(Gbest(4))+cos(Gbest(5))+cos(Gbest(6))+cos(Gbest(7));
V1nesbi=V1/V1star;
V5=abs((1/5)*(cos(5*Gbest(1))+cos(5*Gbest(2))+cos(5*Gbest(3))+cos(5*Gbest(4))+cos(5*Gbest(5))+cos(5*Gbest(6))+cos(5*Gbest(7)))/V1);
V7=abs((1/7)*(cos(7*Gbest(1))+cos(7*Gbest(2))+cos(7*Gbest(3))+cos(7*Gbest(4))+cos(7*Gbest(5))+cos(7*Gbest(6))+cos(7*Gbest(7)))/V1);
V11=abs((1/11)*(cos(11*Gbest(1))+cos(11*Gbest(2))+cos(11*Gbest(3))+cos(11*Gbest(4))+cos(11*Gbest(5))+cos(11*Gbest(6))+cos(11*Gbest(7)))/V1);
V13=abs((1/13)*(cos(13*Gbest(1))+cos(13*Gbest(2))+cos(13*Gbest(3))+cos(13*Gbest(4))+cos(13*Gbest(5))+cos(13*Gbest(6))+cos(13*Gbest(7)))/V1);
V17=abs((1/17)*(cos(17*Gbest(1))+cos(17*Gbest(2))+cos(17*Gbest(3))+cos(17*Gbest(4))+cos(17*Gbest(5))+cos(17*Gbest(6))+cos(17*Gbest(7)))/V1);
V19=abs((1/19)*(cos(19*Gbest(1))+cos(19*Gbest(2))+cos(19*Gbest(3))+cos(19*Gbest(4))+cos(19*Gbest(5))+cos(19*Gbest(6))+cos(19*Gbest(7)))/V1);
V23=abs((1/23)*(cos(23*Gbest(1))+cos(23*Gbest(2))+cos(23*Gbest(3))+cos(23*Gbest(4))+cos(23*Gbest(5))+cos(23*Gbest(6))+cos(23*Gbest(7)))/V1);
V25=abs((1/25)*(cos(25*Gbest(1))+cos(25*Gbest(2))+cos(25*Gbest(3))+cos(25*Gbest(4))+cos(25*Gbest(5))+cos(25*Gbest(6))+cos(25*Gbest(7)))/V1);
V29=abs((1/29)*(cos(29*Gbest(1))+cos(29*Gbest(2))+cos(29*Gbest(3))+cos(29*Gbest(4))+cos(29*Gbest(5))+cos(29*Gbest(6))+cos(29*Gbest(7)))/V1);
V31=abs((1/31)*(cos(31*Gbest(1))+cos(31*Gbest(2))+cos(31*Gbest(3))+cos(31*Gbest(4))+cos(31*Gbest(5))+cos(31*Gbest(6))+cos(31*Gbest(7)))/V1);
V35=abs((1/35)*(cos(35*Gbest(1))+cos(35*Gbest(2))+cos(35*Gbest(3))+cos(35*Gbest(4))+cos(35*Gbest(5))+cos(35*Gbest(6))+cos(35*Gbest(7)))/V1);
V37=abs((1/37)*(cos(37*Gbest(1))+cos(37*Gbest(2))+cos(37*Gbest(3))+cos(37*Gbest(4))+cos(37*Gbest(5))+cos(37*Gbest(6))+cos(37*Gbest(7)))/V1);
V41=abs((1/41)*(cos(41*Gbest(1))+cos(41*Gbest(2))+cos(41*Gbest(3))+cos(41*Gbest(4))+cos(41*Gbest(5))+cos(41*Gbest(6))+cos(41*Gbest(7)))/V1);
V43=abs((1/43)*(cos(43*Gbest(1))+cos(43*Gbest(2))+cos(43*Gbest(3))+cos(43*Gbest(4))+cos(43*Gbest(5))+cos(43*Gbest(6))+cos(43*Gbest(7)))/V1);
V47=abs((1/47)*(cos(47*Gbest(1))+cos(47*Gbest(2))+cos(47*Gbest(3))+cos(47*Gbest(4))+cos(47*Gbest(5))+cos(47*Gbest(6))+cos(47*Gbest(7)))/V1);
V49=abs((1/49)*(cos(49*Gbest(1))+cos(49*Gbest(2))+cos(49*Gbest(3))+cos(49*Gbest(4))+cos(49*Gbest(5))+cos(49*Gbest(6))+cos(49*Gbest(7)))/V1);
 
THD=100.*sqrt((V5^2)+(V7^2)+(V11^2)+(V13^2)+(V17^2)+(V19^2)+(V23^2)+(V25^2)+(V29^2)+(V31^2)+(V35^2)+(V37^2)+(V41^2)+(V43^2)+(V47^2)+(V49^2));
 
Switching_Angles=Gbest(1,1:7)/pi*180;
    fitness_Gbest;



Objective_Function_Value(kk)=fitness_Gbest;
    
    Fundamental_Harmonic(kk)=V1nesbi*100;

    V5th(kk)=V5;
    V7th(kk)=V7;
    V11th(kk)=V11;
    V13th(kk)=V13;
    V17th(kk)=V17; 
    V19th(kk)=V19;
end
 
 
    figure(6);
    semilogy(ma,Objective_Function_Value,'LineWidth',2)
    xlabel('Modulation Index')
    ylabel('ObjectiveFunction Value')
    grid
    
    figure(7);
    plot(ma,Fundamental_Harmonic,'LineWidth',2)


    xlabel('Modulation Index')
    ylabel('Fundamental Harmonic')
    grid
    
    %figure(3);
   % plot(ma,V5th,'LineWidth',2)
   %hold on

    %plot(ma,V7th,'r','LineWidth',2)
    %hold on
   % plot(ma,V11th,'m','LineWidth',2)
%hold on
%plot(ma,V13th,'g','LineWidth',2)
%hold on
%plot(ma,V17th,'k','LineWidth',2)
%hold on
%plot(ma,V19th,'b','LineWidth',2)

    %legend('5th','7th','11th','13th','17th','19th')
    %xlabel('Modulation Index')
    %ylabel('Harmonics')
   % grid

