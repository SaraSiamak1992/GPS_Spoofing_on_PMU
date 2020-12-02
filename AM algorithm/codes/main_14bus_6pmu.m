clc
clear all
close all

global a Nb Mn Z_at num_meas V Z_true W Z Gam_new

%% 4-bus System  --   1 attack  --   delta_teta_spf=pi/4 % spoofing teta  on PMU 1 ----  run Matrix_H_6pmu.m  
Nb=14 ;
Npmu=6 ;
p=1; % the attacked PMU       
a=[0;1 ;0 ;1; 0; 1; 1; 0; 0; 1; 0; 0;0 ; 1] ; %

%% lines connected to bus n
n1=4  ;n2=5 ;n3=4 ; n4=3 ;  n5=2 ;  n6=2  ; 

%% number of measurements
Mn=[2*(n1+1),2*(n2+1),2*(n3+1),2*(n4+1),2*(n5+1),2*(n6+1)];
for i=1:1:Npmu
    sm=sum(Mn(1:i-1));
    num_meas=sum(Mn)
end
%% covariance matrix (sigma)

sg1=0.02 ;  sg2=0.03 ;sg3=0.03 ;sg4=0.02 ; sg5=0.01 ; sg6=0.02 ; 
Sg=[sg1,sg2,sg3,sg4,sg5,sg6];
sig=zeros(num_meas,1);

for k=1:1:Npmu
    sm=sum(Mn(1:k-1));
    sig(sm+1:sm+Mn(k))=Sg(k)*ones(Mn(k),1) ;
end

Sigma=diag(sig);

%% measurement noise
mu=0;
W=zeros(num_meas,1);
for i=1:1:Npmu
    sm=sum(Mn(1:i-1));
    W(sm+1:sm+Mn(i))=normrnd(mu,Sg(i),Mn(i),1);
end

%% Gama
delta_teta_spf=0 % spoofing teta
gama=[cos(delta_teta_spf),-sin(delta_teta_spf);sin(delta_teta_spf),cos(delta_teta_spf)] ;

 Gama=eye(num_meas,num_meas);
 
for i=1:1:Npmu
    if i == p
    sm=sum(Mn(1:i-1))
    
    o=Mn(i)
    j=1
   while j<o
      
      Gam(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama);
      j=j+2
    end
    Gama(sm+1:sm+Mn(i),sm+1:sm+Mn(i))=Gam;
    end
end
Gama




%% H  % num_meas >< 2Nb

H_matrix_6pmu

S_states_14bus
S
 
V=[real(S(1,1));imag(S(1,1));real(S(2,1));imag(S(2,1));real(S(3,1));imag(S(3,1));real(S(4,1));imag(S(4,1));real(S(5,1));imag(S(5,1));real(S(6,1));imag(S(6,1));real(S(7,1));imag(S(7,1))
    real(S(8,1));imag(S(8,1));real(S(9,1));imag(S(9,1));real(S(10,1));imag(S(10,1));real(S(11,1));imag(S(11,1));real(S(12,1));imag(S(12,1));real(S(13,1));imag(S(13,1));real(S(14,1));imag(S(14,1))];
%% %% phasor measurements
Z_true=H*V

Z = Z_true + W

Z_at=Gama*Z_true+W

Z_attk = zeros(num_meas,1);
for k=1:1:Npmu
     sm=sum(Mn(1:k-1))
    Z_attk(sm+1:sm+Mn(k),1)=Z_at(sm+1:sm+Mn(k),1);
end



%% Estimation (with spoofing)
G1=[zeros(2*Nb,2*Nb),zeros(2*Nb,2*Nb),zeros(2*Nb,2*Nb),zeros(2*Nb,2*Nb),zeros(2*Nb,2*Nb),zeros(2*Nb,2*Nb)];

sig=zeros(num_meas,1);

for k=1:1:Npmu
    sm=sum(Mn(1:k-1));
    sig(sm+1:sm+Mn(k))=Sg(k)*ones(Mn(k),1) ;
end
Sigma=diag(sig);
invers_sigm = inv(Sigma);
k=1;
for i=1:1:Nb
    if a(i)==1
        sm=sum(Mn(1:k-1));
        G1(:,(k-1)*2*Nb+1:(k-1)*2*Nb+2*Nb)=a(i)*(H(sm+1:sm+Mn(k),:))'* invers_sigm(sm+1:sm+Mn(k),sm+1:sm+Mn(k))*H(sm+1:sm+Mn(k),:);
        k=k+1;
    end
end
i=1;
G=zeros(2*Nb,2*Nb)
while i<=Npmu
  
    G=G+G1(:,(i-1)*2*Nb+1:(i-1)*2*Nb+2*Nb)
     i=i+1
end

k=1;
 SuM1=zeros(2*Nb,2);
for i=1:1:Nb
    if a(i)==1
      sm=sum(Mn(1:k-1));

      SuM1(:,k)=a(i)*(H(sm+1:sm+Mn(k),:))'*invers_sigm(sm+1:sm+Mn(k),sm+1:sm+Mn(k))*Z_at(sm+1:sm+Mn(k),1);
       k=k+1;
    end
     
end
i=1;
SuM=zeros(2*Nb,1);
while i<=Npmu
  
    SuM=SuM+SuM1(:,i);
     i=i+1;
end
V_hat_ML_at = (inv(G))*SuM

figure
hold on
x=1:1:2*Nb;
bar(x,[V,V_hat_ML_at])
xlabel('n')
ylabel('V & V_hat_ML_at')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% greedy algortithm
Nb=14
a=[0;1 ;0 ;1; 0; 1; 1; 0; 0; 1; 0; 0;0 ; 1] ;
delta_teta_max=input('plz enter delta_teta_max=') ;
Np=input('please enter the number of PMUs attacked ; Np=') ;

b=zeros(Nb,1)

%b_c=zeros(Nb,comb);
 
%% step1   NP=1
NP1=1;
comb=combntns(Nb,NP1)
b_c=zeros(Nb,comb);
for i=1:1:comb
    b(i)=1;
    b_c(i,i)=1;
end
c=zeros(Nb,comb);
%b_c=zeros(Nb,comb);
L=0;

bb=zeros(Nb,comb);
 for i=1:1:comb
      b_bar=b_c(:,i)
      k=1;
      
            while k<=Nb 
                   if b_bar(k)<=a(k)
                       bb(k,i)=b_bar(k);
                      
                   else
                      k=Nb+1;
                   end
                              if k==Nb
                              L=L+1  ;
                             c(:,i)=bb(:,i) ;  
                              end
            k=k+1;
           end
 end
c(:,1)=[];
c(:,3-1)=[];
c(:,5-2)=[];
c(:,8-3)=[];
c(:,9-4)=[];
c(:,11-5)=[];
c(:,12-6)=[];
c(:,13-7)=[];


b=c
Delta_Teta=zeros(L,1);
%  maxim=zeros(L,3);
 for k=1:1:L
     
     b(:,k)
     o=find(b(:,k)==1)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
          if o==2
  %initialization  
int1=zeros(Nb,1)
int2=-delta_teta_max.*b(:,k)
int3=delta_teta_max.*b(:,k)
int=[int1,int2,int3]
% 
%norm2= norm(B_ML*v,2)
% X0=[int1,int2,int3]
for l=1:1:3

A=[];
B=[]
Aeq=[]
beq=[]
lb1=-delta_teta_max.*b(:,k)
ub1=delta_teta_max.*b(:,k)
%%%%
X0=int(l)
lb=lb1(o)
ub=ub1(o)

%Gam_new=zeros(num_meas,num_meas);
% f1=fun_4bus1(delta_teta_spf1)
delta_teta_spf1=fmincon(@fun_6bus1,X0,A,B,Aeq,beq,lb,ub)

% Delta_Teta(k)=delta_teta_spf1
f1=fun_6bus1(delta_teta_spf1)
F1=-f1
max_f1=max(F1)
end
    Delta_Teta(k)=delta_teta_spf1
          end
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
          if o==4
              
               %initialization  
int1=zeros(Nb,1)
int2=-delta_teta_max.*b(:,k)
int3=delta_teta_max.*b(:,k)
int=[int1,int2,int3]
% 
%norm2= norm(B_ML*v,2)
% X0=[int1,int2,int3]
for l=1:1:3

A=[];
B=[]
Aeq=[]
beq=[]
lb1=-delta_teta_max.*b(:,k)
ub1=delta_teta_max.*b(:,k)
%%%%
X0=int(l)
lb=lb1(o)
ub=ub1(o)
%Gam_new=zeros(num_meas,num_meas);
delta_teta_spf2=fmincon(@fun_6bus2,X0,A,B,Aeq,beq,lb,ub)
% Delta_Teta(k,l)=delta_teta_spf2 
f2=fun_6bus2(delta_teta_spf2)
F2=-f2  
max_f2=max(F2)
end  
 Delta_Teta(k,l)=delta_teta_spf2 
          end
         
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
          if o==6
              
               %initialization  
int1=zeros(Nb,1)
int2=-delta_teta_max.*b(:,k)
int3=delta_teta_max.*b(:,k)
int=[int1,int2,int3]
% 
%norm2= norm(B_ML*v,2)
% X0=[int1,int2,int3]
for l=1:1:3

A=[];
B=[]
Aeq=[]
beq=[]
lb1=-delta_teta_max.*b(:,k)
ub1=delta_teta_max.*b(:,k)
%%%%
X0=int(l)
lb=lb1(o)
ub=ub1(o)
%Gam_new=zeros(num_meas,num_meas);
delta_teta_spf3=fmincon(@fun_6bus3,X0,A,B,Aeq,beq,lb,ub)
% Delta_Teta(k)=delta_teta_spf3 
f3=fun_6bus3(delta_teta_spf3)
F3=-f3  
max_f3=max(F3)
end  
Delta_Teta(k)=delta_teta_spf3 
          end
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
          if o==7
              
               %initialization  
int1=zeros(Nb,1)
int2=-delta_teta_max.*b(:,k)
int3=delta_teta_max.*b(:,k)
int=[int1,int2,int3]
% 
%norm2= norm(B_ML*v,2)
% X0=[int1,int2,int3]
for l=1:1:3

A=[];
B=[]
Aeq=[]
beq=[]
lb1=-delta_teta_max.*b(:,k)
ub1=delta_teta_max.*b(:,k)
%%%%
X0=int(l)
lb=lb1(o)
ub=ub1(o)
%Gam_new=zeros(num_meas,num_meas);
delta_teta_spf4=fmincon(@fun_6bus4,X0,A,B,Aeq,beq,lb,ub)
% Delta_Teta(k,l)=delta_teta_spf4 
f4=fun_6bus4(delta_teta_spf4)
F4=-f4  
max_f4=max(F4)
end  

 Delta_Teta(k,l)=delta_teta_spf4
          end  
          
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
          if o==10 
              
               %initialization  
int1=zeros(Nb,1)
int2=-delta_teta_max.*b(:,k)
int3=delta_teta_max.*b(:,k)
int=[int1,int2,int3]
% 
%norm2= norm(B_ML*v,2)
% X0=[int1,int2,int3]
for l=1:1:3

A=[];
B=[]
Aeq=[]
beq=[]
lb1=-delta_teta_max.*b(:,k)
ub1=delta_teta_max.*b(:,k)
%%%%
X0=int(l)
lb=lb1(o)
ub=ub1(o)
%Gam_new;
delta_teta_spf5=fmincon(@fun_6bus5,X0,A,B,Aeq,beq,lb,ub)
% Delta_Teta(k,l)=delta_teta_spf5 
f5=fun_6bus5(delta_teta_spf5)
F5=-f5  
max_f5=max(F5)
end  
Delta_Teta(k,l)=delta_teta_spf5 
          end  
          
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
          if o==14
              
               %initialization  
int1=zeros(Nb,1)
int2=-delta_teta_max.*b(:,k)
int3=delta_teta_max.*b(:,k)
int=[int1,int2,int3]
% 
%norm2= norm(B_ML*v,2)
% X0=[int1,int2,int3]
for l=1:1:3

A=[];
B=[]
Aeq=[]
beq=[]
lb1=-delta_teta_max.*b(:,k)
ub1=delta_teta_max.*b(:,k)
%%%%
X0=int(l)
lb=lb1(o)
ub=ub1(o)
%Gam_new;
delta_teta_spf6=fmincon(@fun_6bus6,X0,A,B,Aeq,beq,lb,ub)
% Delta_Teta(k,l)=delta_teta_spf6 
f6=fun_6bus6(delta_teta_spf6)
F6=-f6  
max_f6=max(F6)
end  
Delta_Teta(k,l)=delta_teta_spf6 
          end  
    
 end

% compare
%delta_teta_spf=[delta_teta_spf1,delta_teta_spf2,delta_teta_spf3,delta_teta_spf4,delta_teta_spf5,delta_teta_spf6]

delta_teta_spf=[Delta_Teta(1,1),Delta_Teta(2,1),Delta_Teta(3,1),Delta_Teta(4,1),Delta_Teta(5,1),Delta_Teta(6,1)]
F=[max_f1;max_f2;max_f3;max_f4;max_f5;max_f6]
maxim=max(F)
o1=find(F==maxim)
disp('the most vulnerable PMU is :')
disp(o1)
pmu=[2,4,6,7,10,14]
pmu1=pmu(o1)

b1=zeros(Nb,1);
for k=1:1:14
    if k==pmu1
        b1(k)=1
    end
end
b1

%% new Gama

 Gam_new=zeros(num_meas,num_meas);
 
g1=cos(delta_teta_spf(o1));
g2=sin(delta_teta_spf(o1));
gama1=[g1,-g2;g2,g1] ;


 for i=1
    
         
      sm=sum(Mn(1:i-1));

     o=Mn(i);
     j=1;
                     while j<=o
              
                          Gam_new(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama1);
                       j=j+2
                      end
          
     
    
 end 

 
 
%% step 2
b=zeros(Nb,1)
NP2=1;
comb=combntns(Nb,NP2);
b_c=zeros(Nb,comb);
t=find(b1==1);
b_c(t,:)=1;
for i=1:1:comb
   if i~=t
       b(i)=1;
       b_c(i,i)=1
   end
end
c=zeros(Nb,comb);
% b_c=zeros(Nb,comb);
L=0;

 for i=1:1:Nb
      b_bar=b_c(:,i);
      k=1;
      j=1;
            while k<=Nb 
                   if b_bar(k)<=a(k)
                   k=k+1;
                   else
            break
                   end
            c(:,i)=b_bar;
                   if j==Nb
                     L=L+1  ;
                    b_c(:,i)=c(:,i)  ; 
                   end
             j=j+1;
         
           end
 end
 k=0
  for i=1:1:comb
     if a(i)==0
         b_c(:,i-k)=[]
         k=k+1
     end
  end
 c=b_c
  [row,clm]=size(c)
 k=0; 
 for i=1:1:clm
     sm=sum(c(:,i-k))
     if sm==1
         c(:,i)=[]
         k=k+1
     end
 end
 

 for k=1:1:clm-1
    
     c(:,k)
     o=find(c(:,k)==1)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
          if o==2 & o~=pmu1
  %initialization  
int1=zeros(Nb,1)
int2=-delta_teta_max.*b(:,k)
int3=delta_teta_max.*b(:,k)
int=[int1,int2,int3]
% 
%norm2= norm(B_ML*v,2)
% X0=[int1,int2,int3]
for l=2:1:3

A=[];
B=[]
Aeq=[]
beq=[]
lb1=-delta_teta_max.*b(:,k)
ub1=delta_teta_max.*b(:,k)
%%%%
X0=int(l)
lb=lb1(o)
ub=ub1(o)
Gam_new;
% f1=fun_4bus1(delta_teta_spf1)
delta_teta_spf1=fmincon(@fun_6bus1,X0,A,B,Aeq,beq,lb,ub)

Delta_Teta(k,l)=delta_teta_spf1
f1=fun_6bus1(delta_teta_spf1)
F1=-f1
max_f1=max(F1)
end
          end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
          if o==4 & o~=pmu1
              
               %initialization  
int1=zeros(Nb,1)
int2=-delta_teta_max.*b(:,k)
int3=delta_teta_max.*b(:,k)
int=[int1,int2,int3]
% 
%norm2= norm(B_ML*v,2)
% X0=[int1,int2,int3]
for l=1:1:3

A=[];
B=[]
Aeq=[]
beq=[]
lb1=-delta_teta_max.*b(:,k)
ub1=delta_teta_max.*b(:,k)
%%%%
X0=int(l)
lb=lb1(o)
ub=ub1(o)
Gam_new;
delta_teta_spf2=fmincon(@fun_6bus2,X0,A,B,Aeq,beq,lb,ub)
Delta_Teta(k,l)=delta_teta_spf2 
f2=fun_6bus2(delta_teta_spf2)
F2=-f2  
max_f2=max(F2)
end  
          end
          
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
          if o==6 & o~=pmu1
              
               %initialization  
int1=zeros(Nb,1)
int2=-delta_teta_max.*b(:,k)
int3=delta_teta_max.*b(:,k)
int=[int1,int2,int3]
% 
%norm2= norm(B_ML*v,2)
% X0=[int1,int2,int3]
for l=1:1:3

A=[];
B=[]
Aeq=[]
beq=[]
lb1=-delta_teta_max.*b(:,k)
ub1=delta_teta_max.*b(:,k)
%%%%
X0=int(l)
lb=lb1(o)
ub=ub1(o)
Gam_new;
delta_teta_spf3=fmincon(@fun_6bus3,X0,A,B,Aeq,beq,lb,ub)
Delta_Teta(k,l)=delta_teta_spf3 
f3=fun_6bus3(delta_teta_spf3)
F3=-f3  
max_f3=max(F3)
end  
          end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
          if o==7 & o~=pmu1
              
               %initialization  
int1=zeros(Nb,1)
int2=-delta_teta_max.*b(:,k)
int3=delta_teta_max.*b(:,k)
int=[int1,int2,int3]
% 
%norm2= norm(B_ML*v,2)
% X0=[int1,int2,int3]
for l=1:1:3

A=[];
B=[]
Aeq=[]
beq=[]
lb1=-delta_teta_max.*b(:,k)
ub1=delta_teta_max.*b(:,k)
%%%%
X0=int(l)
lb=lb1(o)
ub=ub1(o)
Gam_new;
delta_teta_spf4=fmincon(@fun_6bus4,X0,A,B,Aeq,beq,lb,ub)
Delta_Teta(k,l)=delta_teta_spf4 
f4=fun_6bus4(delta_teta_spf4)
F4=-f4  
max_f4=max(F4)
end  
          end  
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
          if o==10 & o~=pmu1
              
               %initialization  
int1=zeros(Nb,1)
int2=-delta_teta_max.*b(:,k)
int3=delta_teta_max.*b(:,k)
int=[int1,int2,int3]
% 
%norm2= norm(B_ML*v,2)
% X0=[int1,int2,int3]
for l=1:1:3

A=[];
B=[]
Aeq=[]
beq=[]
lb1=-delta_teta_max.*b(:,k)
ub1=delta_teta_max.*b(:,k)
%%%%
X0=int(l)
lb=lb1(o)
ub=ub1(o)
Gam_new;
delta_teta_spf5=fmincon(@fun_6bus5,X0,A,B,Aeq,beq,lb,ub)
Delta_Teta(k,l)=delta_teta_spf5 
f5=fun_6bus5(delta_teta_spf5)
F5=-f5  
max_f5=max(F5)
end  
          end  
          
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
          if o==14 & o~=pmu1
              
               %initialization  
int1=zeros(Nb,1)
int2=-delta_teta_max.*b(:,k)
int3=delta_teta_max.*b(:,k)
int=[int1,int2,int3]
% 
%norm2= norm(B_ML*v,2)
% X0=[int1,int2,int3]
for l=1:1:3

A=[];
B=[]
Aeq=[]
beq=[]
lb1=-delta_teta_max.*b(:,k)
ub1=delta_teta_max.*b(:,k)
%%%%
X0=int(l)
lb=lb1(o)
ub=ub1(o)
Gam_new;
delta_teta_spf6=fmincon(@fun_6bus6,X0,A,B,Aeq,beq,lb,ub)
Delta_Teta(k,l)=delta_teta_spf6 
f6=fun_6bus6(delta_teta_spf6)
F6=-f6  
max_f6=max(F6)
end  
          end  
    
 end

% compare
delta_teta_spf2=[delta_teta_spf1,delta_teta_spf2,delta_teta_spf3,delta_teta_spf4,delta_teta_spf5,delta_teta_spf6]
F=[max_f1;max_f2;max_f3;max_f4;max_f5;max_f6]
F(o1)=0
maxim=max(F)
o2=find(F==maxim)
disp('the second vulnerable PMU is :')
disp(o2)
pmu2=pmu(o2);
disp('bus number=')
disp(pmu2)
b2=b1;
for k=1:1:14
    if k==pmu2 
        b2(k)=1;
    end
end
b2

