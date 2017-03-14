function [hfeas,OmegaVal]=LMI_IFAC17a_Prop2(A,B,C,K1,K2,h,q,alpha,sigma,nu)
% This MATLAB program checks the feasibility of the LMIs from Proposition 2 of the paper 
% A. Selivanov and E. Fridman, "Simple conditions for sampled-data stabilization by using artificial delay," in 20th IFAC World Congress, 2017. 

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)

% Input: 
% A, B, C       - parameters of the system (1) 
% K1, K2        - controller gains form (11) 
% h             - the sampling period
% q             - artificial delay from (11) 
% alpha         - desired decay rate 
% sigma         - event-triggering parameter from (13)
% nu            >0 - tuning parameter 

% Output: 
% hfeas         - equals h if feasible, 0 otherwise. Such output is convenient for fminsearch(); 
% OmegaVal      - value of the decision matrix Omega

%% Notations
D=A+B*(K1+K2)*C-q*h*B*K2*C*A; 
[l,n]=size(C); 

%% Decision variables 
P=sdpvar(n); 
P2=sdpvar(n,n,'f'); 
P3=sdpvar(n,n,'f'); 
Q=sdpvar(l); 
W=sdpvar(l); 
U=sdpvar(l); 
Omega=sdpvar(l); 
if l == 1
    Omega2=sdpvar; 
end

%% The LMIs for Psi
% tau=0 
Psi0=blkvar; 
Psi0(1,1)=2*alpha*P+P2'*D+D'*P2; 
Psi0(1,2)=P-P2'+D'*P3; 
Psi0(1,3)=P2'*B*K2; 
Psi0(1,4)=P2'*B*K2; 
Psi0(1,5)=P2'*B*K1; 
Psi0(1,6)=P2'*B*K2; 
Psi0(1,7)=sigma*C'*Omega; 
Psi0(2,2)=-P3-P3'+h*C'*U*C+(q*h)^2*A'*C'*Q*C*A+h^2*exp(2*alpha*h)*C'*W*C; 
Psi0(2,3)=P3'*B*K2; 
Psi0(2,4)=P3'*B*K2; 
Psi0(2,5)=P3'*B*K1; 
Psi0(2,6)=P3'*B*K2; 
Psi0(3,3)=-4/(q*h)^2*exp(-2*alpha*q*h)*Q; 
Psi0(4,4)=-pi^2/4*exp(-2*alpha*q*h)*W; 
Psi0(5,5)=-Omega; 
Psi0(7,7)=-sigma*Omega; 
if l == 1
    Psi0(1,8)=sigma*C'*Omega2; 
    Psi0(2,8)=-sigma*q*h*C'*Omega2; 
    Psi0(3,8)=sigma*Omega2; 
    Psi0(4,8)=sigma*Omega2; 
    Psi0(6,6)=-Omega2; 
    Psi0(8,8)=-sigma*Omega2; 
else
    Psi0(1,8)=sigma*C'*Omega; 
    Psi0(2,8)=-sigma*q*h*C'*Omega; 
    Psi0(3,8)=sigma*Omega; 
    Psi0(4,8)=sigma*Omega; 
    Psi0(6,6)=-nu*Omega; 
    Psi0(8,8)=-sigma/nu*Omega; 
end
Psi0=sdpvar(Psi0); 

% tau=h 
Psih=blkvar; 
Psih(1,1)=2*alpha*P+P2'*D+D'*P2; 
Psih(1,2)=P-P2'+D'*P3; 
Psih(1,3)=P2'*B*K2; 
Psih(1,4)=P2'*B*K2; 
Psih(1,5)=h*P2'*B*K1; 
Psih(1,6)=P2'*B*K1; 
Psih(1,7)=P2'*B*K2; 
Psih(1,8)=sigma*C'*Omega; 
Psih(2,2)=-P3-P3'+(q*h)^2*A'*C'*Q*C*A+h^2*exp(2*alpha*h)*C'*W*C; 
Psih(2,3)=P3'*B*K2; 
Psih(2,4)=P3'*B*K2; 
Psih(2,5)=h*P3'*B*K1; 
Psih(2,6)=P3'*B*K1; 
Psih(2,7)=P3'*B*K2; 
Psih(3,3)=-4/(q*h)^2*exp(-2*alpha*q*h)*Q; 
Psih(4,4)=-pi^2/4*exp(-2*alpha*q*h)*W; 
Psih(5,5)=-h*exp(-2*alpha*h)*U; 
Psih(5,8)=sigma*h*Omega; 
Psih(6,6)=-Omega; 
Psih(8,8)=-sigma*Omega; 
if l==1
    Psih(1,9)=sigma*C'*Omega2; 
    Psih(2,9)=-sigma*q*h*C'*Omega2; 
    Psih(3,9)=sigma*Omega2; 
    Psih(4,9)=sigma*Omega2; 
    Psih(7,7)=-Omega2; 
    Psih(9,9)=-sigma*Omega2; 
else
    Psih(1,9)=sigma*C'*Omega; 
    Psih(2,9)=-sigma*q*h*C'*Omega; 
    Psih(3,9)=sigma*Omega; 
    Psih(4,9)=sigma*Omega; 
    Psih(7,7)=-nu*Omega; 
    Psih(9,9)=-sigma/nu*Omega; 
end
Psih=sdpvar(Psih); 
    
LMIs=[P>=0, Q>=0, W>=0, U>=0, Omega>=0, Psi0<=0, Psih<=0];
if l == 1
    LMIs=[LMIs, Omega2>=0]; 
end

%% Solution of LMIs
options=sdpsettings('solver','lmilab','verbose',0); 
sol=optimize(LMIs,[],options); 

hfeas=0; OmegaVal=[]; 
if sol.problem == 0 
    primal=check(LMIs); 
    if min(primal)>0 
        hfeas=h; 
        OmegaVal=value(Omega); 
    end
else
    yalmiperror(sol.problem) 
end