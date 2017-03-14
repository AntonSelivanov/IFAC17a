function hfeas=LMI_IFAC17a_Prop1(A,B,C,K1,K2,h,q,alpha)
% This MATLAB program checks the feasibility of the LMIs from Proposition 1 of the paper 
% A. Selivanov and E. Fridman, "Simple conditions for sampled-data stabilization by using artificial delay," in 20th IFAC World Congress, 2017. 

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)

% Input: 
% A, B, C       - parameters of the system (1) 
% K1, K2        - controller gains form (3) 
% h             - the sampling period
% q             - artificial delay from (3) 
% alpha         - desired decay rate 

% Output: 
% hfeas         - equals h if feasible, 0 otherwise. Such output is convenient for fminsearch(); 

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

%% The LMIs for Phi
% tau=0 
Phi0=blkvar; 
Phi0(1,1)=2*alpha*P+P2'*D+D'*P2; 
Phi0(1,2)=P-P2'+D'*P3; 
Phi0(1,3)=P2'*B*K2; 
Phi0(1,4)=P2'*B*K2; 
Phi0(2,2)=-P3-P3'+h*C'*U*C+(q*h)^2*A'*C'*Q*C*A+h^2*exp(2*alpha*h)*C'*W*C; 
Phi0(2,3)=P3'*B*K2; 
Phi0(2,4)=P3'*B*K2; 
Phi0(3,3)=-4/(q*h)^2*exp(-2*alpha*q*h)*Q; 
Phi0(4,4)=-pi^2/4*exp(-2*alpha*q*h)*W;     
Phi0=sdpvar(Phi0); 

% tau=h 
Phih=blkvar; 
Phih(1,1)=2*alpha*P+P2'*D+D'*P2; 
Phih(1,2)=P-P2'+D'*P3; 
Phih(1,3)=P2'*B*K2; 
Phih(1,4)=P2'*B*K2; 
Phih(1,5)=h*P2'*B*K1; 
Phih(2,2)=-P3-P3'+(q*h)^2*A'*C'*Q*C*A+h^2*exp(2*alpha*h)*C'*W*C; 
Phih(2,3)=P3'*B*K2; 
Phih(2,4)=P3'*B*K2; 
Phih(2,5)=h*P3'*B*K1; 
Phih(3,3)=-4/(q*h)^2*exp(-2*alpha*q*h)*Q; 
Phih(4,4)=-pi^2/4*exp(-2*alpha*q*h)*W; 
Phih(5,5)=-h*exp(-2*alpha*h)*U; 
Phih=sdpvar(Phih); 

LMIs=[P>=0,Q>=0,W>=0,U>=0,Phi0<=0,Phih<=0];

%% Solution of LMIs
options=sdpsettings('solver','lmilab','verbose',0); 
sol=optimize(LMIs,[],options); 

hfeas=0; 
if sol.problem == 0
    primal=check(LMIs); 
    if min(primal)>0
        hfeas=h; 
    end
else
    yalmiperror(sol.problem) 
end
