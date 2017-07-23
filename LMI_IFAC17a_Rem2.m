function hfeas=LMI_IFAC17a_rem2(A,B,C,K1,K2,h,q,alpha)
% This MATLAB program checks the feasibility of the LMIs from Remark 2 of the paper 
% A. Selivanov and E. Fridman, "Simple conditions for sampled-data stabilization by using artificial delay," in 20th IFAC World Congress, 2017, pp. 13837–13841. 

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)

% Input: 
% A             - is a cell array of vetrices A^j from Remark 2
% B, C          - parameters of the system (1) 
% K1, K2        - controller gains form (3) 
% h             - the sampling period
% q             - artificial delay from (3) 
% alpha         - desired decay rate 

% Output: 
% hfeas         - equals h if feasible, 0 otherwise. Such output is convenient for fminsearch(); 

if ~iscell(A)
    A={A}; 
end
[l,n]=size(C); 

%% Decision variables 
P2=sdpvar(n,n,'f'); 
P3=sdpvar(n,n,'f'); 
Q=sdpvar(l); 
W=sdpvar(l); 
U=sdpvar(l); 
LMIs=[Q>=0,W>=0,U>=0];

for j=1:length(A) % loop over polytope vertices
    P=sdpvar(n); 
    D=A{j}+B*(K1+K2)*C-q*h*B*K2*C*A{j};     
    
    %% The LMIs for Xi
    % tau=0 
    Xi0=blkvar; 
    Xi0(1,1)=2*alpha*P+P2'*D+D'*P2; 
    Xi0(1,2)=P-P2'+D'*P3; 
    Xi0(1,3)=P2'*B*K2; 
    Xi0(1,4)=P2'*B*K2; 
    Xi0(2,2)=-P3-P3'+h*C'*U*C+h^2*exp(2*alpha*h)*C'*W*C; 
    Xi0(2,3)=P3'*B*K2; 
    Xi0(2,4)=P3'*B*K2; 
    Xi0(2,5)=q*h*A{j}'*C'*Q; 
    Xi0(3,3)=-4/(q*h)^2*exp(-2*alpha*q*h)*Q; 
    Xi0(4,4)=-pi^2/4*exp(-2*alpha*q*h)*W;     
    Xi0(5,5)=-Q; 
    Xi0=sdpvar(Xi0); 

    % tau=h 
    Xih=blkvar; 
    Xih(1,1)=2*alpha*P+P2'*D+D'*P2; 
    Xih(1,2)=P-P2'+D'*P3; 
    Xih(1,3)=P2'*B*K2; 
    Xih(1,4)=P2'*B*K2; 
    Xih(1,5)=h*P2'*B*K1; 
    Xih(2,2)=-P3-P3'+h^2*exp(2*alpha*h)*C'*W*C; 
    Xih(2,3)=P3'*B*K2; 
    Xih(2,4)=P3'*B*K2; 
    Xih(2,5)=h*P3'*B*K1; 
    Xih(2,6)=q*h*A{j}'*C'*Q; 
    Xih(3,3)=-4/(q*h)^2*exp(-2*alpha*q*h)*Q; 
    Xih(4,4)=-pi^2/4*exp(-2*alpha*q*h)*W; 
    Xih(5,5)=-h*exp(-2*alpha*h)*U; 
    Xih(6,6)=-Q; 
    Xih=sdpvar(Xih); 
    
    LMIs=[LMIs, P>=0, Xi0<=0, Xih<=0];  %#ok<*AGROW>
end

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
