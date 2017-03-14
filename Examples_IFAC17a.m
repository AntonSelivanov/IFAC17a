% This MATLAB program checks the feasibility of the LMIs from Propositions 1, 2 and Remark 2 of the paper 
% A. Selivanov and E. Fridman, "Simple conditions for sampled-data stabilization by using artificial delay," in 20th IFAC World Congress, 2017. 

%% Example 1
B=[0; 1]; C=[1 0]; K1=-.35; K2=.1; 
q=3; 
alpha=1e-3; 

% Proposition 1
A=[0 1; .1 0]; 

h0=.1; 
hmin=fminsearch(@(h) LMI_IFAC17a_Prop1(A,B,C,K1,K2,h,q,alpha),h0); 
hmax=fminsearch(@(h) -LMI_IFAC17a_Prop1(A,B,C,K1,K2,h,q,alpha),h0); 
if hmax ~= hmin
    disp(['Example 1, Proposition 1: h in [' num2str(hmin) ',' num2str(hmax) ']']); 
else
    disp('Example 1, Proposition 1: Not Feasible'); 
end

% Remark 2
A={[0 1; -.1 0],[0 1; .1 0]}; % Uncertain matrix 

h0=.1; 
hmin=fminsearch(@(h) LMI_IFAC17a_Rem2(A,B,C,K1,K2,h,q,alpha),h0); 
hmax=fminsearch(@(h) -LMI_IFAC17a_Rem2(A,B,C,K1,K2,h,q,alpha),h0); 
if hmax ~= hmin
    disp(['Example 1, Remark 2: h in [' num2str(hmin) ',' num2str(hmax) ']']); 
else
    disp('Example 1, Remark 2: Not Feasible'); 
end

%% Example 2
A=[0 1 0; 1 1 1; 1 0 -1]; B=[0; 1; 0]; C=[1 0 0]; K1=-17; K2=13; 
alpha=0; 
q=5; 
sigma=9e-4; 
nu=.5; % irrelevant for l=1

h0=.06; 
hmax=fminsearch(@(h) -LMI_IFAC17a_Prop2(A,B,C,K1,K2,h,q,alpha,sigma,nu),h0); 
if hmax ~= hmin
    disp(['Example 2, Proposition 2: hmax=' num2str(hmax)]); 
    [~,Omega]=LMI_IFAC17a_Prop2(A,B,C,K1,K2,hmax,q,alpha,sigma,nu); 
else
    disp('Example 2, Proposition 2: Not Feasible'); 
end