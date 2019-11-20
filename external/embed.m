% USAGE
%  [V,tau] = embed(S,tau,ne,skip);
% FUNCTION:
%    creates a time-delay embeded vector-signal from a 1D signal 
% INPUT:
%  S - 1D signal;
%  tau - <-60> embedding delay, if tau<0 minimal correlation method is used in the range [1,-tau]
%  ne - number od embedding dimensions
%  skip - have to ask, no need to know
% OUTPUT:
%  V - [embedding,time] vector signal
%  tau - the embedding delay (if automatically computed ), if <0 no zero-cross is detected
% VERSION:
%  Stiliyan, 12.08.2005


function [V,tau1] = embed(S,tau,ne,skip);
if nargin<3 skip=1; end; 
nt=size(S,2);
nc=size(S,1);
tau1=tau;
if tau<0
    tau=abs(tau);
    XC=xcov(S); XC=XC(nt+1:end);
    L=find(sign(XC(1:end-1))-sign(XC(2:end))~=0);
    if (length(L)>0) && (L(1)<tau) tau=L(1); tau1=tau; else [m,tau]=min(abs(XC(1:tau)));tau1=-tau; end;
    clear L XC;
end;
nv=floor((nt-1-(ne-1)*tau)/skip)+1;
%V(1:ne,1:nv)=0;
for k=1:ne 
    S1=S(1+(k-1)*tau:skip:end);
    S1=S1(1:nv);
    V(k,:)=reshape(S1,[1,nv]); 
end;
