% USAGE
%  C= h2delay(S1,S2,T,n,method);
% FUNCTION:
%  gives H2 assymetric association index for delayed signals
% INPUTS:
%  S1,S2 : 1D signals;
%  T : array of S1 delays, (T>0 means S1 is after S2); default is [-100:100]
%  n : number of bins for S2, default [100]
%  method: string with a distance function method [h2_m] 
% OUTPUT:
%   C - 1D association function
% USES:
%   h2_c.m
% VERSION:
%  Stiliyan Kalitzin 24.08.01

function C = h2delay(S1, S2, T, n, method)

if nargin < 2;          help h2delay;       return; end

if nargin < 3;          T = [];                     end
if isempty(T);          T = -100:100;               end

if nargin < 4;          n = [];                     end
if isempty(n);          n = 100;                    end

if nargin < 5;          method = [];                end
if isempty(method);     method = 'h2_m';            end

% S1=detrend(S1); 
% S2=detrend(S2); 
C=zeros(1,numel(T)); 
for k=1:length(T)
   if T(k)>=0 
       eval(['c=',method,'((S1(:,T(k)+1:end)),(S2),n);']);
   else
       eval(['c=',method,'((S1),S2(:,-T(k)+1:end),n);']); 
   end
if isempty(c)
    continue;
end
if isnan(c)
    continue;
end
C(k)=c; 
end

% function P = pci(A,B,varargin)
% m=min(length(A),length(B));
% A=A(1:m); B=B(1:m);
% a=sum(abs(A.*B)); 
% if a>0 
%     P=sum(A.*conj(B))./a; 
% else P=0; 
% end 