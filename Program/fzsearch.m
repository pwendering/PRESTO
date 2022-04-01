function [d,varargout]  =  fzsearch(r,p,n,cas)

% fzsearch (R, P, N, CAS) finds the best or predetermined approximate
% matching between the substrings of a string R (reference) and a string P
% (pattern). The Levenshtein distance is used as a measure of matching.
% Levenshtein distance is the minimal number of substitutions, deletions
% and insertions of characters required to convert string A to string B.
% D  =  fzsearch(R) computes NaN.
% D  =  fzsearch(R,P) computes the best matching between substrings of the 
% and a pattern.
% D is a vector, where D(1) is a distance and D(2), D(3)... are indexes
% of the ends of the best matching substrings.
% [D,A]  =  fzsearch(R,P) computes as above and a cell array A of the best 
% matching substrings.
% D  =  fzsearch(R,P,N) computes the match between substrings of R and P in 
% interval from the best match (say, M) to M + N. D is a (N + 1) cell array.
% Each cell contains a vector which first number is a distance B + K, K  = 0..N 
% and others are indexes of the ends of substrings with the same distance from P.
% D  =  fzsearch (R, P, N, CAS) is the same of the previous but the case is
% ignored for CAS > 0
% D  =  fzsearch('','') computes D  =  0.
% D  =  fzsearch('',P) computes D(1)  =  numel(P) and D(2)  =  0.
% D  =  fzsearch(R,'') computes D(1)  =  1 and D(i)  =  i-1, i  =  2..numel(R)+1.
% Example
% reference = '1713512737451262';
% pattern = '2345';
% [d,A]=fzsearch(reference,pattern);
% fprintf('A distance of the best matching: %2.0f\n',d(1))
% disp('Indexes of the ends of substrings:')
% disp(d(2:end))
% disp('A set of substrings:')
% disp(A)
% A distance of the best matching:  2
% Indexes of the ends of substrings:
%      5    12
% A set of substrings:
%     '135'
%     '35'
%     '273745'
%     '73745'
%     '3745'
%     '745'
%     '45'
varargout(1) = {''};
if nargout == 2 &&  nargin ~= 2
  warning...
    ('Calculation of the set of substrings is possible only for 2 input arguments')  
end
switch nargin
  case 1
    d = NaN;
    return
  case 2    
    n = 0;  
  case 4    
    if cas > 0
      p = upper(p);
      r = upper(r);
    end
end
if isempty(p)&&isempty(r)
  d = 0;
  return 
end
if isempty(p)
  d = [1,(1:numel(r))];
  return
end
if isempty(r)
  d = [numel(p),0];
  return
end
dl=stdis(r,p,0); 
m = min(dl(end,2:end)); % calculation of distance of the best match
if nargin > 2 
  d(n+1) = {0}; % memory allocation
  for k = 0:n
% Calculation of indexes for distance m + k     
    d(k+1) = {[m+k,find(dl(end,2:end) == m+k)]};
  end
else
  d = [m,find(dl(end,2:end) == m)];
end
if nargout == 2 && nargin == 2
% Calculation of the set of substrings
  z(100,1) = {''};
  i = 0;
  for k = 2:numel(d)
    s = max(1,d(k)-numel(p)-d(1)+1); 
    for j = 1:2*d(1)+1 
      c = r(s+j-1:d(k));
      dl = stdis(c,p,1);
      if  dl(end,end )== d(1)
        i = i+1;
        z(i) = {c};
      end
    end
  end
  varargout(1) = {z(1:i)};
end
  
