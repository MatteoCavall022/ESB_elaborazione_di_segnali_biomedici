function aa=daf(x,ntlag,fc,T);
% DAF - Discrete Ambiguity Function
% Usage:   aa=daf(x,ntlag,fc,T);
% 
% Inputs:  x      input signal
%          ntlag  number of time lag for the autocorrelation sequence estimate
%          fc     sampling rate in Hz
%          T      duration of the signal in seconds
% Output:  aa     matrix containing the DAF

%	author(s):	G. Gagliati
%			P. Bonato, 4-22-1996, revised

y=act(x,ntlag);
y=y';
aa=[];
aa=fft(y);
aa=aa';
aa=(1/fc)*(1/T)*aa;
fprintf('Discrete Ambiguity function estimate: done!\n');
