function y=cct(x1,x2,ntlag)
%act Instantaneous Crosscorrelation Function
%	y=cct(x1,x2,nfft) estimates the instantaneous autocorrelation function
%	of an input signal x. The number of time lags considered to estimate
%	the instantaneous autocorrelation function is equal to ntlag. It
%	corresponds to the number of rows of the output matrix unless
%	ntlag is greater than the number of samples of the input signal.
%	Each row of y contains the products x[n+k]*conj(x[n-k]) for a 
%	value of k, where n is related to the time axis and K is related 
%	to the time lag axis. ntlag should be smaller than the number of
%	samples of the signal vector x. If the number of points of the
%	input vector is smaller than ntlag, the number of time lags
%	is limited to the number of samples of the input signal.

%	input variable(s):	x	input signal
%				ntlag	number of time lags

%	output variable(s):	y	instantaneous autocorrelation function

%	author(s):	G. Gagliati
%			P. Bonato, 4-22-1996, revised


n=max(size(x1));
row=min(ntlag,n);
y=zeros(row,n);		% initialize the output matrix y

y(1,:)=x1.*conj(x2);	% zero lag instantaneous autocorrelation function 
for i=1:row/2-1
	y(i+1,:)=[x1(i+1:n) zeros(1,i)].*conj([zeros(1,i) x2(1:n-i)]);
			% positive lag instantaneous autocorrelation function
	y(row-i+1,:)=conj(y(i+1,:));
			% negative lag instantaneous autocorrelation function
end
i=fix(row/2);
y(i+1,:)=[x1(i+1:n) zeros(1,i)].*conj([zeros(1,i) x2(1:n-i)]);
			% midtime lag instantaneous autocorrelation function
if i*2+1==row
	y(i+2,:)=conj(y(i+1,:));
			% if row is odd we also have a conj element of y(i+1,:)
end

fprintf(' Instantaneous Autocorrelation function estimate: done!\n')


