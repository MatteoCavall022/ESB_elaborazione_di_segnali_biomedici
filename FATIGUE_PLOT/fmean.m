function [fmean,fmedian] = fmean(x, fsamp, epoch_len)
x=x-mean(x); 
NFFT =fsamp*epoch_len;

%per la frequenza media
[P, f] = pwelch(x, rectwin(fsamp*epoch_len), 0, NFFT ,fsamp);
fmean=sum(f.*P)/sum(P);

%per la frequenza mediana
[Px,fm]=pwelch(x,rectwin(fsamp*epoch_len),0,NFFT,fsamp);
A=sum(Px)/2;
k=1;
S=0;

while S<=A
   S=S+Px(k);
   k=k+1;
end
   fmedian=fm(k);
end

% - se voglio in output entrambe [fmean,fmedian] = fmean(x, fsamp, epoch_len)
% - se voglio in output fmean --> fmean_output = fmean(x, fsamp, epoch_len)
% - voglio in output fmedian --> [~, fmedian]=fmean(x, fsamp, epoch_len)
% per fare tilde alt e 126 sul tastierino