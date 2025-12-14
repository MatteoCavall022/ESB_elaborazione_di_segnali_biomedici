close all
 clear;
 clc;

 fprintf('\n Removal of line artifacts from ECG - Different filters \n\n');
 fprintf(' Press a key to choose an ecg50xxx.bin file \n\n');
 pause

%----------------
% Load data files
%----------------
 [file,path]=uigetfile('*.bin','Load');
 filename=sprintf('%s%s',path,file);
 h=fopen(filename,'r');
 x=fread(h,inf,'float');
 fclose(h);
%------------------------------------------
% Load of the reference signal ecgor001.bin
%------------------------------------------
 h=fopen('ecgor001.bin','r');
 xor=fread(h,inf,'float');
 fclose(h);
 z=0.01;
 B=2;
 fc=50;
 fs=512;
 ft=50;
 [b,a]=filtro_ricorsivo_interf_rete(z,B,fc,fs,ft,x,xor);
