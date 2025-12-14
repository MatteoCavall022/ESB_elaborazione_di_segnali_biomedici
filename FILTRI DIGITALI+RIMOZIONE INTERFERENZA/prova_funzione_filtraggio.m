 close all
 clear;
 clc;
 
 fprintf('\n Rimozione del tremore muscolare dal tracciato ECG mediante differenti filtri \n\n');
 fprintf('Premere un tasto e caricare uno dei file ecgTMxxx.bin \n\n');
 pause

%-----------------------
% Load of the data files
%-----------------------
 [file,path]=uigetfile('*.bin','Load');
 filename=sprintf('%s%s',path,file);
 h=fopen(filename,'r');
 xtm=fread(h,inf,'float');
 fclose(h);
 
%------------------------------------------
% Load of the reference signal ecgor001.bin
%------------------------------------------
 h=fopen('ecgor001.bin','r');
 xor=fread(h,inf,'float');
 fclose(h);
fs=512;
typefilter='filtraggio_blocchi_IIR';
 [a,b,ytm] = filtri_digitali(xor,xtm,fs, typefilter);
 