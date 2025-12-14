function kcw=kercw(sig,rowaa,colaa,fc);
kcw=[];
if ceil(colaa/2)==colaa/2
	asseteta=-fc/2:fc/colaa:fc/2-fc/colaa;
else
	asseteta=-(fc/2-fc/(colaa*2)):fc/colaa:fc/2-fc/(2*colaa);
end
if ceil(rowaa/2)==rowaa/2
i=1;
for assetau=-rowaa/fc:2/fc:(rowaa-1)/fc;
kcw(i,:)=exp(-(assetau^2*asseteta.^2)/sig);
i=i+1;
end
else
i=1;
for assetau=-(rowaa-1)/fc:2/fc:(rowaa-1)/fc;
kcw=exp(-(assetau.^2*asseteta.^2)/sig);
i=i+1;
end
end
%end
