function [] = allatonce()
allN=[40,80];
allK=[4];

list=zeros(max(size(allN))*max(size(allK)),3);
count=0;
for nn=1:max(size(allN))
   for kk=1:max(size(allK))
       count=count+1;
       list(count,:)=...
           [allN(nn),allK(kk),galerkin_schr_1dimension(allN(nn),allK(kk))]
   end
end

list

end

