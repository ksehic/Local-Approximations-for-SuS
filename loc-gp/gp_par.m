function [T11R,Tr1,Y11,Rpls11] = gp_par(i,Ndef1,id1,X01,X,Y1,epsilon)

[Rpls11,Tr1] = pls1(X(id1(1:Ndef1+i),:),Y1(id1(1:Ndef1+i)),epsilon);

Y11 = Y1(id1(1:Ndef1+i));

XmeanLoc = mean(X(id1(1:Ndef1+i),:));

T11 = X01-XmeanLoc;

T11R = T11*Rpls11;

return