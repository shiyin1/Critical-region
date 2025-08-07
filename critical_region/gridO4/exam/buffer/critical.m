clear
kaft=load('k_aft.dat');
mpi2min=load('mpi2min.dat');
msg2min=load('msg2min.dat');
Vertex2pisg=load('Vertex2pisg.dat');
Vertex3sg=load('Vertex3sg.dat');
piML=load('pi_ML_k.dat');
piQL=load('pi_QL_k.dat');
sgML1=load('sg_ML1_k.dat');
sgML2=load('sg_ML2_k.dat');
sgQL=load('sg_QL_k.dat');
%n1=length(k);
%etaphik=reshape(etaphi,[81,n1]);
%Zphik=reshape(Z_phi,[81,n1]);
%V0=reshape(V0k,[9001,n1]);
%plot(phi,etaphiferk(:,(60:69)))
figure(1)
plot(kaft,piML)
title('piML')
%figure(2)
%plot(kaft,piQL)
%title('piQL')
figure(3)
plot(kaft,sgML1)
title('sgML1')
figure(4)
plot(kaft,sgML2)
title('sgML2')
%figure(5)
%plot(kaft,sgQL)
%title('sgQL')
%plot(kaft,msg2min,'r')
%plot(kaft,sgML1,'b')
%plot(k,phi0)
%plot(k,mpi)
%hold on
%plot(k,msg)
hold off
%legend('T=146,mu=0','')
%for i=1:1
%   figure(1)
%plot(phi./sqrt(Z_phi1(:,i)),dVdrho01(:,i))
%end
%figure(2)
%plot(phi,10*(phi.^2/2)+3.33677*10E-308*(phi.^2/2).^80)
%plot(phi,Z_phi2(:,1:34))
%for i=1:34
%plot(phi,phi./sqrt(Z_phi2(:,i)))
%hold on
%end
