clear
load('a.dat')
delta_phi=1e-9;
phi=0:delta_phi:1;
c=1e-34;
for i=1:10:200
    for j=1:maxim
        summation=0.5*a(1);
        for k=2:gridnum
            summation=summation+a(k)*cos((k-1)*acos((2*phi(j)-rhomin-rhomax)/(rhomax-rhomin)));
        end
        V0(j)=summation;
    end
    V0=V0-c*i*phi;
    logc(i)=log(c*i);
    phi0(i)=find(V0==min(V0))*delta_phi;
    logphi0(i)=log(phi0(i));
end
plot(logc,logphi0);
