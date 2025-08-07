clear
load('V.dat')
load('phi.dat')
V_save=reshape(V,[90,260]);
phi_i=(0:0.0001:100)';
for num=1:260
    V0(:,num)=interp1(phi,V_save(:,num),phi_i,'spline');
    phi0(num)=find(V0(:,num)==min(V0(:,num)))*0.0001;
end
