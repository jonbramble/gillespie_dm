
ca = 0; cd = 10000; cm = 0;
k1 = 0.01; k2 = 1.0; k3 = 1; k4 = 0.01;

rmx = 120000; % total number of steps
x0 = [ca,cd,cm];
ak = [k1,k2,k3,k4];
nu = [0,-1,1;0,1,-1;1,0,-1;-1,0,1];
akx = [0,1,0;0,0,1;0,0,1;1,0,0];

[ta,xa] = ssa_dm(x0,ak,akx,rmx,nu);

plot(ta,xa)