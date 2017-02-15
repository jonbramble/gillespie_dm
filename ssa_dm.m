% Stochastic Simulation Algorithms

% x0 initial value of state vector
% ak propensity factors - ie probability of reactions
% akx propensity matrix - relates propensity factors to their components in x
% rmax - maximum number of interactions
% nu - stoichiometry matrix

function [tt,xx] = ssa_dm(x0,ak,akx,rmax,nu)
  
  rnums = rand(2,rmax); % 2 random numbers for each step
  tt = zeros(rmax+1,1);
  xx = zeros(rmax+1,3);

  x=x0; %set the initial state of the system
  rk = 0; % reaction counter
  t = 0; % current time
  
  ajp = zeros(1,4);
  
  while rk < rmax
      rk = rk + 1;
      r = akx*x';
      aa = ak.*r';  %this is the propensity array
      
      a0 = sum(aa);
      
      r1 = rnums(1,rk);
      r2 = rnums(2,rk);
    
      dt = log(1/r1)/a0; % get the nxt time point
    
      csaa = 0;
      for j=1:4
       csaa = csaa + aa(j);  
       ts = (csaa >  r2*a0);       
       ajp(j) = ts;
      end
    
      ss = find(ajp,1);
      nu_state = nu(ss,:);        %the chosen reaction - chose the row
   
      t = t + dt;                 %update the time
      tt(rk+1) = t;               %track the time
    
      x = x + nu_state;           %update the dynamic state
      xx(rk+1,:) = x;             %track the state
  end
 
end

