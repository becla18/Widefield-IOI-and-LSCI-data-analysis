function t=get_tc(K,exp_time, varargin)

  if(nargin>2)
      beta=varargin{1};
  else
      beta=1;
  end
  tc_tmp=logspace(-10,-1,50001); 
  x = exp_time./tc_tmp;
  K_tmp=sqrt(beta*(exp(-2*x)-1+2*x)./(2*x.^2));
  
  if(ndims(K)==2 & min(size(K))==1)
      t=interp1(K_tmp,tc_tmp,K);
  else
      Kfoo=reshape(K,[1 prod(size(K))]);
      t=interp1(K_tmp,tc_tmp,Kfoo);
      t=reshape(t,size(K));
  end
  