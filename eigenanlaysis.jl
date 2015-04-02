function pedigree_eigenanalysis(X;q=[0.025,0.975],pratio=1.0)
  if ndims(X) == 3
    n,m,iter = size(X);
  else
    n,m = size(X);
    iter = 1;
  end

  eigen = zeros(Float64,(m,iter));
  Xnorm = Array(Float64,(n,m))
  for i in 1:iter
    for j in 1:m
      sd = std(X[:,j,i]);
      sd = (sd == 0) ? 1 : sd;
      Xnorm[:,j] = ( X[:,j,i] - mean(X[:,j,i]) ) / sd;
    end
    foo = fit(PCA,Xnorm',pratio=pratio);
    #foo = fit(PCA,float(X[:,:,i]'));
    eigen[1:outdim(foo),i] = principalvars(foo);
  end

  if iter > 1
    eigmean = mean(eigen,2);
    eigstd = std(eigen,2);
    eigci = Array(Float64,(m,2));
    for i in 1:m
      eigci[i,:] = quantile(vec(eigen[i,:]),q);
    end
    return eigen,eigmean,eigstd,eigci
  else
    return eigen
  end

end
