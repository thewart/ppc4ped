function pedigree_eigenanalysis(X;q=[0.025,0.975],pratio=0.9995)
  if ndims(X) == 3
    n,m,iter = size(X);
  else
    n,m = size(X);
    iter = 1;
  end

  eigen = zeros(Float64,(m,iter));
  for i in 1:iter
    Xnorm = ( X[:,:,i] .- mean(X[:,:,i],1) ) ./ std(X[:,:,i],1);
    foo = fit(PCA,float(Xnorm'),pratio=pratio);
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
