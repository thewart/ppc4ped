function pedigree_genogibbs(geno,ped,iter=100,mafprior=(1,1),ϵprior=(1,20),
                            z0=[],maf0=[],ϵ0=0.1)

  function construct_zlik(ϵ)
    zlik = [1 1 1 1;          # geno[i]=-1 (not genotyped)
            ϵ/2 ϵ/2 ϵ/2 1-ϵ;  # geno[i]= 0
            ϵ/2 1-ϵ 1-ϵ ϵ/2;  # geno[i]= 1
            1-ϵ ϵ/2 ϵ/2 ϵ/2]; # geno[i]= 2
    return(zlik')
  end

  n,m = size(geno);
  foundloc = find(ped[:,2:3]' .== 0);
  nfound = length(foundloc);
  ifound,jfound = ind2sub((2,n),foundloc);
  missloc = geno .< 0;
  nmiss = sum(missloc);

  maf0 = Array(Float64,g);
  if isempty(maf0)
    for l in 1:m
      maf0[l] = mean(geno[~missloc[:,l],l])/2
    end
  end

  if isempty(z0)
       z0 = Array(Int64,(2,n,m))
#        z0[find(missloc),:] = rand(Bernoulli(maf0),sum(missloc)*2);
#        z0[find(geno.==0),:] = 0;
#        z0[find(geno.==2),:] = 1;
#        z0[find(geno.==1),:] = [1 0; 0 1][rand(Bernoulli(0.5),sum(geno.==1))+1,:];
    foo, z0[:,:,l] = pedigree_genosim(ped,maf0[l]);
  end

 zstates = Int64[1 1;
                 0 1;
                 1 0;
                 0 0];

  sireof,damof = pedigree_childfinder(ped);
  sireany = sum(sireof,1)[:] .> 0;
  damany = sum(damof,1)[:] .> 0;

  #initialize gibbs chain
  z = Array(Int64,(2,n,m,iter+1));
  ϵ = Array(Float64,iter+1);
  maf = Array(Float64,m,iter+1);
  badz = falses(iter);

  pk = ones(Float64,4);

  z[:,:,:,1] = znow = z0;
  ϵ[1] = ϵ0;
  maf[:,1] = maf0;

  pzflib = [[0,1,0,1]  #sire has genotype 0
            [.5,.5,.5,.5] #sire has genotype 1
            [1,0,1,0]]; #sire has genotype 2
  pzmlib = [[0,0,1,1]  #dam has genotype 0
            [.5,.5,.5,.5] #dam has genotype 1
            [1,1,0,0]]; #dam has genotype 2
  pklib = [[0,.5,.5,1] #child has mat/paternal major allele
            [1,.5,.5,0]]; #child has mat/paternal minor allele

  nma = zeros(Int64,m);

  for t in 1:iter
    zlik = construct_zlik(ϵ[t]);

    for i in 1:n
      sire = ped[i,2]
      dam = ped[i,3]
      for j in 1:m

    #calculate genotype probability given base MAF and/or parental genotypes

      if sire > 0 #father in pedigree; sampled from father's alleles
        pzf = pzflib[:,znow[1,sire,j] + znow[2,sire,j]];
      else #father not in pedigree; sample from base population
        pzf = [maf[j,t], 1-maf[j,t], maf[j,t], 1-maf[j,t]];
      end

      if dam > 0 #mother in pedigree; sampled from father's alleles
        pzm = pzmlib[:,znow[1,dam,j] + znow[2,dam,j]];
      else #mother not in pedigree; sample from base population
        pzm = [maf[j,t], maf[j,t], 1-maf[j,t], 1-maf[j,t]];
      end

      pk[:] = 1;
    #calculate offspring genotype probabilities (if there are any)
      if sireany[i]
        for k in find(sireof[:,j])
         pk = pk .* pklib[:,znow[1,k]+1];
        end
      elseif damany[i]
        for k in find(damof[:,j])
         pk = pk .* pklib[:,znow[2,k]+1];
        end
      end

    #calculate conditional posteriors
      post = pk .* pzm .* pzf .* zlik[:,geno[i,j]+2];

#     #check for impossible z
#       if sum(post) == 0
#         post[:] = 1;
#         badz[t] = true;
#       end

    #sample new allele
      z[:,i,j,t+1] = zstates[rand(Categorical(post./sum(post))),:];
    #simplify bookeeping by keeping most recent zs represented separately from iteration history
      znow[:,i,j] = z[:,i,j,t+1]

        #update founder minor allele counts
        if ped[i,2]==0
          nma[j] += znow[1,i,j];
        end
        if ped[i,3]==0
          nma[j] += znow[2,i,j];
        end

        #update genotyping error counts
        if (znow[1,i,j] + znow[2,i,j]) != geno[i,j]
          nerr += 1;
        end
    end
    end

  #update minor allele frequency and allele error frequency
    for j in 1:m
      maf[j,t+1] = rand( Beta(nma[j] + mafprior[1], nfound - nma[j] + mafprior[2]) );
      nma[j] = 0; #while we're at it reset maf
    end

  #update genotype error probability
    ϵ[t+1] = rand( Beta(nerr + ϵprior[1],n*m - nmiss - nerr + ϵprior[2]) );

  end
  genosim = reshape(sum(z[:,:,2:(iter+1)],1),(n,m,iter));
  return genosim,maf[:,2:(iter+1)],ϵ[2:(iter+1)],z[:,:,2:(iter+1)],badz

end

function pedigree_childfinder(ped)
  n = size(ped)[1];
  sireof = Array(Bool,(n,n));
  damof = Array(Bool,(n,n));
  for i in 1:n, j in 1:n
    sireof[i,j] = i == ped[j,2];
    damof[i,j] = i == dam[j,3];
  end
end
