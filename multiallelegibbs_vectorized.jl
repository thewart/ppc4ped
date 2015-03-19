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
  foundloc = ped[:,2:3] .== 0;
  nfound = sum(foundloc);
  missloc = geno .< 0;
  nmiss = sum(missloc);

  maf0 = Array(Float64,g);
  if isempty(maf0)
    for l in 1:m
      maf0[l] = mean(geno[~missloc[:,l],l])/2
    end
  end

  if isempty(z0)
       z0 = Array(Int64,(n,2,m))
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

# #calculate likelihood of assayed genotypes under all possible true genotypes
# #these only need to be calculated once, then filled in with the latest error prob.
# #do not depend on any other parameters.
#   zster = ones((4,n))
#   zlik = Array(Float64,(4,n))
#   for i in 1:n
#     if geno[i] >= 0
#       for j in 1:4
#         zster[j,i] = (sum(zstates[j,:]) != geno[i])
#       end
#     end
#   end

  #initialize gibbs chain
  z = Array(Int64,(n,2,m,iter+1));
  ϵ = Array(Float64,iter+1);
  maf = Array(Float64,m,iter+1);
  badz = falses(iter);

  z[:,:,:,1] = znow = z0;
  ϵ[1] = ϵ0;
  maf[:,1] = maf0;

  pzflib = [[0,1,0,1]  #sire has genotype 0
            [.5,.5,.5,.5] #sire has genotype 1
            [1,0,1,0]]; #sire has genotype 2
  pzmlib = [[0,0,1,1]  #dam has genotype 0
            [.5,.5,.5,.5] #dam has genotype 1
            [1,1,0,0]]; #dam has genotype 2
  pzklib = [[0,.5,.5,1] #child has mat/paternal major allele
            [1,.5,.5,0]]; #child has mat/paternal minor allele

  for t in 1:iter

    zlik = construct_zlik(ϵ[t]);
    for i in 1:n

    #calculate genotype probability given base MAF and/or parental genotypes
      sire = ped[i,2]
      dam = ped[i,3]

      if sire > 0 #father in pedigree; sampled from father's alleles
        sirezg = sum(znow[sire,:,:],2)[:];
        pzf = pzflib[:,sirezg];
      else #father not in pedigree; sample from base population
        pzf = [maf[:,t]'; 1-maf[:,t]'; maf[:,t]'; 1-maf[:,t]'];
      end

      if dam > 0 #mother in pedigree; sampled from father's alleles
        damzg = sum(znow[dam,:,:],2)[:];
        pzm = pzmlib[:,damzg];
      else #mother not in pedigree; sample from base population
        pzm = [maf[:,t]'; maf[:,t]'; 1-maf[:,t]'; 1-maf[:,t]'];
      end

      pk = ones(Float64,(4,m));
    #calculate offspring genotype probabilities (if there are any)
      if any(ped[:,2] .== i)
        kid = find( ped[:,2] .== i)
        for j in kid
         #pk = pk .* (znow[j,1] == 1 ? [1,0.5,0.5,0] : [0,0.5,0.5,1]);
          p
        end
      elseif any(ped[:,3] .== i)
        kid = find( ped[:,3] .== i)
        for j in kid
         pk = pk .* (znow[j,2] == 1 ? [1,0.5,0.5,0] : [0,0.5,0.5,1]);
        end
      end

    #calculate conditional posteriors
      post = pk .* pzm .* pzf .* zlik[:,geno[i]+2]

    #check for impossible z
      if sum(post) == 0
        post[:] = 1;
        badz[t] = true;
      end
    #sample new
      z[i,:,t+1] = zstates[rand(Categorical(post./sum(post))),:];
    #simplify bookeeping by keeping most recent zs represented separately from iteration history
      znow[i,:] = z[i,:,t+1]
    end

  #update minor allele frequency
    nma = sum(znow[find(foundloc)]);
    maf[t+1] = rand( Beta(nma + mafprior[1], nfound - nma + mafprior[2]) );

  #update genotype error probability
    nerr = sum( ( sum(znow,2) .!= geno )[find(~missloc)] );
    ϵ[t+1] = rand( Beta(nerr + ϵprior[1],n - nmiss - nerr + ϵprior[2]) );

  end
  genosim = reshape(sum(z[:,:,2:(iter+1)],2),(n,iter));
  return genosim,maf[2:(iter+1)],ϵ[2:(iter+1)],z[:,:,2:(iter+1)],badz

end

