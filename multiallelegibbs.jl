function pedigree_genogibbs(geno,ped,iter=100,thin=1,mafprior=(1,1),ϵprior=(1,20),
                            z0=[],maf0=[],ϵ0=0.1)

  function construct_zlik(ϵ)
    zlik = [1 1 1 1;          # geno[i]=-1 (not genotyped)
            ϵ/2 ϵ/2 ϵ/2 1-ϵ;  # geno[i]= 0
            ϵ/2 1-ϵ 1-ϵ ϵ/2;  # geno[i]= 1
            1-ϵ ϵ/2 ϵ/2 ϵ/2]; # geno[i]= 2
    return(zlik')
  end

  if ndims(geno) > 1
    n,m = size(geno);
  else
    n = size(geno)[1];
    m = 1;
  end
  foundloc = find(ped[:,2:3]' .== 0);
  nfound = length(foundloc);
  ifound,jfound = ind2sub((2,n),foundloc);
  missloc = geno .< 0;
  nmiss = sum(missloc);
  sire = ped[:,2];
  dam = ped[:,3];

  if isempty(maf0)
    maf0 = Array(Float64,m);
    for j in 1:m
      maf0[j] = mean(geno[~missloc[:,j],j])/2
    end
  end

  if isempty(z0)
    z0 = Array(Int8,(2,n,m))
#        z0[find(missloc),:] = rand(Bernoulli(maf0),sum(missloc)*2);
#        z0[find(geno.==0),:] = 0;
#        z0[find(geno.==2),:] = 1;
#        z0[find(geno.==1),:] = [1 0; 0 1][rand(Bernoulli(0.5),sum(geno.==1))+1,:];
    for j in 1:m
      z0[:,:,j] = pedigree_genosim(ped,maf0[j])[2]';
    end
  end

 zstates = Int8[1 1;
                 0 1;
                 1 0;
                 0 0];

  sireof,damof = pedigree_childfinder(ped);
  sireany = sum(sireof,1)[:] .> 0;
  damany = sum(damof,1)[:] .> 0;

  #initialize gibbs chain
  saveiter = [thin:thin:iter];
  nsave = length(saveiter);
  iter = maximum(saveiter);

  z = Array(Int8,(2,n,m,nsave));
  ϵ = Array(Float64,nsave);
  maf = Array(Float64,m,nsave);
#  badz = falses(iter);

  pk = ones(Float64,4);

  znow = z0;
  ϵnow = ϵ0;
  mafnow = maf0;

  pzflib = hcat([0,1,0,1],  #sire has genotype 0
            [.5,.5,.5,.5], #sire has genotype 1
            [1,0,1,0]); #sire has genotype 2
  pzmlib = hcat([0,0,1,1],  #dam has genotype 0
            [.5,.5,.5,.5], #dam has genotype 1
            [1,1,0,0]); #dam has genotype 2
  pklib = hcat([0,.5,.5,1], #child has mat/paternal major allele
            [1,.5,.5,0]); #child has mat/paternal minor allele

  nma = zeros(Int64,m);
  post = Array(Float64,4);
  nerr = 0;

  for t in 1:iter
    zlik = construct_zlik(ϵnow);

    for i in 1:n
      for j in 1:m

      #calculate genotype probability given base MAF and/or parental genotypes
        if sire[i] > 0 #father in pedigree; sampled from father's alleles
          pzf = pzflib[:,znow[1,sire[i],j] + znow[2,sire[i],j]+1];
        else #father not in pedigree; sample from base population
          pzf = [mafnow[j], 1-mafnow[j], mafnow[j], 1-mafnow[j]];
        end

        if dam[i] > 0 #mother in pedigree; sampled from father's alleles
          pzm = pzmlib[:,znow[1,dam[i],j] + znow[2,dam[i],j]+1];
        else #mother not in pedigree; sample from base population
          pzm = [mafnow[j], mafnow[j], 1-mafnow[j], 1-mafnow[j]];
        end

        pk[:] = 1;
    #calculate offspring genotype probabilities (if there are any)
        if sireany[i]
          for k in find(sireof[:,i])
            for s in 1:4
              pk[s] = pk[s] * pklib[s,znow[1,k,j]+1];
            end
          end
        elseif damany[i]
          for k in find(damof[:,i])
            for s in 1:4
              pk[s] = pk[s] * pklib[s,znow[2,k,j]+1];
            end
          end
        end

      #calculate conditional posteriors

        for s in 1:4
          post[s] = pk[s] * pzm[s] * pzf[s] * zlik[s,geno[i,j]+2];
        end

#     #check for impossible z
#       if sum(post) == 0
#         post[:] = 1;
#         badz[t] = true;
#       end
        #sample new allele
        znow[:,i,j] = zstates[rand(Categorical(post./sum(post))),:];
        #simplify bookeeping by keeping most recent zs represented separately from iteration history

        #update founder minor allele counts
        if ped[i,2]==0
          nma[j] += znow[1,i,j];
        end
        if ped[i,3]==0
          nma[j] += znow[2,i,j];
        end

        #update genotyping error counts
        if ( (znow[1,i,j] + znow[2,i,j]) != geno[i,j] ) & ( geno[i,j] != -1 )
          nerr += 1;
        end

      end
    end

    #update minor allele frequency and allele error frequency
    for j in 1:m
      mafnow[j] = rand( Beta(nma[j] + mafprior[1], nfound - nma[j] + mafprior[2]) );
      nma[j] = 0; #while we're at it reset maf
    end

    #update genotype error probability
    ϵnow = rand( Beta(nerr + ϵprior[1],n*m - nmiss - nerr + ϵprior[2]) );
    nerr = 0; #reset error probability

    nsamp = findfirst(saveiter .== t);
    if nsamp > 0
      maf[:,nsamp] = mafnow;
      ϵ[nsamp] = ϵnow;
      z[:,:,:,nsamp] = znow;
    end

  end
  genosim = reshape(sum(z,1),(n,m,nsave));
  return genosim,maf,ϵ,z

end

function pedigree_childfinder(ped)
  n = size(ped)[1];
  sireof = Array(Bool,(n,n));
  damof = Array(Bool,(n,n));
  for i in 1:n, j in 1:n
    sireof[i,j] = j == ped[i,2];
    damof[i,j] = j == ped[i,3];
  end
  return sireof,damof
end
