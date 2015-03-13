function pedigree_genosim(ped,maf)
  n = size(ped,1)
  geno = Array(Int64,n,2)
  dset = Int64[]

#base population probabilities
  bpop = Distributions.Bernoulli(maf)

# Identify individuals with only genotyped and unknown/founder
  while length(dset) < n
    mset = find( (indexin(ped[:,2],dset) .> 0) | (ped[:,2] .==0) & (indexin(ped[:,3],dset) .> 0) | (ped[:,3] .==0) )
    mset = setdiff(mset,dset)

    for i in mset
      sire = ped[i,2]
      dam = ped[i,3]

      geno[i,1] = (sire > 0) ? geno[sire,rand(1:2)] : rand(bpop)
      geno[i,2] = (dam > 0) ? geno[dam,rand(1:2)] : rand(bpop)

    end

    dset = append!(dset,mset)
  end

  return sum(geno,2),geno
end

function pedigree_pedsim(init,ngen,popcap)

  gensize = min(init * 2.^[0:(ngen-1)],popcap);
  genbound = cumsum(gensize);

  ped = Array(Int64,(last(genbound),3));
  sire = Array(Int64,last(genbound));
  dam = Array(Int64,last(genbound));

  grng = [1,genbound[1]];
  sire[grng[1]:grng[2]] = 0;
  dam[grng[1]:grng[2]] = 0;

  for i in 2:ngen
    oldrng = grng[1:2];
    grng[1] = genbound[i-1]+1;
    grng[2] = genbound[i];

    #even = female, odd = male
    sire[grng[1]:grng[2]] = rand(oldrng[1]:2:oldrng[2],gensize[i])
    dam[grng[1]:grng[2]] = rand( (oldrng[1]+1):2:oldrng[2],gensize[i]);
  end

  ped[:,1] = [1:last(genbound)];
  ped[:,2] = sire;
  ped[:,3] = dam;
  return ped
end
