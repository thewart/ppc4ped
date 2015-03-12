function pedigree_genosim(ped,maf)
  n = size(ped,1)
  geno = Array(Int64,n,2)
  dset = Int64[]

#base population probabilities
  bpop = Bernoulli(maf)

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

  return sum(geno,2)
end

