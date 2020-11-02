# only consider lens below 99% percentile
get_fragment_length_dist = function(gene_models, rowID, num_thre = 10,
                             bam_path, strandmode = 0, quant){
  print("start estimating mean and sd of fragment length ...")
  flag = 1

  # B. Fant: Replacing a mclapply with a standard for loop
  # mclapply had a nasty habit of hanging and ended up anaylizing way too many reads and in the wrong way
  # W/ the v34 GENCODE annotation you end up analyzing 472 genes, very reasonable for a "for" loop, even for R that hates them
  # It's one of those weird examples where good coding does you a disservice, so we replace it w/ cookie-cutter stuff
  # NOTE: in this state it is not parallelized nor is it particularly able to be parallelized. However, after benchmark it takes ~150-250 sec tu run (=4 min tops)
  # so we won't bother abt it too much, it remains plenty fast.
  # 
  fglens = c()
  ngenes_counting = 0

  sel = sample(1:length(rowID), length(rowID))
  rowID = rowID[sel]

  for (rowid in rowID){
      if (length(fglens < 10000)){
        cgene = gene_models[[rowid]]
        readTxcoords = AIDE:::get_reads(cgene, num_thre = num_thre, bam_path, strandmode = strandmode)
        if (!is.null(readTxcoords)){
          fglen_cgene = unlist(readTxcoords$ends - readTxcoords$starts + 1)
          fglens = append(fglens, fglen_cgene)
          ngenes_counting = ngenes_counting + 1
        }
      } else{
        break
      }
    }

  print(paste("number of genes used:", ngenes_counting))
  if(length(fglens) < 2000){
    print("do not have enough reads for estimation ...")
    print("using default values ...")
    return(list(mean = 250, sd = 80, cutoff = NA))
  }

  fglens = fglens[fglens < quantile(fglens, quant)]
  print(paste("number of fragments:", length(fglens)))
  print(paste("longest fragments:", max(fglens)))
  print(paste("shortest fragments:", min(fglens)))

  mean = mean(fglens)
  sd = sd(fglens)
  cutoff = range(fglens, na.rm = TRUE)
  return(list(mean = mean, sd = sd, cutoff = cutoff, flag = flag))
}
