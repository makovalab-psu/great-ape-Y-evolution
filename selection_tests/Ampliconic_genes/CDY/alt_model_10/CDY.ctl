seqfile  = CDY.phy 	         * sequence data file name
treefile = CDY.tree10.nwk	 	 * tree structure file name
outfile  = CDY.mlc 	                 * main result file name

 noisy = 9     * 0,1,2,3,9: how much rubbish on the screen
verbose = 1     * 1: detailed output, 0: concise output
runmode = 0     * 0: user tree;  1: semi-automatic;  2: automatic
                * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

  seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
    clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
    model = 2   * models for codons:
                * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
  NSsites = 0   * 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete;
                * 4:freqs; 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:3normal
    icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
    
fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
    kappa = 2   * initial or fixed kappa
fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate
    omega = 2   * initial or fixed omega, for codons or codon-based AAs

    fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = .0  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * different alphas for genes
        ncatG = 4   * # of categories in the dG or AdG models of rates

        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)
       method = 0   * 0: simultaneous; 1: one branch at a time
