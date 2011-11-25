%example Matlab script for EF-1alpha rate shift visualization

%load the rate shift data
[rawdata,branchdata,sitedata,score]=read_stochmap_out2('../examples/EF1A.shift.out');
branchscores=branchdata(:,2)-branchdata(:,3);%conditional-prior
sitescores=sitedata(:,2)-sitedata(:,3);%conditional-prior

%generate a Dendroscope file with edges coloured by conditional-prior rate shifts
colortreeedges2(branchdata(:,5),branchscores,'../examples/EF1alpha_Feb11_fullnames.dendro','../examples/EF1alpha_shiftcolour.dendro',[min(branchscores),max(branchscores)],[],[],0,0);

%generate b-values to colour by sign of conditional-prior rate shifts on
%structure in Pymol
fname='../examples/1IJF_EFalphaplusyeast.aln';
gapchar='-';
target='A_Aerper';
[gapcvectplus,gapcvectminus,seq,gapsigns]=addgaps2(fname,target,sitescores,gapchar,sign(sitescores));
dlmwrite('../examples/siteshiftminus.txt',gapcvectminus);
dlmwrite('../examples/siteshiftplus.txt',gapcvectplus);


