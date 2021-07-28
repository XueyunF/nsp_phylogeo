# Library importation
import sys
import os
import getopt
import numpy
from numpy import array
import moments
import pylab
import matplotlib
from datetime import datetime
import Optimize_Functions
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D
import Model_2pop_folded
import Model_2pop_unfolded


#Help function
def usage():
    """ Function for help """
    print("# This script allow you to test different demographic models on your genomic data using moments on python 3 \n"+
    "# It uses the optimisation routine propose by Portick et al, 2017 for dadi and adapted to moments by Momigliano et al., 2020 \n\n"+
    "# This is an exemple of the most complete command line :\n"+
    "# python Script_moments_model_optimisaton.py -f SNP_file -y Pop1 -x Pop2 -p 30,30 -m IM_b -s True -r 1 -z -n 6 \n\n"+
    "# that will performed the optimisaton for the IM model with a folded SFS (-s True) with singleton masked (-z) obtained from the file named SNP_file.data \n\n"+
    "# -h --help : Display the help you are looking at.\n"+
    "# -f SNPs input data file with .data as extention \n"+
    "# -y --population1 : Take the name of the first population in the sfs (y-axis)\n"+
    "# -x --population2 : Take the name of the second population in the sfs (x-axis)\n"+
    "# -p --projection : tak two number for the size of the projecition"+
    "# -m --model : model to test, which are found in Model_2pop and comprise five basic (SI, IM, AM, SC, 2EP) and their extention:\n"+
    "# SC_b for SC with bottleneck SC_ae for ancestral expension SC_2N for Hill-Robertson effect SC_2M for heterogeneous gene flow and combinaton of thereof  \n"+
    "# For more information on models see paper by  Momigliano et al. 2020\n"+
    "# -s --folded : True or False: True if SFS is folded, or False if SFS unflolded using an outgroup \n"+
    "# -r --replicat : replicat number \n" +
    "# -z : mask the singletons (no argument to provide) \n"+
    "# -n --rounds:  number of rounds of optimization it must be 2 to 6 \n"+
    "########################## Enjoy ###########################")
    return()
	

def takearg(argv):
   fs_fil_name = ''
   popname1=''
   popname2=''
   projection=''
   model=''
   folded=''
   replicat='1'
   rounds='3'
   masked=False
   try:
      opts, args = getopt.getopt(argv,"hf:x:y:p:m:s:r:n:z",["fs_fil_name=","popname1=","popname2=","projection=","model=","folded=","replicat=","rounds=","masked"])
   except getopt.GetoptError:
    print("ERROR : INCORRECT INPUT PROVIDED.\nDescription of the script usage bellow: \n\n")
    print(usage())
    sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         usage()
         sys.exit()
      elif opt in ("-f", "--fs_fil_name"):
         fs_fil_name = arg
      elif opt in ("-x", "--popname1"):
         popname1 = arg
      elif opt in ("-y", "--popname2"):
         popname2 = arg
      elif opt in ("-p", "--projection"):
         projection = arg.split(',')
      elif opt in ("-m", "--model"):
         model = arg
      elif opt in ("-r", "--replicat"):
         replicat = arg
      elif opt in ("-s", "--folded"):
         folded = arg
      elif opt in ("-n", "--rounds"):
         rounds = arg
      elif opt in ("-z", "--masked"):
         masked = True

   print ('Input file is ', fs_fil_name)
   print ('popname1 is ', popname1)
   print ('popname2 is ', popname2)
   print ('projection of', projection,' is folded? ',folded, 'and is masked? ', masked)
   print ('Replicate of ', model,' number ',replicat, 'nrounds', rounds)

   return(fs_fil_name, popname1, popname2, projection,model, folded,replicat, rounds,masked)
if __name__ == "__takearg__":
   takearg(sys.argv[1:])

###########################
###      actual script         #####
###########################

# Load parameters
fs_fil_name, popname1, popname2, projection,model, folded,replicat, rounds, masked= takearg(sys.argv[1:]) 

print ("nrounds",rounds,'\n')
# Load the data
#Create python dictionary from snps file
input_file=fs_fil_name+'.data'
dd = moments.Misc.make_data_dict(input_file)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=[popname2, popname1]
print (pop_ids)
#**************
#projection sizes, in ALLELES not individuals
proj=list(map(int,projection))

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
if folded == 'True':
    data = moments.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)
if folded == 'False':
    data = moments.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = True)

#fig = pylab.figure(1)
#moments.Plotting.plot_single_2d_sfs(data, vmin=0.1)
#fig.savefig(fs_fil_name+'.svg')

if masked:
	data.mask[1,0] = True
	data.mask[0,1] = True


#Summary statistic of the SFS
print ("\nSFS summary \n")
print ("sample sizes", data.sample_sizes)
sfs_sum = numpy.around(data.S(), 2)
print ("Sum of SFS = ", sfs_sum)
Fst = numpy.around(data.Fst(),4)
print ("FST =", Fst, '\n')

print ("nrounds",rounds,'\n')

#Optimisaton set up
if folded == 'True':
    prefix = fs_fil_name+'_folded_'+replicat
if folded == 'False':
    prefix = fs_fil_name+'_unfolded_'+replicat

if rounds == '2':
    reps = [30,30] 
    maxiters = [20,30]
    folds = [3,2]
	
if rounds == '3':
    reps = [30,20,20] 
    maxiters = [10,20,30]
    folds = [3,2,1]
	
if rounds == '4':
    reps = [30,20,20,20]
    maxiters = [10,20,20,30]
    folds = [3,2,2,1]

if rounds == '5':
    reps = [30,20,20,20,30]
    maxiters = [10,20,20,20,30]
    folds = [3,2,2,2,1]

if rounds == '6':
    reps = [30,20,20,20,20,30]
    maxiters = [10,20,20,20,30,30]
    folds = [3,2,2,2,2,1]

#################
###   Model SI   ###
#################


if model == 'SI':
    if folded == "True":
        p_labels = "nu1, nu2, T1"
        upper = [100,100,10]
        lower = [1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI", Model_2pop_folded.SI,rounds, 3, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if folded == "False":
        p_labels = "nu1, nu2, T1, O"
        upper = [100,100,10, 0.999]
        lower = [1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI", Model_2pop_unfolded.SI,rounds, 4, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")


if model == 'SI_NeC':
    if folded == "True":
        p_labels = "nu1, nu2, nu1b, nu2b, T1, T2"
        upper = [100,100,100,100,10,10]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_NeC", Model_2pop_folded.SI_NeC,rounds, 6, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if folded == "False":
        p_labels = "nu1, nu2,nu1b, nu2b, T1, T2, O"
        upper = [100,100,100,100,10,10, 0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SI_NeC", Model_2pop_unfolded.SI_NeC,rounds, 7, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
		
		
#################
###   Model IM   ###
#################

if model == 'IM':
    if folded == "True":
        p_labels = "nu1, nu2, T1, m12, m21"
        upper = [100,100,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM", Model_2pop_folded.IM,rounds, 5, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if folded == "False":
        p_labels = "nu1, nu2, T1, m12, m21,O"
        upper = [100,100,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM", Model_2pop_unfolded.IM,rounds, 6, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model == 'IM_NeC':
    if folded == "True":
        p_labels = "nu1, nu2, nu1b, nu2b, T1, T2, m12, m21"
        upper = [100,100,100,100,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_NeC", Model_2pop_folded.IM_NeC,rounds, 8, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if folded == "False":
        p_labels = "nu1, nu2,  nu1b, nu2b, T1, T2, m12, m21, O"
        upper = [100,100,100,100,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "IM_NeC", Model_2pop_unfolded.IM_NeC,rounds, 9, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")


#################
###  Model SC   ###
#################

if model =="SC":
    if folded == "True":
        p_labels = "nu1, nu2, T1, T2,m12, m21"
        upper = [100,100,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC", Model_2pop_folded.SC,rounds, 6, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if folded == "False":
        p_labels = "nu1, nu2, T1, T2,m12, m21,O"
        upper = [100,100,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC", Model_2pop_unfolded.SC,rounds, 7, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model =="SC_NeC":
    if folded == "True":
        p_labels = "nu1, nu2, nu1b, nu2b,T1, T2,m12, m21"
        upper = [100,100,100,100,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_NeC", Model_2pop_folded.SC_NeC,rounds, 8, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if folded == "False":
        p_labels = "nu1, nu2,nu1b,nu2b, T1, T2,m12, m21,O"
        upper = [100,100,100,100,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "SC_NeC", Model_2pop_unfolded.SC_NeC,rounds, 9, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

#################
###  Model AM   ###
#################

if model =="AM":
    if folded == "True":
        p_labels = "nu1, nu2, T1, T2,m12, m21"
        upper = [100,100,10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM", Model_2pop_folded.AM,rounds, 6, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if folded == "False":       
        p_labels = "nu1, nu2, T1, T2,m12, m21,O"
        upper = [100,100,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM", Model_2pop_unfolded.AM,rounds, 7, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model =="AM_NeC":
    if folded == "True":
        p_labels = "nu1, nu2, nu1b, nu2b, T1, T2,m12, m21"
        upper = [100,100,100, 100, 10,10,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_NeC", Model_2pop_folded.AM_NeC,rounds, 8, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if folded == "False":       
        p_labels = "nu1, nu2, nu1b, nu2b, T1, T2,m12, m21,O"
        upper = [100,100,100,100,10,10,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "AM_NeC", Model_2pop_unfolded.AM_NeC,rounds, 9, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

#################
###  Model two_ep  ###
#################
if model =="two_ep":
    if folded == "True":
        p_labels = "nu1, nu2, T1, T2,m12_0,m21_0,m12,m21"
        upper = [100,100,10,10,200,200,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep", Model_2pop_folded.two_ep,rounds, 8, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if folded == "False":       
        p_labels = "nu1, nu2, T1, T2,m12_0,m21_0,m12,m21,O"
        upper = [100,100,10,10,200,200,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep", Model_2pop_unfolded.two_ep,rounds, 9, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")

if model =="two_ep_NeC":
    if folded == "True":
        p_labels = "nu1, nu2, nu1b, nu2b, T1, T2,m12_0,m21_0,m12,m21"
        upper = [100,100,100,100,10,10,200,200,200,200]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_NeC", Model_2pop_folded.two_ep_NeC,rounds, 10, data_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")
    if folded == "False":       
        p_labels = "nu1, nu2, nu1b, nu2b, T1, T2,m12_0,m21_0,m12,m21,O"
        upper = [100,100,100,100,10,10,200,200,200,200,0.999]
        lower = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]
        Optimize_Functions.Optimize_Routine(data, prefix, "two_ep_NeC", Model_2pop_unfolded.two_ep_NeC,rounds, 11, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
        if done: print ("\n" + model + "folded = "+folded+ " is done\n")


