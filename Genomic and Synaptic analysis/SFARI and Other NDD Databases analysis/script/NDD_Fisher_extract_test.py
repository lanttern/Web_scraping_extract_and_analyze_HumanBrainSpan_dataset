import scipy as sp
import scipy.stats
"""
NDD_data = input("Please enter contingency data: ")
print "You data are: " 
print NDD_data
(Oddsratio, p_value) = sp.stats.fisher_exact(NDD_data)
print "Results of Fisher exact test: "
print "Oddsratio: " + str(Oddsratio) + " p_value: " + str(p_value)
"""
MET_pulldown_genome = [[7 - 1,65],[77,29923]]
full_genome = [[16 - 1,1237],[77,29923]]
sub_genome = [[62 -1 ,322],[616,29384]]
full_sub_genome = [[62 -1,322],[105,1148]]
"""
# for SynDB stastics
MET_pulldown_synapse = [[9 - 1,22],[276,2973]]
full_synapse = [[70-1,561],[276,2973]]
sub_synapse = [[41-1,151],[276,2973]]
full_sub_synapse = [[41-1,151],[70,561]]
synapse_genome = [[276-1,2973],[616,29384]]
"""
# for hPSD database
MET_pulldown_synapse = [[5-1,26],[30,3219]]
full_synapse = [[12-1,619],[30,3219]]
sub_synapse = [[16-1,49],[54,694]]
full_sub_synapse = [[16-1,49],[25-1,158]]
synapse_genome = [[54-1,694],[616,29384]]
data = [MET_pulldown_genome, full_genome,sub_genome,full_sub_genome,MET_pulldown_synapse, full_synapse,sub_synapse,full_sub_synapse, synapse_genome]
for num in data:
    print num, sp.stats.fisher_exact(num)
