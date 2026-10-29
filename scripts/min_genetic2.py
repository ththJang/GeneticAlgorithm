#!/usr/bin/python
import os, sys
import random
import copy
#import commands
#from collections import deque
from operator import itemgetter


############################################################################### function
#def evaluation(lst): 
## input format : a list of variables comprising force field parameter.
## output : force error summation for each case at the end and the list itself
## evaluate the force error according to the input
## do the MD calculation with "run 0"
## note that ref_force should be already exist in the same directory
#	num=len(lst)
#	os.system('~/program/FIND_min_emdext/emdext 1 1 1 -1 in.lammps qm_0/solute.cube HW OW SOLUTE 1.0 1.0 0.0 > lammps.equil.out') # MD run 0 
#	# take the force output in QM/MM calculation in *.lammpstrj
#	os.system('./lammps_force.sh')

#	(exitstatus, error)=commands.getstatusoutput('./error_calc.py ref_force lammps_force')
#	lst.append(error)
#	return lst
def mutation(ori_val,rangeofx,extend):
# input format : original value, a list including the range of domain of elements.
# output : value after mutated or conserved
# have a completely different value of input.
# in case of that random() has lower than the set value.
	do_mutate=0.1
	new_ran=0
	if random.random()<do_mutate:
		new_ran=random.random()
#		ori_val=0.9*rangeofx*(1-int(extend))+new_ran*rangeofx*extend
		ori_val=rangeofx*(1-extend+2*extend*new_ran)
	return ori_val
	
def next_gen(chromosomes):
# input format : a list of all chromosomes
# output : half of original list + crossovered list
# sorting and store only probable lists
	population=len(chromosomes) # number of consisting lists(population) it is multiples of 4.
	kill=population*0.5
	new=sorted(chromosomes, key=itemgetter(-1)) # sorting in the increasing order of force error
	for i in range(len(new[0])):
		print("%f " % new[0][i],end = "")
	print()
#	print new[0] # this is the case for the minimum , format : variables w/ error
	for i in range(int(kill)): # half of the chromosomes remained.
		del new[-1]
	for i in range(len(new)):   # get rid of the errors
		del new[i][-1]
	auto_selected=selection(new) 
	lst=crossover(auto_selected) # recover the population
	return lst
def crossover(lst):
# input format : half list of the original set of chromosomes
# output : get back to full set
# crossover occurs randomly rate of 50% 
#	NAtoms=len(lst[0])/3 # len(lst[0])/x , x for buckingham
	NAtoms=len(lst[0]) # # of QMLMPTYPE
	newchromosome=len(lst)
	for i in range(len(lst)):
		lst.append(copy.deepcopy(lst[i]))
	for i in range(int(len(lst)/4)):
		for k in range(int(NAtoms)):
			alpha=random.random()
			beta=random.random()            # do crossover since only one parameter for each pair
			newx1=beta*lst[newchromosome+2*i][k]+(1-beta)*lst[(newchromosome+1)+2*i][k]
			newx2=(1-beta)*lst[newchromosome+2*i][k]+beta*lst[(newchromosome+1)+2*i][k]
			lst[newchromosome+2*i][k]=copy.deepcopy(newx1)
			lst[(newchromosome+1)+2*i][k]=copy.deepcopy(newx2)
#			if alpha<0.333: # since buckingham potential is specified by three variables
#				j=0
#			elif alpha<0.666:
#				j=1
#			else:
#				j=2
#			newx1=beta*lst[6+2*i][j+3*k]+(1-beta)*lst[7+2*i][j+3*k] # 6, 7 for crossover in case of that the population is 12
#			newx2=(1-beta)*lst[6+2*i][j+3*k]+beta*lst[7+2*i][j+3*k] # crossover between new chromosomes
#			lst[6+2*i][j+3*k]=newx1
#			lst[7+2*i][j+3*k]=newx2

	return lst
def selection(chromosomes):
# input format : half list of the original set of chromosomes
# output : according to the force error, most probable chromosome among the top half will be survived statistically.
# to be automatically selected
	sum=0
	ele_sum=0
	travel=len(chromosomes)
	prob=[]
	new=[]
	for i in range(1,travel+1):
		sum+=i
	sum+=travel # prob revision
	ele_sum=travel # prob revision
	for i in range(travel):
		ele_sum+=travel-i
		prob.append(float(ele_sum)/sum)
	for j in range(travel):
		i=0
		beta=random.random()
		while prob[i]<beta:
			i+=1
		new.append(chromosomes[i])
	return new
###################################################################   main
###################################################################   using chromosomes_w_error making new generation
chromosomes_file=sys.argv[1]
chromosomes=[]
de=[]
de_float=[]
for line in open(chromosomes_file,'r'):         # parse input file 12*(#shot)X(#ATOMTYPE+1error)
	de=line.split()
	for i in range(len(de)):
		de_float.append(float(de[i]))
	chromosomes.append(de_float)
	de_float=[]
child=next_gen(chromosomes)
#mutation_range=[0.3,2,12.886,0.3,2.4,12.886,0.3,2.65,12.886,0.3,2.36,12.886]			# system dependent
#mutation_range=[1000,1000,500,500,2000]			# system dependent
#mutation_range = [1 for k in range(len(de)-1)] # len -1 for error
#mutation_range = [1, 12, 3, 20, -30, 13, 2, 1, 10, 10, 10, -5]
#mutation_range = [0.9, 10.2511, 2.0214, 15, -22.5, 12.7434, 1.1345, 0.7861, 3.5071, 1.1693, 3.6922, -0.8921]
mutation_range = [0.916237, 14.736512, 2.073908, 10.256803, -43.310242, 13.056735, 1.132739, 0.449305, 3.556502, 0.869538, 0.0, 0.932314]
atom_order=len(mutation_range)
for i in range(len(child)):
	for j in range(int(len(child[0]))): # child[0]/x x is a number of variables for each atom i.e. x=3 for buckingham
		child[i][j]=mutation(child[i][j],mutation_range[(j)%atom_order],1.0)

#for i in range(len(child)):
#	for j in range(int(len(child[0])/3)): # child[0]/x x is a number of variables for each atom i.e. x=3 for buckingham
#		child[i][3*j]=mutation(child[i][3*j],mutation_range[(3*j)%atom_order],1)
#	for j in range(int(len(child[0])/3)):
#		child[i][3*j+1]=mutation(child[i][3*j+1],mutation_range[(3*j+1)%atom_order],0.1)
#	for j in range(int(len(child[0])/3)):
#		child[i][3*j+2]=mutation(child[i][3*j+2],mutation_range[(3*j+2)%atom_order],0.1)
for i in range(len(child)):
	for j in range(len(child[i])):
		print("%f " % child[i][j],end = "")
	print()
#for i in range(len(child)):
#	child[i]=mutation(child[i],[5,10,5,10])
#for i in range(len(child)):
#	print child[i]
#chromosomes=[[1,3],[3,6],[5,9],[7,12]]
#child=next_gen(chromosomes)
#for i in range(len(child)):
#	mutation(child[i],[3,5])
#print child
