#!/usr/bin/env python3
"""
WE ARE STILL PROTING TO PYTHON3!
This is  a simplified version of dna_param.py
It does not need VMD, but is easy to use, for instance to change DNA sequence in a structure.

It has one more function, change_dna_seq_in_pdb 

3DNA SHOULD BE INSTALLED AND IN PATH.
X3DNA ENVIRONMENT VARIABLE SET

CURVES FUNCTIONS IN DEVELOPMENT

"""
import os
import subprocess


import pandas as pd
import re

import uuid

# from scipy.spatial import KDTree
# from scipy.spatial import cKDTree
import numpy as np
from collections import OrderedDict

__author__="Alexey Shaytan"

TEMP='/tmp/dna_tools'
os.system("mkdir -p %s"%TEMP)
P_X3DNA_DIR=''
P_X3DNA_analyze='analyze'
P_X3DNA_find_pair='find_pair'
P_X3DNA_rebuild='rebuild'
P_X3DNA_x3dna_utils='x3dna_utils'


P_CURVES='/Users/alexeyshaytan/bin/Cur+'
P_CURVES_LIB='/Users/alexeyshaytan/soft/curves+/standard'


#os.environ['X3DNA']=P_X3DNA_DIR

def X3DNA_find_pair(DNA_pdb):
	""" Runs find_pair program from X3DNA and returns a path to unique file with defined pairs

	This is needed to supply this path to 3DNA_analyze,
	so that is will know what pairs to expect in data.
	The need to do find pairs separately arises if we analyze
	MD trajectory, where base pairing may vary, but we still
	want the same number of data points in output.

	Parameters
	----------
	DNA_pdb - path to pdb

	Return
	----------
	returns a unique string - which is the file name of find_pair outfile
	in the TEMP directory.
	"""

	#At first we need to makup a couple of unique file names
	unique=str(uuid.uuid4())
	pdb = unique+'.pdb'
	outf = unique


	print("Writing coords to "+pdb)
	os.system("cp "+DNA_pdb+' '+TEMP+'/'+pdb)

	cmd=P_X3DNA_find_pair+' '+pdb+' '+outf
    #print(cmd)
	p = subprocess.Popen(cmd,shell=True,cwd=TEMP,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True)
	# wait for the process to terminate
	out, err = p.communicate()
	#errcode = p.returncode
	print("OUT:",out)
	print("ERR:",err)
	return(outf)


def X3DNA_analyze(DNA_pdb,ref_fp_id):
	"""Performs the analysis using X3DNA

	Parameters
	----------
	DNA_pdb - see above
	ref_fp_id - this is output id from X3DNA_find_pair function,
	which was obtained for the structure that will be considered as a reference
	to determine which bases are paired.

	Return
	--------
	PANDAS data frame of the following format:
	Rows are numbered sequentially and correspond to the output of X3DNA (user has to check how X3DNA handled numbering in specific structure).
	This output established the numbering of base-pairs (usually this coincide with the numbering of the first strand).
	The only newance is with numbering of sugar and backbone parameters:
	in the X3DNA output both stands are output sequentially in 5'-3' manner.
	Again this might depend on the format of the supplied PDB - so the user is advised to check,
	that with his file behavior is the same.
	We assign to all column names subscrip _1 for the first stand, and _2 to the second.
	Moreover, we renumber the second stand in a 3'-5' orientation (IMPORTANT)!!!!!
	All the names of the returned parameters correspond to their names in X3DNA.
	Additional columns:
	Pairing - 1 if X3DNA sees a base pair there with respect to reference (even if if is non standart pairing), 0 if not.
	x,y,z - the centers of reference frames of individual base pairs.
	BPnum - numer of base pair from 1 to N
	"""

	#Now we still run find_pairs on this frame to check if any base pairing was lost
	#and simultaneously to output pdb
	cur_fp_id=X3DNA_find_pair(DNA_pdb)
	pdb = cur_fp_id+'.pdb'
	inp = cur_fp_id

	#we have to do substitution in ref_fp_id file and copy it
	#so it will process new file using original base pair information
	#sed 's/pattern/replacement/g'

	cmd='sed "s/'+ref_fp_id+'/'+cur_fp_id+'/g" '+ref_fp_id+'>'+cur_fp_id+'.fr' #fr - dreived from reference
	p = subprocess.Popen(cmd,shell=True,cwd=TEMP,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True)
	out, err = p.communicate()
	print('OUT:',out,err)
	
	#Now we can run X3DNA_analyze
	cmd=P_X3DNA_analyze+' '+cur_fp_id+'.fr'
	p = subprocess.Popen(cmd,shell=True,cwd=TEMP,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True)
	out, err = p.communicate()
	print('OUT:',out,err)

################################################
#Extract base pairing, bp centers, bp params, bp step params
################################################

	#####Base pairing (might be some got unpaired with respect to reference)
	
	#####Extract centers of base pairs from ref_frames.dat
	df_rf=parse_ref_frames(TEMP+'/'+'ref_frames.dat')
	# print(df_rf)
	###Extract base pair and base pair step parameters
	df_bp=parse_bases_param(TEMP+'/'+'bp_step.par')
	# print(df_bp)
	###Extract base pairing by comparing reference and current
	df_pairing=check_pairing(TEMP+'/'+ref_fp_id,TEMP+'/'+cur_fp_id)
	# print df_pairing
################################################
	#Special call to X3DNA_analyze that will get sugar and backbone params
	#This call stangly deletes some files from previous call
	#So we need to extract base-pair and ref frames info before
	cmd=P_X3DNA_analyze+' -t=backbone.tor '+pdb
	# print cmd
	p = subprocess.Popen(cmd,shell=True,cwd=TEMP,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# wait for the process to terminate
	out, err = p.communicate()
	errcode = p.returncode
	print('OUT:',out,err)
####################################
##Now let's get torsion parameters
	df_tor=parse_tor_param(TEMP+'/backbone.tor')
	# print(df_tor)
	#Now we concatenate all the data frames
	df_res=pd.concat([df_rf,df_bp,df_pairing,df_tor],axis=1)
	df_res['BPnum']=range(1,len(df_res)+1)
	df_res=df_res.reset_index(drop=True)
	return(df_res)


def X3DNA_analyze_bp_step(DNA_pdb,ref_fp_id):
	"""Performs the analysis using X3DNA and output only bp_step

	Parameters
	----------
	DNA_atomsel - see above
	ref_fp_id - this is output id from X3DNA_find_pair function,
	which was obtained for the structure that will be considered as a reference
	to determine which bases are paired.

	Return
	--------
	PANDAS data frame of the following format:
	Rows are numbered sequentially and correspond to the output of X3DNA (user has to check how X3DNA handled numbering in specific structure).
	This output established the numbering of base-pairs (usually this coincide with the numbering of the first strand).
	
	All the names of the returned parameters correspond to their names in X3DNA.
	Additional columns:
	Pairing - 1 if X3DNA sees a base pair there with respect to reference (even if if is non standart pairing), 0 if not.
	x,y,z - the centers of reference frames of individual base pairs.
	BPnum - numer of base pair from 1 to N
	"""

	#Now we still run find_pairs on this frame to check if any base pairing was lost
	#and simultaneously to output pdb
	cur_fp_id=X3DNA_find_pair(DNA_pdb)
	pdb = cur_fp_id+'.pdb'
	inp = cur_fp_id

	#we have to do substitution in ref_fp_id file and copy it
	#so it will process new file using original base pair information
	#sed 's/pattern/replacement/g'

	cmd='sed "s/'+ref_fp_id+'/'+cur_fp_id+'/g" '+ref_fp_id+'>'+cur_fp_id+'.fr' #fr - dreived from reference
	p = subprocess.Popen(cmd,shell=True,cwd=TEMP,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = p.communicate()
	print('OUT:'+out+err)
	
	#Now we can run X3DNA_analyze
	cmd=P_X3DNA_analyze+' '+cur_fp_id+'.fr'
	p = subprocess.Popen(cmd,shell=True,cwd=TEMP,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = p.communicate()
	print('OUT:'+out+err)

################################################
#Extract base pairing, bp centers, bp params, bp step params
################################################


	###Extract base pair and base pair step parameters
	df_bp=parse_bases_param(TEMP+'/'+'bp_step.par')
	# print(df_bp)

################################################

####################################

	df_res=df_bp
	df_res['BPnum']=range(1,len(df_res)+1)
	df_res=df_res.reset_index(drop=True)
	return(df_res)




def parse_ref_frames(file):
	"""
	Parses ref_frames.dat file from X3DNA output
	and returns a PANDAS data frame

	"""
	print("Processing ", file)
	#let's get the data from file
	# let's strip white spaces at first and replace with tabs
	tmpfile=file.replace('.dat','.tmp')
	with open(file,'r') as f:
		with open(tmpfile,'w') as tf:
			f.readline()
			f.readline()
			for line in f:
				tf.write('\t'.join(line.partition('#')[0].split())+'\n')
				next(f,'')
				next(f,'')
				next(f,'')
				next(f,'')

		#now we will read data frames, 
		df=pd.read_csv(tmpfile, sep=u'\t',skiprows=0,header=None)
		df.columns=['x','y','z']

	# bdf.to_csv('../analysis_data/dna_bl_big_df.csv')
	return(df)

def parse_bases_param(file):
	"""
	Parse bp_step.par file as output by X3DNA
	Get a data frame with base and base-pair step
	parameters

	Note that in each line of the data frame there are parameters
	for some base pair and base pair step that preceeds (!) this base-pair.
	I.e. in the first line no base pair step parameters are specified!

	offest - offset for DNA numbering.
	"""
	print("Processing ", file)
	#let's get the data from file
	# let's strip white spaces at first and replace with tabs
	tmpfile=file+'tmp'
	with open(file,'r') as f:
		with open(tmpfile,'w') as tf:
			for line in f:
				tf.write('\t'.join(line.split())+'\n')



	df=pd.read_csv(tmpfile, sep=u'\t',skiprows=2)
	#The first column is of the type A-T or A+T, the latter indicating something ...
	# df=df.rename(columns={'#':'BP'})
#Paired or unpaired
	# for i in range(len(df)):
	# 	if(i in bp_list):
	# 		df.iloc[i,0]=1
	# 	else:
	# 		df.iloc[i,0]=0
	# print df.head()
	df.iloc[0,7:]=np.nan # since base bair step params are not defined in first row
	
	#This is changed to add base pair names 29 June 2015 - hope this would not break the VMD nucleosome analysis scripts.
	# df=df.drop(df.columns[0], axis=1)
	new_columns = df.columns.values
	new_columns[0]='BPname'
	df.columns=new_columns

	return(df)



def check_pairing(ref,cur):
	"""
	Functions compairs two files output by 3DNA find_pair
	for same structures in different conformations
	and looks what base pairs are present/lost
	in cur with respect to reference

	This function is not well tested!!!
	"""
	bp_list_ref=list()

	with open(ref,'r') as inf:
		for line in inf:
			if(re.search('\.\.\.\.>C:\.*(-?\d+)_:',line)):
				bp_list_ref.append(int(re.search('\.\.\.\.>C:\.*(-?\d+)_:',line).group(1)))
	print("Reference BP list")
	print(bp_list_ref)

	bp_list_cur=list()

	with open(cur,'r') as inf:
		for line in inf:
			if(re.search('\.\.\.\.>C:\.*(-?\d+)_:',line)):
				bp_list_cur.append(int(re.search('\.\.\.\.>C:\.*(-?\d+)_:',line).group(1)))
	print("Current BP list")
	print(bp_list_cur)
#Let's construct data frame by comparing
	df_pairing = pd.DataFrame(columns=['Pairing'])
	for i in range(len(bp_list_ref)):
		if(bp_list_ref[i] in bp_list_cur):
			df_pairing=df_pairing.append({'Pairing':1}, ignore_index=True)
		else:
			df_pairing=df_pairing.append({'Pairing':0}, ignore_index=True)
	return(df_pairing)

def parse_tor_param(file):
	"""
	Parse torsion parameters as returned by X3DNA (-t option) (tor-file)

	"""

	torfile=file
	angfile=file+'.bba'
	puckfile=file+'.pck'
	with open(torfile,'r') as f:
		f.readline()
		with open(angfile,'w') as tf:
			for line in f:
				if('****' not in line):
						plist=line.split()
						if(len(plist)>5):
							if((len(plist)!=12) and (plist[0]!='base')):
								plist.insert(3,'no')
						tf.write('\t'.join(plist)+'\n')
				else:
					break
		for line in f:
			if('*****' in line): break
		with open(puckfile,'w') as tf:
			for line in f:
				tf.write('\t'.join(line.split())+'\n')

	dfa=pd.read_csv(angfile, sep=u'\t',skiprows=19)
	dfp=pd.read_csv(puckfile, sep=u'\t',skiprows=18)
	#Now the idea is to slice this data frames into strands
	dfa1=dfa.iloc[0:int(len(dfa)/2),1:]
	dfa2=dfa.iloc[int(len(dfa)/2):,1:]
	n=dfa1.columns
	n1=[]
	n2=[]
	for i in n:
		n1.append(i+'_1')
		n2.append(i+'_2')
	dfa1.columns=n1
	dfa2.columns=n2
	dfa1.index=range(len(dfa1))
	dfa2.index=range(len(dfa2)-1,-1,-1)
	dfa_new=pd.concat([dfa1,dfa2],axis=1)

	dfp1=dfp.iloc[0:int(len(dfp)/2),1:]
	dfp2=dfp.iloc[int(len(dfp)/2):,1:]
	n=dfp1.columns
	n1=[]
	n2=[]
	for i in n:
		n1.append(i+'_1')
		n2.append(i+'_2')
	dfp1.columns=n1
	dfp2.columns=n2
	dfp1.index=range(len(dfp1))
	dfp2.index=range(len(dfp2)-1,-1,-1)
	dfp_new=pd.concat([dfp1,dfp2],axis=1)

	df_new=pd.concat([dfa_new,dfp_new],axis=1)

	return(df_new)

def CURVES_analyze(DNA_pdb,length):
	"""Performs the analysis using Curves+

	Parameters
	----------
	DNA_atomsel - see above
	length - length of one DNA strand.
	Returns
	-------
	Curently returns groove params.
	"""


	unique=str(uuid.uuid4())
	pdb = unique+'.pdb'

	# print("Writing coords to "+pdb)
	# DNA_atomsel.write('pdb',TEMP+'/'+pdb)
	os.system('cp '+DNA_pdb+' '+TEMP+'/'+pdb)

	#Now we can run CURVES+
	cmd=P_CURVES+' <<!\n &inp file=%s, lis=%s,\n lib=%s\n &end\n2 1 -1 0 0\n1:%d\n%d:%d\n!'%(pdb,pdb,P_CURVES_LIB,length,length*2,length+1)
	print(cmd)
	p = subprocess.Popen(cmd,shell=True,cwd=TEMP,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = p.communicate()
	print('OUT:'+out+err)
	#Now let's parse Curves output
	lis=TEMP+'/'+pdb+'.lis'
	return(parse_lis(lis))



def parse_lis(file):
	"""
	Parses CURVES+ lis file to get groove parameters
	"""

	tmpfile=file.replace('.lis','.tmpb')
	tmpfile2=file.replace('.lis','.tmpg')

	#Let's parse the backbone
	with open(file,'r') as f: 
		with open(tmpfile,'w') as t:
			with open(tmpfile2,'w') as t2:
				#loop till we find
				#  (D) Backbone Parameters
				for line in f:
					if(re.search("\(D\)",line)): break
				# Strand 1     Alpha  Beta   Gamma  Delta  Epsil  Zeta   Chi    Phase  Ampli  Puckr  
				for line in f:
					if(re.search("Strand",line)):
						t.write('Strand\t'+'Resid\t'+'\t'.join(line.split()[2:])+'\n')
						f.next() # skip one blank line
						break
				for line in f:
					if(re.search("\d+\)\s+[ATGC]",line)):
						t.write('1'+'\t'+'\t'.join(line.split()[2:]).replace('----','')+'\n')
					else: break
				for line in f:
					if(re.search("Strand",line)):
						f.next() # skip one blank line
						break
				for line in f:
					if(re.search("\d+\)\s+[ATGC]",line)):
						t.write('2'+'\t'+'\t'.join(line.split()[2:]).replace('----','')+'\n')
					else: break
				for line in f:
					if(re.search("\(E\)",line)): break
				for line in f:
					if(re.search("Level",line)):
						t2.write('\t'.join(line.split())+'\n')
						f.next() # skip one blank line
						break
				for line in f:
					if(re.search("\s+\d+",line)):
						t2.write(line[0:8].strip()+'\t'+line[16:22].strip()+'\t'+line[23:30].strip()+'\t'+line[31:38].strip()+'\t'+line[39:46].strip()+'\n')
					else: break
					
	#Uncomment this if you want backbone data, its here ready!!!
	# df_bb=pd.read_csv(tmpfile, sep=u'\t')

	df_gr=pd.read_csv(tmpfile2, sep=u'\t')

	return(df_gr)
	# df_bb.to_csv('../analysis_data/dna_bb_param_cryst.csv')

	# df_gr.to_csv('../analysis_data/dna_gr_param_cryst.csv')


###Here goes the rebuilding mechanism.
#Here is how it works from command line
#find_pair 1kx5.pdb 1kx5.inp
#analyze 1kx5.inp
##change bp_step.par
#x3dna_utils cp_std BDNA
#rebuild -atomic bp_step.par nucl_new.pdb

def gen_bp_step(data_frame,new_seq=None):
	"""
	Generates the bp_step.par file based on a data frame (might be the output of X3DNA_analyze or X3DNA_analyze_bp_step)
	Returns a file name

	Optionally we can place a new sequnce here in new_seq ['A','T',..].
	"""

	unique=str(uuid.uuid4())
	par = unique+'.par'
	full_path=TEMP+'/'+par

	new_df=data_frame[['BPname','Shear','Stretch','Stagger','Buckle','Prop-Tw','Opening','Shift','Slide','Rise','Tilt','Roll','Twist']]

	#put new sequence
	if(new_seq):
		d={'A':'A-T','T':'T-A','G':'G-C','C':'C-G'}
		comp_seq=map(lambda x: d[x],new_seq)
		new_df['BPname']=pd.Series(comp_seq,index=new_df.index)

	# new_columns = new_df.columns.values
	# new_columns[0]='#'
	# new_df.columns=new_columns
	new_df=new_df.fillna(0.0)


	f=open(full_path+'tmp','w')
	# f=open(par,'w')

	f.write(" %d # base-pairs\n"%len(new_df.index))
	f.write('   0 # ***local base-pair & step parameters***\n')
	f.write('#        Shear    Stretch   Stagger   Buckle   Prop-Tw   Opening     Shift     Slide     Rise      Tilt      Roll      Twist\n')
	# new_df.to_string(f,index=False,na_rep='0.000',float_format='{:>10.3f}'.format)
	new_df.to_string(f,index=False,header=False,justify='left',col_space=0,na_rep='0.000',formatters=['{:4s}'.format,'{:>9.3f}'.format,'{:>9.3f}'.format,'{:>9.3f}'.format,'{:>9.3f}'.format,'{:>9.3f}'.format,'{:>9.3f}'.format,'{:>9.3f}'.format,'{:>9.3f}'.format,'{:>9.3f}'.format,'{:>9.3f}'.format,'{:>9.3f}'.format,'{:>9.3f}'.format])
	f.close()
	#Now megahack - we need to strip one space from every row.
	with open(full_path+'tmp','rb') as f:
		with open(full_path,'wb') as f2:
			for line in f:
				f2.write(line[0:])


	return(full_path)


def build_dna(data_frame,pdbfile,new_seq=None):
	"""
	Runs gen_bp_step and then runs the rebuilding of DNA via X3DNA
	and output to pdbfile
	"""
	par_fname=gen_bp_step(data_frame,new_seq)

	cmd=P_X3DNA_x3dna_utils+' cp_std BDNA'
	p = subprocess.Popen(cmd,shell=True,cwd=TEMP,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = p.communicate()
	print('OUT:',out,err)

	cmd=P_X3DNA_rebuild+' -atomic '+par_fname+' '+par_fname+'.pdb'
	p = subprocess.Popen(cmd,shell=True,cwd=TEMP,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = p.communicate()
	print('OUT:',out,err)

	os.system('mv '+par_fname+'.pdb '+pdbfile)

def change_dna_seq_in_pdb(DNA_pdb,NEW_pdb,new_seq):
	"""
	Self-explanatory
	"""

	data_frame=X3DNA_analyze_bp_step(DNA_pdb,X3DNA_find_pair(DNA_pdb))
	build_dna(data_frame,NEW_pdb,new_seq)


	

if __name__ == '__main__':
	print("Kuku")
	X3DNA_find_pair('x')
	help(atomsel)

