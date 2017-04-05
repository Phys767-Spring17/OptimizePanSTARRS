import os as os
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table,join,Column
from decimal import *
import math
from subprocess import Popen


#fits files are composed of various objects, including the quasars sought to be analyzed, 
#they include time, magnitude and error for four bands g,r,i and z. The text files have the
#quasars that are in the fits files. 

#fits files
RA33to35='MD09_DF_RA333to3335.fit'	#3912228 rows and 7186 objects
RA34to45='MD09_DF_RA334to3345.fit'	#4934139 rows and 8809 objects
RA35to55='MD09_DF_RA335to3355.fit'	#2146552 rows and 4561 objects 
RA25to33='MD09_DF_RA3325to333.fit'	#1464642 rows and 3534 objects
RA35to34='MD09_DF_RA3335to334.fit'	#4837498 rows and 8654 objects
RA45to35='MD09_DF_RA3345to335.fit'	#4220240 rows and 7744 objects

#text files containing the quasars being studied in the corresponding fits files (TOTAL 670)
tRA33to35='match_mb_RA333to3335_qso.txt'	#151 quasars
tRA34to45='match_mb_RA334to3345_qso.txt'	#163 quasars
tRA35to55='match_mb_RA335to3355_qso.txt'	#4 quasars
tRA25to33='match_mb_RA3325to333_qso.txt'	#55 quasars
tRA35to34='match_mb_RA3335to334_qso.txt'	#184 quasars
tRA45to35='match_mb_RA3345to335_qso.txt'	#113 quasars

#text files containing the varying quasars being studied in the corresponding fits files (TOTAL 104)
vRA33to35='variables_RA333to3335.txt'	#21 quasars
vRA34to45='variables_RA334to3345.txt'	#24 quasars
vRA35to55='variables_RA335to3355.txt'	#1 quasars
vRA25to33='variables_RA3325to333.txt'	#6 quasars
vRA35to34='variables_RA3335to334.txt'	#33 quasars
vRA45to35='variables_RA3345to335.txt'	#19 quasars

#text files containing the a sample of the varying quasars being studied in the corresponding fits files
sRA33to35='variables_RA333to3335_test.txt'	#2 quasars
sRA34to45='variables_RA334to3345_test.txt'	#2 quasars
sRA35to55='variables_RA335to3355_test.txt'	#1 quasars
sRA25to33='variables_RA3325to333_test.txt'	#1 quasars
sRA35to34='variables_RA3335to334_test.txt'	#2 quasars
sRA45to35='variables_RA3345to335_test.txt'	#2 quasars

#opens fits file and creates an hdulist table with the object number, RA, DEC, time, magnitude values, magnitude errors
def select_fits_columns(filename):
	hdulist=fits.open(filename)
	data=hdulist[1].data
	tbdata=Table([data['BUNDLE'][0],data['RA'][0],data['DEC'][0],data['OBSMJD'][0],data['PSFMAGS'][0],
		data['PSFERRS'][0]],names=('obj','RA','DEC','time','mag','err'))
	return tbdata

#opens a text file composed of strings and returns a list of intergers. ***THE VALUES IN THE TEXT FILE MUST BE INTERGER VALUES***
def textopen(filename):
	txt = open(filename, 'r')
	lines=[]
	for row in txt.readlines():
		lines.append(int(row.split()[0]))
	txt.close()
	return lines

#This function returns a dictionary for ONE band with: q=quasar number; ra=RA; dec=DEC; ti=time (julian days); m=magnitude; e=error; mean=magnitude mean; std=standard deviation of the magnitudes; length=length of rows which correspond number of observations made for that quasar some which are multiple for a given day; con_ti=condensed time that correspond to individual julian day, rather than multiple oservations like the regular time; con_mag=mean magnitude for a day; con_e=error summed in quadrature for each day; and con_lenght=new length the condensed file showing only the number of days the particular quasar was observed for that band
#The function also removes the -999 values that were in the fits file.
#Index values 0, 1, 2, 3 correspond to g, r, i, and z bands.
def create_band_dict(fitsfile,txtfile,index=None):
	table=select_fits_columns(fitsfile)
	txt=textopen(txtfile)
	#lists to create the keys that will go into the dictionary.
	q=[];ra=[];dec=[];ti=[];m=[];e=[];mean=[];std=[];length=[]
	con_time=[];con_mag=[];con_err=[];con_length=[]
	for i in range(len(txt)):
		#lists for the condensed time, magnitude and error values that are created for a given day
		#with multiple observations
		floor_time=[];obs=[];days=[];magn=[];quad_e=[]
		#cut0 selects one quasar in the text file from the objects in the fits file at a time
		cut0=txt[i]==table['obj']
		t=table[cut0]
		time=np.matrix.transpose(t['time'])
		mag=np.matrix.transpose(t['mag'])
		err=np.matrix.transpose(t['err'])
		#cut1 removes the rows with -999 magnitudes
		cut1=mag[index]!=-999.0
		q.append(t['obj'][0])
		ra.append(t['RA'][cut1])
		dec.append(t['DEC'][cut1])
		ti.append(time[index][cut1])
		m.append(mag[index][cut1])
		e.append(err[index][cut1])
		mean.append(np.mean(mag[index][cut1]))
		std.append(np.std(mag[index][cut1]))
		length.append(len(mag[index][cut1]))
		#section that creates the condensed time, magnitude, error, and row length to use one set 
		#of values for each day, rather than multiple value for a given day. Row length correspond
		#to the number of days data was collected for that quasar.
		for i in range(len(time[index][cut1])):
			#creates a list of the floored days (ex a=1.1, 1.7 becomes a=1.0, 1.0) 
			floor_time.append(math.floor(time[index][cut1][i]))
			#creates a list of observation days (ex b=0.0, 1.0, 1.0, 1.0, 2.0 becomes b=0.0, 1.0, 2.0)
			if math.floor(time[index][cut1][i]) not in obs:
				obs.append(math.floor(time[index][cut1][i]))
		#uses the list of observations to select a particular day in oder to average the time,mag,err	
		for i in range(len(np.array(obs))):
			cut2=np.array(obs)[i]==np.array(floor_time)
			days.append(np.mean(time[index][cut1][cut2]))
			magn.append(np.mean(mag[index][cut1][cut2]))
			error=err[index][cut1][cut2]
			quad_e.append((np.sqrt(np.sum(error*error)))/len(error))
		con_time.append(days)
		con_mag.append(magn)
		con_err.append(quad_e)
		con_length.append(len(days))
	d={'QUASAR':np.array(q),'RA':np.array(ra),'DEC':np.array(dec),'TIME':np.array(ti),'MAG':np.array(m),
		'ERR':np.array(e),'MEAN':np.array(mean),'STD':np.array(std),'LENGTH':np.array(length),
		'TIME_con':np.array(con_time),'MAG_con':np.array(con_mag),'ERR_con':np.array(con_err),
		'LENGTH_con':np.array(con_length)}
	return d

#This function uses the 'create_band_dict' function to return a dictionary sorted by bands g,r,i,z with all the information  in the 'create_band_dict' function
def create_quasar_dict(fitsfile,txtfile):
	d={}
	d['g']=create_band_dict(fitsfile,txtfile,0)
	d['r']=create_band_dict(fitsfile,txtfile,1)
	d['i']=create_band_dict(fitsfile,txtfile,2)
	d['z']=create_band_dict(fitsfile,txtfile,3)
	return d

#This function creates ascii files with unique names that contain the time, magnitude and error in a table format for each quasar in the set of text and fits files. It also removes the created file. To keep the ascii file, comment out 'os.remove commands' that you'd like. ***BEWARE CREATING ASCII FILES FOR ALL FILES WILL RESULT 2680 FILES. SEE TEXT FILE COMMENT FOR NUMBER OF QUASARS IN EACH FILE BEFORE CREATING ASCII FILES***. It also creates a text file 'name' of the various ascii files created by this function without the .dat 
def create_quasar_savetxt(fitsfile,txtfile,name="TEST"):
	d=create_quasar_dict(fitsfile,txtfile)
	txt = open(name, "w")
	for i in range(len(d['g']['QUASAR'])):
		t_g=Table([d['g']['TIME_con'][i],d['g']['MAG_con'][i],d['g']['ERR_con'][i]],names=('TIME','MAG','ERR'))	
		t_r=Table([d['r']['TIME_con'][i],d['r']['MAG_con'][i],d['r']['ERR_con'][i]],names=('TIME','MAG','ERR'))
		t_i=Table([d['i']['TIME_con'][i],d['i']['MAG_con'][i],d['i']['ERR_con'][i]],names=('TIME','MAG','ERR'))
		t_z=Table([d['z']['TIME_con'][i],d['z']['MAG_con'][i],d['z']['ERR_con'][i]],names=('TIME','MAG','ERR'))
		txt.write('RA%d_'%d['g']['RA'][0][i]+'DEC%.3f_'%d['g']['DEC'][0][i]+'_QSO'+str(d['g']['QUASAR'][i])+'_g\n')
		txt.write('RA%d_'%d['g']['RA'][0][i]+'DEC%.3f_'%d['g']['DEC'][0][i]+'_QSO'+str(d['g']['QUASAR'][i])+'_r\n')
		txt.write('RA%d_'%d['g']['RA'][0][i]+'DEC%.3f_'%d['g']['DEC'][0][i]+'_QSO'+str(d['g']['QUASAR'][i])+'_i\n')
		txt.write('RA%d_'%d['g']['RA'][0][i]+'DEC%.3f_'%d['g']['DEC'][0][i]+'_QSO'+str(d['g']['QUASAR'][i])+'_z\n')	
		np.savetxt('RA%d_'%d['g']['RA'][0][i]+'DEC%.3f_'%d['g']['DEC'][0][i]+'_QSO'+
			str(d['g']['QUASAR'][i])+'_g.dat',t_g,fmt=('%.32f','%.32f','%.32f'))
		np.savetxt('RA%d_'%d['r']['RA'][0][i]+'DEC%.3f_'%d['r']['DEC'][0][i]+'_QSO'+
			str(d['g']['QUASAR'][i])+'_r.dat',t_r,fmt=('%.32f','%.32f','%.32f'))
		np.savetxt('RA%d_'%d['i']['RA'][0][i]+'DEC%.3f_'%d['i']['DEC'][0][i]+'_QSO'+
			str(d['g']['QUASAR'][i])+'_i.dat',t_i,fmt=('%.32f','%.32f','%.32f'))
		np.savetxt('RA%d_'%d['z']['RA'][0][i]+'DEC%.3f_'%d['z']['DEC'][0][i]+'_QSO'+
			str(d['g']['QUASAR'][i])+'_z.dat',t_z,fmt=('%.32f','%.32f','%.32f'))
		#when not commented out, removes files from operating system
		'''os.remove('RA%d_'%d['g']['RA'][0][i]+'DEC%.3f_'%d['g']['DEC'][0][i]+'_QSO'+
			str(d['g']['QUASAR'][i])+'_g.dat')
		os.remove('RA%d_'%d['r']['RA'][0][i]+'DEC%.3f_'%d['r']['DEC'][0][i]+'_QSO'+
			str(d['r']['QUASAR'][i])+'_r.dat')
		os.remove('RA%d_'%d['i']['RA'][0][i]+'DEC%.3f_'%d['i']['DEC'][0][i]+'_QSO'+
			str(d['i']['QUASAR'][i])+'_i.dat')
		os.remove('RA%d_'%d['z']['RA'][0][i]+'DEC%.3f_'%d['z']['DEC'][0][i]+'_QSO'+
			str(d['z']['QUASAR'][i])+'_z.dat')'''
	txt.close()

#create_quasar_savetxt(RA33to35,sRA33to35,"sample_RA33to35_ascii_list")
#create_quasar_savetxt(RA34to45,sRA34to45,"sample_RA34to45_ascii_list")
#create_quasar_savetxt(RA35to55,sRA35to55,"sample_RA35to55_ascii_list")
#create_quasar_savetxt(RA25to33,sRA25to33,"sample_RA25to33_ascii_list")
#create_quasar_savetxt(RA35to34,sRA35to34,"sample_RA35to34_ascii_list")
#create_quasar_savetxt(RA45to35,sRA45to35,"sample_RA45to35_ascii_list")
'''
create_quasar_savetxt(RA33to35,vRA33to35,"vary_RA33to35_ascii_list.txt")
create_quasar_savetxt(RA34to45,vRA34to45,"vary_RA34to45_ascii_list.txt")
create_quasar_savetxt(RA35to55,vRA35to55,"vary_RA35to55_ascii_list.txt")
create_quasar_savetxt(RA25to33,vRA25to33,"vary_RA25to33_ascii_list.txt")
create_quasar_savetxt(RA35to34,vRA35to34,"vary_RA35to34_ascii_list.txt")
create_quasar_savetxt(RA45to35,vRA45to35,"vary_RA45to35_ascii_list.txt")'''

#takes in a light curve and values of tau and sigma, processes it through cole's Zoghbi code and returns a tuple composed of tau & sigma entered and likelihood value of the pair. ***LIGHT CURVE MUST NOT HAVE HEADERS***
def cole(filename, tau, sigma):
    fd = os.popen("%s %s %g %g" % ('./Zoghbi2013Cole', filename, tau, sigma))
    lines = fd.readlines()
    fd.close()
    lh = float(lines[-1])
    #print lh
    return (tau,sigma,lh)

#creates a table with index, log tau, tau, sigma and likelihood by entering a light curve, beginning log tau, ending log tau, number of values that need to be examined within the log tau range (ex if log tau is 1-4 and you enter 5, this will divide the range of log tau to 1, 1.6, 2.2, 2.8, 3.4, 4), beginning sigma, ending sigma and number of values that need to be examined within the sigma range (similar to the log tau divisions)
def likelihood_table(filename,logtau_start,logtau_stop,tau_div,sig_start,sig_stop,sig_div,save_filename,save='yes'):
	logtau=np.linspace(logtau_start,logtau_stop,tau_div)
	tau=10**logtau	#cole function accepts values in tau rather than log tau
	sig=np.linspace(sig_start,sig_stop,sig_div)
    ntau=len(tau) #number of tau's
    nsig=len(sig) #number of sigma's
	#creates the tau and sigma pairs to feed into cole's code
	param=[]
	for i in range(ntau):
		for j in range(nsig):
			param.append([tau[i],sig[j]])
	#creates tau, sigma, and likelihoods matrix from coles code given the range of tau and sigma pairs
	TSL=[]
	for i in range(len(param)):
		TSL.append(cole(filename,param[i][0],param[i][1]))
	#creates table with index, log tau (days), tau (days), sigma [mag/(days^0.5)], and likelihood
	T, S, L =np.array(zip(*TSL)) #creates arrays for tau, sigma and likelihod returned by cole function
	index=np.arange(0,len(T),1)
	t=Table([index,np.round(np.log10(T),3),np.round(T,2),np.round(S,4),L],
		names=('INDEX','LOG_TAU','TAU','SIGMA','LIKELIHOOD'),dtype=('i8','f64','f64','f64','f64'))
	if save=='yes':
		ascii.write(t,save_filename+'_tau'+str(tau_div)+'_sig'+str(sig_div)+'_table.dat')
	if save=='no':
		ascii.write(t,save_filename+'_tau'+str(tau_div)+'_sig'+str(sig_div)+'_table.dat')		
		os.remove(save_filename+'_tau'+str(tau_div)+'_sig'+str(sig_div)+'_table.dat')
	return t

#creates likelihood tables from files like those made by create_quasar_savetxt function
#lcs=light curves; lc=light curve
def likelihood(filename,logtau_start,logtau_stop,tau_div,sig_start,sig_stop,sig_div,save='yes'):
	txt = open(filename, 'r')
	lcs=[]
	for row in txt.readlines():
		lcs.append(row.split()[0])
	txt.close()

	for lc in range(len(lcs)):
		print lcs[lc]
		likelihood_table(lcs[lc]+'.dat',logtau_start,logtau_stop,tau_div,sig_start,sig_stop,sig_div,lcs[lc],save)
	print ''
'''
likelihood("sample_RA33to35_ascii_list.txt",1,4,100,0.001,0.021,100,'yes')
likelihood("sample_RA34to45_ascii_list.txt",1,4,100,0.001,0.021,100,'yes')
likelihood("sample_RA35to55_ascii_list.txt",1,4,100,0.001,0.021,100,'yes')
likelihood("sample_RA25to33_ascii_list.txt",1,4,100,0.001,0.021,100,'yes')
likelihood("sample_RA35to34_ascii_list.txt",1,4,100,0.001,0.021,100,'yes')
likelihood("sample_RA45to35_ascii_list.txt",1,4,100,0.001,0.021,100,'yes')

likelihood("vary_RA33to35_ascii_list.txt",1,4,10,0.001,0.021,10,'yes')
likelihood("vary_RA34to45_ascii_list.txt",1,4,10,0.001,0.021,10,'yes')
likelihood("vary_RA35to55_ascii_list.txt",1,4,10,0.001,0.021,10,'yes')
likelihood("vary_RA25to33_ascii_list.txt",1,4,10,0.001,0.021,10,'yes')
likelihood("vary_RA35to34_ascii_list.txt",1,4,10,0.001,0.021,10,'yes')
likelihood("vary_RA45to35_ascii_list.txt",1,4,10,0.001,0.021,10,'yes')

likelihood("vary_RA33to35_ascii_list.txt",1,4,25,0.001,0.021,25,'yes')
likelihood("vary_RA34to45_ascii_list.txt",1,4,25,0.001,0.021,25,'yes')
likelihood("vary_RA35to55_ascii_list.txt",1,4,25,0.001,0.021,25,'yes')
likelihood("vary_RA25to33_ascii_list.txt",1,4,25,0.001,0.021,25,'yes')
likelihood("vary_RA35to34_ascii_list.txt",1,4,25,0.001,0.021,25,'yes')
likelihood("vary_RA45to35_ascii_list.txt",1,4,25,0.001,0.021,25,'yes')

likelihood("vary_RA33to35_ascii_list.txt",1,4,50,0.001,0.021,50,'yes')
likelihood("vary_RA34to45_ascii_list.txt",1,4,50,0.001,0.021,50,'yes')
likelihood("vary_RA35to55_ascii_list.txt",1,4,50,0.001,0.021,50,'yes')
likelihood("vary_RA25to33_ascii_list.txt",1,4,50,0.001,0.021,50,'yes')
likelihood("vary_RA35to34_ascii_list.txt",1,4,50,0.001,0.021,50,'yes')
likelihood("vary_RA45to35_ascii_list.txt",1,4,50,0.001,0.021,50,'yes')'''

#opens ascii file (.dat) and creates a table with column names.
#If ascii file has no header then skiprows should equal 0, if header exists then skiprows should be 1.
def create_ascii_likelihood_table(filename,skiprows):
	data = np.loadtxt(filename,skiprows=skiprows)
	col0=[];col1=[];col2=[];col3=[];col4=[]
	for i in range(len(data)):
		col0.append(data[i][0])
		col1.append(data[i][1])
		col2.append(data[i][2])
		col3.append(data[i][3])
		col4.append(data[i][4])
	col0=np.array(col0);col1=np.array(col1);col2=np.array(col2);col3=np.array(col3);col4=np.array(col4)
	tbdata=Table([col0,col1, col2,col3, col4], 
		names=('INDEX','LOG_TAU','TAU','SIGMA','LIKELIHOOD'),dtype=('i8','f64','f64','f64','f64'))
	return tbdata



