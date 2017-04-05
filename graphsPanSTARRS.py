import os as os
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table,join,Column
import matplotlib.pyplot as plt
from decimal import *
import math
from subprocess import Popen
from likelihood import *

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

def likelihood_stats(filename,print_stats='no'):
	t=create_ascii_likelihood_table(filename,1)
	cut_Lmax=t['LIKELIHOOD']==np.amax(t['LIKELIHOOD'])
	cut_conf=t['LIKELIHOOD']>=(np.amax(t['LIKELIHOOD'])-(2.3/2))
	L_max=np.amax(t['LIKELIHOOD']) 		#maximum likelihood
	L_min=np.amin(t['LIKELIHOOD'])		#minimum likelihood
	L_conf=t[cut_conf]					#confidence of likelihood for 1 sigma given 2 parameters
	L_logtau=t['LOG_TAU'][cut_Lmax][0]	#max Likelihood log tau
	L_tau=t['TAU'][cut_Lmax][0]			#max Likelihood tau
	L_sig=t['SIGMA'][cut_Lmax][0]		#max Likelihood sigma
	logT_min=np.amin(L_conf['LOG_TAU'])
	logT_max=np.amax(L_conf['LOG_TAU'])
	Tmin=np.amin(L_conf['TAU'])
	Tmax=np.amax(L_conf['TAU'])
	Smin=np.amin(L_conf['SIGMA'])
	Smax=np.amax(L_conf['SIGMA'])
	if print_stats=='yes':
		print ''
		print filename
		print ''
		print 'Maximum Likelihood:..........................	'+str(L_max)	
		print 'Max Likelihood log tau (days):.............	'+str(L_logtau)
		print 'Max Likelihood tau (days):..................	'+str(L_tau)
		print 'Max Likelihood sigma [mag/(days^0.5)]:'+str(L_sig)
		print '------------------------------------------'
		print 'Minimum Likelihood: '+str(L_min)
		print ''
		print 'Confidence Minimum Log tau (days): '+str(logT_min)
		print 'Confidence Maximum Log tau (days): '+str(logT_max)
		print 'Confidence Range Log tau (days): '+str(L_logtau)+'-'+str(L_logtau-logT_min)+'; '+str(L_logtau)+'+'+str(logT_max-L_logtau)
		print ''
		print 'Confidence Minimum tau (days): '+str(Tmin)
		print 'Confidence Maximum tau (days): '+str(Tmax)
		print 'Confidence Range tau (days): '+str(L_tau)+'-'+str(L_tau-Tmin)+'; '+str(L_tau)+'+'+str(Tmax-L_tau)
		print ''
		print 'Confidence Minimum sigma [mag/(days^0.5)]: '+str(Smin)
		print 'Confidence Maximum sigma [mag/(days^0.5)]: '+str(Smax)
		print 'Confidence Range sigma [mag/(days^0.5)]: '+str(L_sig)+'-'+str(L_sig-Smin)+'; '+str(L_sig)+'+'+str(Smax-L_sig)
		print ''
	return (L_max,L_logtau,L_tau,L_sig,logT_min,logT_max,Tmin,Tmax,Smin,Smax)


def likelihood_band_stats_table(filename,qsofile,print_stats='no'):
	txt = open(filename, 'r')
	tbs=[]
	for row in txt.readlines():
		tbs.append(row.split()[0])
	txt.close()
	
	qso=textopen(qsofile)
	
	L_stats=[]
	for tb in range(len(tbs)):
		L_stats.append(likelihood_stats(tbs[tb]+'_tau50_sig50_table.dat',print_stats))
	L,L_logtau,L_tau,L_sig,logT_min,logT_max,Tmin,Tmax,Smin,Smax=np.array(zip(*L_stats))

	t=Table([np.array(qso),L,L_logtau,L_tau,L_sig,logT_min,logT_max,Tmin,Tmax,Smin,Smax],
		names=('QUASAR','LIKELIHOOD','LOG_TAU_LH','TAU_LH','SIGMA_LH','LOG_TAU_Cmax','LOG_TAU_Cmin', 			'TAU_Cmax','TAU_Cmin','SIGMA_Cmax','SIGMA_Cmin'), 
		dtype=('i8','f64','f64','f64','f64','f64','f64','f64','f64','f64','f64'))
	return t

def create_entire_likelihood_band_stats_table(file1,txt1,file2,txt2,file3,txt3,file4,txt4,file5,txt5,file6,txt6,print_stats='no'):
	t1=likelihood_band_stats_table(file1,txt1,print_stats)
	t2=likelihood_band_stats_table(file2,txt2,print_stats)
	t3=likelihood_band_stats_table(file3,txt3,print_stats)
	t4=likelihood_band_stats_table(file4,txt4,print_stats)
	t5=likelihood_band_stats_table(file5,txt5,print_stats)
	t6=likelihood_band_stats_table(file6,txt6,print_stats)
	for i in range(len(t2)):
		t1.add_row(t2[i])
	for i in range(len(t3)):
		t1.add_row(t3[i])
	for i in range(len(t4)):
		t1.add_row(t4[i])
	for i in range(len(t5)):
		t1.add_row(t5[i])
	for i in range(len(t6)):
		t1.add_row(t6[i])
	index=np.arange(0,len(t1),1)
	t1.add_column(Column(index,name='INDEX'), index=0)
	return t1

def create_entire_likelihood_stats_dic():
	d={}
	d['g']=create_entire_likelihood_band_stats_table(agRA33to35, vRA33to35, agRA34to45, vRA34to45, agRA35to55, vRA35to55, agRA25to33, vRA25to33, agRA35to34, vRA35to34, agRA45to35, vRA45to35)
	d['r']=create_entire_likelihood_band_stats_table(arRA33to35, vRA33to35, arRA34to45, vRA34to45, arRA35to55, vRA35to55, arRA25to33, vRA25to33, arRA35to34, vRA35to34, arRA45to35, vRA45to35)
	d['i']=create_entire_likelihood_band_stats_table(aiRA33to35, vRA33to35, aiRA34to45, vRA34to45, aiRA35to55, vRA35to55, aiRA25to33, vRA25to33, aiRA35to34, vRA35to34, aiRA45to35, vRA45to35)
	d['z']=create_entire_likelihood_band_stats_table(azRA33to35, vRA33to35, azRA34to45, vRA34to45, azRA35to55, vRA35to55, azRA25to33, vRA25to33, azRA35to34, vRA35to34, azRA45to35, vRA45to35)
	return d

def Likelihood_scatgraph(filename,vmin=5,vmax=5,titleop='no',title=None):
	t=create_ascii_likelihood_table(filename,1)
	Lmax=np.amax(t['LIKELIHOOD'])
	Lmin=np.amin(t['LIKELIHOOD'])
	points=np.arange(0,len(t['LIKELIHOOD']),1)
	plt.scatter(points,t['LIKELIHOOD'],c='b', marker='.')
	plt.xlabel(r'$\tau$'+' '+'$(days)$'+' '+'$and$'+' '+'$\sigma$'+' '+'$(mag/day^{1/2})$'+' '+'$Parameter$'+' '+'$Pairs$',fontsize=25,fontweight='bold')
	plt.ylabel('$Likelihood$',fontsize=25,fontweight='bold')
	plt.xticks(np.arange(0,len(points),1000),points[0:len(points):1000],fontsize=16)
	plt.yticks(fontsize=16)
	plt.axis([-500,10500,Lmin-vmin,Lmax+vmax],fontsize=30,fontweight='bold')
	plt.legend(loc=4)
	plt.show()

#Likelihood_scatgraph(QSO4565_g)
'''
Likelihood_scatgraph(QSO616_g,100,100)
Likelihood_scatgraph(QSO616_r,100,100)
Likelihood_scatgraph(QSO616_i,100,100)
Likelihood_scatgraph(QSO616_z,100,100)
'''
#Likelihood_scatgraph(QSO557_g,100,100)
'''
plt.hist(qso557_g100['LIKELIHOOD'], bins=15, range=(80,110))
plt.axis([112,78,0,1600])
plt.show()

plt.hist(qso557_g10['LIKELIHOOD'], bins=8, range=(90,110))
plt.axis([112,88,0,20])
plt.show()

plt.hist(qso557_g25['LIKELIHOOD'], bins=10, range=(80,110))
plt.axis([112,78,0,150])
plt.show()'''

def Likelihood_heatgraph(filename,logtau_start,logtau_stop,tau_div,sig_start,sig_stop,sig_div,xtick_step=1,ytick_step=1):
	t=create_ascii_likelihood_table(filename,1)
	Lmax=np.amax(t['LIKELIHOOD'])
	logtau=np.round(np.linspace(logtau_start,logtau_stop,tau_div),3)
	tau=10**logtau
	sig=np.round(np.linspace(sig_start,sig_stop,sig_div),4)
	ntau=len(tau)
	nsig=len(sig)
	contour_con=[Lmax-(2.3/2),Lmax]
	contour_max=[Lmax-0.001,Lmax]
	contour=[np.round(Lmax-5000), np.round(Lmax-2000), np.round(Lmax-1000), np.round(Lmax-500), np.round(Lmax-200), 		np.round(Lmax-100), np.round(Lmax-50), np.round(Lmax-20), np.round(Lmax-10), Lmax-(11.83/2), 			  			Lmax-(6.18/2)]
	LH=np.reshape(t['LIKELIHOOD'],(ntau,nsig))
	plt.figure()	
	cpc = plt.contour(LH,10,levels=contour_con,colors='k',linewidths=3)
	cpm = plt.contour(LH,10,levels=contour_max,colors='white',linewidths=5)
	cp = plt.contour(LH,10,levels=contour,colors='gray',linewidths=2)
	plt.clabel(cpc, inline=False,fontsize=15,colors='white',fontweight='bold')
	plt.clabel(cp, inline=False,fontsize=13.5,colors='k',fontweight='bold')
	plt.imshow(LH, interpolation='nearest',cmap=plt.get_cmap('jet'),origin='lower')
	plt.xticks(np.arange(0,nsig,xtick_step),sig[0:nsig:xtick_step],fontsize=30)
	plt.yticks(np.arange(0,ntau,ytick_step),np.round(logtau,1)[0:ntau:ytick_step],fontsize=30)
	plt.xlabel('$\sigma$'+' $(mag/day^{1/2})$',fontsize=50,fontweight='bold')
	plt.ylabel('$log$ '+r'$\tau$'+' $(days)$',fontsize=50,fontweight='bold')
	plt.gca().xaxis.set_label_coords(.5,-.05)
	plt.gca().yaxis.set_label_coords(-.08,.5)
	plt.colorbar()
	plt.show()
'''
Likelihood_heatgraph(QSO616_g25,1,4,25,0.001,0.021,25)
Likelihood_heatgraph(QSO616_r25,1,4,25,0.001,0.021,25)
Likelihood_heatgraph(QSO616_i25,1,4,25,0.001,0.021,25)
Likelihood_heatgraph(QSO616_z25,1,4,25,0.001,0.021,25)

Likelihood_heatgraph(QSO616_g10,1,4,10,0.001,0.021,10)
Likelihood_heatgraph(QSO616_r10,1,4,10,0.001,0.021,10)
Likelihood_heatgraph(QSO616_i10,1,4,10,0.001,0.021,10)
Likelihood_heatgraph(QSO616_z10,1,4,10,0.001,0.021,10)

Likelihood_heatgraph(QSO644_g,1,4,100,0.001,0.021,100,15,10)
Likelihood_heatgraph(QSO644_r,1,4,100,0.001,0.021,100,10,10)
Likelihood_heatgraph(QSO644_i,1,4,100,0.001,0.021,100,10,10)
Likelihood_heatgraph(QSO644_z,1,4,100,0.001,0.021,100,10,10)

Likelihood_heatgraph(QSO616_g,1,4,100,0.001,0.021,100,10,10)
Likelihood_heatgraph(QSO616_r,1,4,100,0.001,0.021,100,10,10)
Likelihood_heatgraph(QSO616_i,1,4,100,0.001,0.021,100,10,10)
Likelihood_heatgraph(QSO616_z,1,4,100,0.001,0.021,100,10,10)

Likelihood_heatgraph(QSO4565_g,1,4,100,0.001,0.021,100)
Likelihood_heatgraph(QSO4565_r,1,4,100,0.001,0.021,100)
Likelihood_heatgraph(QSO4565_i,1,4,100,0.001,0.021,100)
Likelihood_heatgraph(QSO4565_z,1,4,100,0.001,0.021,100)

plt.figure()
#cp = plt.contour(LH,10,levels=contour,colors='gray',linewidths=2)
cp = plt.contour(LH,10,colors='gray',linewidths=2)
cpc = plt.contour(LH,10,levels=contour_con,colors='k',linewidths=3)
cpm = plt.contour(LH,10,levels=contour_max,colors='k',linewidths=20)
plt.clabel(cp, inline=True,fontsize=13.5,colors='k',fontweight='bold')
plt.clabel(cpc, inline=False,fontsize=15,colors='white',fontweight='bold')
plt.imshow(LH, interpolation='nearest',cmap=plt.get_cmap('jet'),origin='lower')
plt.xticks(np.arange(0,nsig,10),sig[0:nsig:10],fontsize=17)
plt.yticks(np.arange(0,ntau,10),np.round(logtau,1)[0:ntau:10],fontsize=17)
plt.xlabel('$\sigma$'+' $(mag/day^{1/2})$',fontsize=27,fontweight='bold')
plt.ylabel('$log$ '+r'$\tau$'+' $(days)$',fontsize=27,fontweight='bold')
plt.title('RA334_DEC0.613__QSO4565_g',y=1.05,fontsize=17,fontweight='bold')
plt.gca().xaxis.set_label_coords(.5,-.09)
plt.gca().yaxis.set_label_coords(-.09,.5)
plt.colorbar()
plt.show()

#likelihood_stats(QSO616_g)

cut_max=q4565_g['LIKELIHOOD']==np.amax(q4565_g['LIKELIHOOD'])
cut_min=q4565_g['LIKELIHOOD']==np.amin(q4565_g['LIKELIHOOD'])
cut_conf=q4565_g['LIKELIHOOD']>=(np.amax(q4565_g['LIKELIHOOD'])-(2.3/2))
#cut_conf=np.logical_and(q4565_g['LIKELIHOOD']>(np.amax(q4565_g['LIKELIHOOD'])-0.5),q4565_g['LIKELIHOOD']<=(np.amax(q4565_g['LIKELIHOOD'])))
#cut_conf=np.logical_and(q4565_g['LIKELIHOOD']>(np.amax(q4565_g['LIKELIHOOD'])-0.5),q4565_g['LIKELIHOOD']<(np.amax(q4565_g['LIKELIHOOD'])+0.5))

conf_q4565_g=q4565_g[cut_conf]

print q4565_g[cut_max]
print ''
print q4565_g[cut_min]
print ''
print conf_q4565_g[0:15]
print conf_q4565_g[15:30]
print ''

points=np.arange(0,len(q4565_g['LIKELIHOOD']),1)

print q4565_g['LIKELIHOOD']

plt.plot(points,q4565_g['LIKELIHOOD'],c='b', marker='.')
plt.xlabel(r'$\tau$'+' $  (days)  and$ '+'$\sigma$'+' $(mag/day^{1/2})    Parameter    Pairs$',fontsize=25,fontweight='bold')
plt.ylabel('$Likelihood$',fontsize=25,fontweight='bold')
plt.xticks(np.arange(0,len(points),1000),points[0:len(points):1000],fontsize=16)
plt.yticks(fontsize=16)
#plt.title('Likelihood for Tau and Sigma Pairs',y=1.05,fontsize=17,fontweight='bold')
plt.axis([-500,10500,-5,15],fontsize=30,fontweight='bold')
plt.legend(loc=4)
plt.show()

contour=[]
for i in range (len(np.arange(np.amin(TSL_100['LIKELIHOOD']),150,60))):
	contour.append(np.arange(np.round(np.amin(TSL_100['LIKELIHOOD'])),150,60)[i])
for i in range (len(np.arange(150,np.amax(TSL_100['LIKELIHOOD']),7))):
	contour.append(np.arange(150,np.amax(TSL_100['LIKELIHOOD']),7)[i])
np.array(contour)[:] = np.array(contour)[::-1]

contour_con=[np.amax(q4565_g['LIKELIHOOD'])-(2.3/2),np.amax(q4565_g['LIKELIHOOD'])]
contour_max=[np.amax(q4565_g['LIKELIHOOD'])-0.00001,np.amax(q4565_g['LIKELIHOOD'])]

logtau=np.round(np.linspace(1,4,100),3)
tau=10**logtau
sig=np.round(np.linspace(0.001,0.021,100),4)
ntau=len(tau)
nsig=len(sig)

LH=np.reshape(q4565_g['LIKELIHOOD'],(ntau,nsig))
#print LH

plt.figure()
#cp = plt.contour(LH,10,levels=contour,colors='gray',linewidths=2)
cp = plt.contour(LH,10,colors='gray',linewidths=2)
cpc = plt.contour(LH,10,levels=contour_con,colors='k',linewidths=3)
cpm = plt.contour(LH,10,levels=contour_max,colors='k',linewidths=20)
plt.clabel(cp, inline=True,fontsize=13.5,colors='k',fontweight='bold')
plt.clabel(cpc, inline=False,fontsize=15,colors='white',fontweight='bold')
plt.imshow(LH, interpolation='nearest',cmap=plt.get_cmap('jet'),origin='lower')
plt.xticks(np.arange(0,nsig,10),sig[0:nsig:10],fontsize=17)
plt.yticks(np.arange(0,ntau,10),np.round(logtau,1)[0:ntau:10],fontsize=17)
plt.xlabel('$\sigma$'+' $(mag/day^{1/2})$',fontsize=27,fontweight='bold')
plt.ylabel('$log$ '+r'$\tau$'+' $(days)$',fontsize=27,fontweight='bold')
plt.title('RA334_DEC0.613__QSO4565_g',y=1.05,fontsize=17,fontweight='bold')
plt.gca().xaxis.set_label_coords(.5,-.09)
plt.gca().yaxis.set_label_coords(-.09,.5)
plt.colorbar()
plt.show()

plt.figure()
cp = plt.contour(LH,10,cmap=plt.get_cmap('jet'),linewidths=3)
plt.colorbar()
cpc = plt.contour(LH,10,levels=contour_con,colors='k',linewidths=3)
cpm = plt.contour(LH,10,levels=contour_max,colors='k',linewidths=20)
plt.clabel(cp, inline=True,fontsize=13.5,colors='k',fontweight='bold')
plt.clabel(cpc, inline=False,fontsize=13.5,colors='white',fontweight='bold')
plt.xticks(np.arange(0,nsig,10),sig[0:nsig:10],fontsize=17)
plt.yticks(np.arange(0,ntau,10),np.round(logtau,1)[0:ntau:10],fontsize=17)
plt.xlabel('$\sigma$'+' $(mag/day^{1/2})$',fontsize=27,fontweight='bold')
plt.ylabel('$log$ '+r'$\tau$'+' $(days)$',fontsize=27,fontweight='bold')
plt.title('RA334_DEC0.613__QSO4565_g',y=1.05,fontsize=17,fontweight='bold')
plt.gca().xaxis.set_label_coords(.5,-.09)
plt.gca().yaxis.set_label_coords(-.09,.5)
plt.show()
'''
def multi_heat(fileg, filer, filei, filez, logtau_start, logtau_stop, tau_div, sig_start, sig_stop, sig_div, xtick_step=1, ytick_step=1, confidence_font=7, title=' '):
	logtau=np.round(np.linspace(logtau_start,logtau_stop,tau_div),3)
	tau=10**logtau
	sig=np.round(np.linspace(sig_start,sig_stop,sig_div),4)
	ntau=len(tau)
	nsig=len(sig)
	tg=create_ascii_likelihood_table(fileg,1)
	tr=create_ascii_likelihood_table(filer,1)
	ti=create_ascii_likelihood_table(filei,1)
	tz=create_ascii_likelihood_table(filez,1)
	Lmaxg=np.amax(tg['LIKELIHOOD'])
	Lmaxr=np.amax(tr['LIKELIHOOD'])
	Lmaxi=np.amax(ti['LIKELIHOOD'])
	Lmaxz=np.amax(tz['LIKELIHOOD'])
	contour_maxg=[Lmaxg-0.001,Lmaxg]
	contour_maxr=[Lmaxr-0.001,Lmaxr]
	contour_maxi=[Lmaxi-0.001,Lmaxi]
	contour_maxz=[Lmaxz-0.001,Lmaxz]
	contour_cong=[Lmaxg-(2.3/2),Lmaxg]
	contour_conr=[Lmaxr-(2.3/2),Lmaxr]
	contour_coni=[Lmaxi-(2.3/2),Lmaxi]
	contour_conz=[Lmaxz-(2.3/2),Lmaxz]
	contour_g=[np.round(Lmaxg-5000), np.round(Lmaxg-2000), np.round(Lmaxg-1000), np.round(Lmaxg-500), np.round(Lmaxg-200), 			np.round(Lmaxg-100), np.round(Lmaxg-50), np.round(Lmaxg-20), np.round(Lmaxg-10), Lmaxg-(11.83/2), 			Lmaxg-(6.18/2)]
	contour_r=[np.round(Lmaxr-5000), np.round(Lmaxr-2000), np.round(Lmaxr-1000), np.round(Lmaxr-500), np.round(Lmaxr-200), 			np.round(Lmaxr-100), np.round(Lmaxr-50), np.round(Lmaxr-20), np.round(Lmaxr-10), Lmaxr-(11.83/2), 			Lmaxr-(6.18/2)]
	contour_i=[np.round(Lmaxi-5000), np.round(Lmaxi-2000), np.round(Lmaxi-1000), np.round(Lmaxi-500), np.round(Lmaxi-200), 			np.round(Lmaxi-100), np.round(Lmaxi-50), np.round(Lmaxi-20), np.round(Lmaxi-10), Lmaxi-(11.83/2), 			Lmaxi-(6.18/2)]
	contour_z=[np.round(Lmaxz-5000), np.round(Lmaxz-2000), np.round(Lmaxz-1000), np.round(Lmaxz-500), np.round(Lmaxz-200), 			np.round(Lmaxz-100), np.round(Lmaxz-50), np.round(Lmaxz-20), np.round(Lmaxz-10), Lmaxz-(11.83/2), 			Lmaxz-(6.18/2)]
	LHg=np.reshape(tg['LIKELIHOOD'],(ntau,nsig))
	LHr=np.reshape(tr['LIKELIHOOD'],(ntau,nsig))
	LHi=np.reshape(ti['LIKELIHOOD'],(ntau,nsig))
	LHz=np.reshape(tz['LIKELIHOOD'],(ntau,nsig))
	
	plt.subplot(2,2,1)
	cpcg = plt.contour(LHg,levels=contour_cong,colors='k',linewidths=3)
	cpmg = plt.contour(LHg,levels=contour_maxg,colors='white',linewidths=5)
	cpg = plt.contour(LHg,levels=contour_g,colors='grey',linewidths=2)
	plt.clabel(cpcg, inline=False,fontsize=confidence_font,colors='white',fontweight='bold')
	plt.clabel(cpg, inline=False,fontsize=13,colors='k',fontweight='bold')
	plt.imshow(LHg, interpolation='nearest',cmap=plt.get_cmap('jet'),origin='lower')
	plt.xticks(np.arange(0,nsig,xtick_step),'',fontsize=30)
	plt.yticks(np.arange(0,ntau,ytick_step),np.round(logtau,1)[0:ntau:ytick_step],fontsize=30)
	plt.ylabel('$log$ '+r'$\tau$'+' $(days)$',fontsize=50,fontweight='bold')
	plt.title(title,y=1.05,fontsize=17,fontweight='bold')
	plt.colorbar()	
	
	plt.subplot(2,2,2)
	cpcr = plt.contour(LHr,levels=contour_conr,colors='k',linewidths=3)
	cpmr = plt.contour(LHr,levels=contour_maxr,colors='white',linewidths=5)
	cpr =  plt.contour(LHr,levels=contour_r,colors='gray',linewidths=2)
	plt.clabel(cpcr, inline=False,fontsize=confidence_font,colors='white',fontweight='bold')
	plt.clabel(cpr, inline=False,fontsize=13,colors='k',fontweight='bold')
	plt.imshow(LHr, interpolation='nearest',cmap=plt.get_cmap('jet'),origin='lower')
	plt.xticks(np.arange(0,nsig,xtick_step),'',fontsize=30)
	plt.yticks(np.arange(0,ntau,ytick_step),'',fontsize=30)
	plt.colorbar()	

	plt.subplot(2,2,3)
	cpci = plt.contour(LHi,levels=contour_coni,colors='k',linewidths=3)
	cpmi = plt.contour(LHi,levels=contour_maxi,colors='white',linewidths=5)
	cpi =  plt.contour(LHi,levels=contour_i,colors='gray',linewidths=2)
	plt.clabel(cpci, inline=False,fontsize=confidence_font,colors='white',fontweight='bold')
	plt.clabel(cpi, inline=False,fontsize=13,colors='k',fontweight='bold')
	plt.imshow(LHi, interpolation='nearest',cmap=plt.get_cmap('jet'),origin='lower')
	plt.xticks(np.arange(0,nsig,xtick_step),sig[0:nsig:xtick_step],fontsize=30)
	plt.yticks(np.arange(0,ntau,ytick_step),np.round(logtau,1)[0:ntau:ytick_step],fontsize=30)
	plt.xlabel('$\sigma$'+' $(mag/day^{1/2})$',fontsize=50,fontweight='bold')
	plt.ylabel('$log$ '+r'$\tau$'+' $(days)$',fontsize=50,fontweight='bold')
	plt.gca().xaxis.set_label_coords(1.3,-.07)
	plt.colorbar()
	
	plt.subplot(2,2,4)
	cpcz = plt.contour(LHz,levels=contour_conz,colors='k',linewidths=3)
	cpmz = plt.contour(LHz,levels=contour_maxz,colors='white',linewidths=5)
	cpz =  plt.contour(LHz,levels=contour_z,colors='gray',linewidths=2)
	plt.clabel(cpcz, inline=False,fontsize=confidence_font,colors='white',fontweight='bold')
	plt.clabel(cpz, inline=False,fontsize=13,colors='k',fontweight='bold')
	plt.imshow(LHz, interpolation='nearest',cmap=plt.get_cmap('jet'),origin='lower')
	plt.xticks(np.arange(0,nsig,xtick_step),sig[0:nsig:xtick_step],fontsize=30)
	plt.yticks(np.arange(0,ntau,ytick_step),'',fontsize=30)
	plt.colorbar()
	plt.show()

'''
multi_heat(QSO2601_g25,QSO2601_r25,QSO2601_i25,QSO2601_z25,1,4,25,0.001,0.021,25,4,2)
multi_heat(QSO331_g25,QSO331_r25,QSO331_i25,QSO331_z25,1,4,25,0.001,0.021,25,4,2)
multi_heat(QSO715_g25,QSO715_r25,QSO715_i25,QSO715_z25,1,4,25,0.001,0.021,25,4,2)
multi_heat(QSO669_g25,QSO669_r25,QSO669_i25,QSO669_z25,1,4,25,0.001,0.021,25,4,2)
multi_heat(QSO1863_g25,QSO1863_r25,QSO1863_i25,QSO1863_z25,1,4,25,0.001,0.021,25,4,2)
multi_heat(QSO240_g25,QSO240_r25,QSO240_i25,QSO240_z25,1,4,25,0.001,0.021,25,4,2)
multi_heat(QSO2511_g25,QSO2511_r25,QSO2511_i25,QSO2511_z25,1,4,25,0.001,0.021,25,4,2)
multi_heat(QSO584_g25,QSO584_r25,QSO584_i25,QSO584_z25,1,4,25,0.001,0.021,25,4,2)
multi_heat(QSO557_g25,QSO557_r25,QSO557_i25,QSO557_z25,1,4,25,0.001,0.021,25,4,2)
multi_heat(QSO644_g25,QSO644_r25,QSO644_i25,QSO644_z25,1,4,25,0.001,0.021,25,4,2)
multi_heat(QSO616_g25,QSO616_r25,QSO616_i25,QSO616_z25,1,4,25,0.001,0.021,25,4,2)

multi_heat(QSO2601_g10,QSO2601_r10,QSO2601_i10,QSO2601_z10,1,4,10,0.001,0.021,10,2,1)
multi_heat(QSO331_g10,QSO331_r10,QSO331_i10,QSO331_z10,1,4,10,0.001,0.021,10,2,1)
multi_heat(QSO715_g10,QSO715_r10,QSO715_i10,QSO715_z10,1,4,10,0.001,0.021,10,2,1)
multi_heat(QSO669_g10,QSO669_r10,QSO669_i10,QSO669_z10,1,4,100,0.001,0.021,100,15,10)
multi_heat(QSO1863_g10,QSO1863_r10,QSO1863_i10,QSO1863_z10,1,4,10,0.001,0.021,10,2,1)
multi_heat(QSO240_g10,QSO240_r10,QSO240_i10,QSO240_z10,1,4,10,0.001,0.021,10,2,1)
multi_heat(QSO2511_g10,QSO2511_r10,QSO2511_i10,QSO2511_z10,1,4,10,0.001,0.021,10,2,1)
multi_heat(QSO584_g10,QSO584_r10,QSO584_i10,QSO584_z10,1,4,10,0.001,0.021,10,2,1)
multi_heat(QSO557_g10,QSO557_r10,QSO557_i10,QSO557_z10,1,4,10,0.001,0.021,10,2,1)
multi_heat(QSO644_g10,QSO644_r10,QSO644_i10,QSO644_z10,1,4,10,0.001,0.021,10,2,1)
multi_heat(QSO616_g10,QSO616_r10,QSO616_i10,QSO616_z10,1,4,10,0.001,0.021,10,2,1)'''

#multi_heat(QSO2601_g,QSO2601_r,QSO2601_i,QSO2601_z,1,4,100,0.001,0.021,100,15,10)
#multi_heat(QSO331_g,QSO331_r,QSO331_i,QSO331_z,1,4,100,0.001,0.021,100,15,10)
#multi_heat(QSO715_g,QSO715_r,QSO715_i,QSO715_z,1,4,100,0.001,0.021,100,15,10)
#multi_heat(QSO669_g,QSO669_r,QSO669_i,QSO669_z,1,4,100,0.001,0.021,100,15,10)
#multi_heat(QSO1863_g,QSO1863_r,QSO1863_i,QSO1863_z,1,4,100,0.001,0.021,100)
#multi_heat(QSO240_g,QSO240_r,QSO240_i,QSO240_z,1,4,100,0.001,0.021,100,15,10)
#multi_heat(QSO2511_g,QSO2511_r,QSO2511_i,QSO2511_z,1,4,100,0.001,0.021,100,15,10)
#multi_heat(QSO584_g,QSO584_r,QSO584_i,QSO584_z,1,4,100,0.001,0.021,100,15,10)
#multi_heat(QSO557_g,QSO557_r,QSO557_i,QSO557_z,1,4,100,0.001,0.021,100,15,10)
multi_heat(QSO644_g,QSO644_r,QSO644_i,QSO644_z,1,4,100,0.001,0.021,100,30,15)
#multi_heat(QSO616_g,QSO616_r,QSO616_i,QSO616_z,1,4,100,0.001,0.021,100,15,10)


#multi_heat(QSO3336_g,QSO3336_r,QSO3336_i,QSO3336_z,1,4,100,0.001,0.021,100,15,10)
#multi_heat(QSO5096_g,QSO5096_r,QSO5096_i,QSO5096_z,1,4,100,0.001,0.021,100,15,10)
#multi_heat(QSO4574_g,QSO4574_r,QSO4574_i,QSO4574_z,1,4,100,0.001,0.021,100,15,10)
#multi_heat(QSO6946_g,QSO6946_r,QSO6946_i,QSO6946_z,1,4,100,0.001,0.021,100,15,10)
#multi_heat(QSO3129_g,QSO3129_r,QSO3129_i,QSO3129_z,1,4,100,0.001,0.021,100,15,10)
#multi_heat(QSO6754_g,QSO6754_r,QSO6754_i,QSO6754_z,1,4,100,0.001,0.021,100,15,10)
#multi_heat(QSO3525_g,QSO3525_r,QSO3525_i,QSO3525_z,1,4,100,0.001,0.021,100,15,10)
#multi_heat(QSO6730_g,QSO6730_r,QSO6730_i,QSO6730_z,1,4,100,0.001,0.021,100,15,10)
#multi_heat(QSO1736_g,QSO1736_r,QSO1736_i,QSO1736_z,1,4,100,0.001,0.021,100,15,10)
#multi_heat(QSO4565_g,QSO4565_r,QSO4565_i,QSO4565_z,1,4,100,0.001,0.021,100,15,10)

def entire_likelihood_heatmap(fileg, filer, filei, filez,logtau_start,logtau_stop,tau_div,sig_start,sig_stop,sig_div,xtick_step=1,ytick_step=1,confidence_font=7):
	txtg = open(fileg, 'r')
	txtr = open(filer, 'r')
	txti = open(filei, 'r')
	txtz = open(filez, 'r')

	g=[];r=[];i=[];z=[]
	for rowg in txtg.readlines():
		g.append(rowg.split()[0])
	for rowr in txtr.readlines():
		r.append(rowr.split()[0])
	for rowi in txti.readlines():
		i.append(rowi.split()[0])
	for rowz in txtz.readlines(): 
		z.append(rowz.split()[0])
	txtg.close()
	txtr.close()
	txti.close()
	txtz.close()

	tb=Table([np.array(g),np.array(r),np.array(i),np.array(z)],
		names=('G','R','I','Z'))

	for t in range(len(tb)):
		#print tb[t]
		multi_heat(tb['G'][t]+'_tau'+str(tau_div)+'_sig'+str(sig_div)+'_table.dat', tb['R'][t]+'_tau'+str(tau_div)+'_sig'+str(sig_div)+'_table.dat', tb['I'][t]+'_tau'+str(tau_div)+'_sig'+str(sig_div)+'_table.dat', tb['Z'][t]+'_tau'+str(tau_div)+'_sig'+str(sig_div)+'_table.dat', logtau_start, logtau_stop, tau_div, sig_start, sig_stop, sig_div, xtick_step, ytick_step,confidence_font, str(tb['G'][t]))
'''
#entire_likelihood_heatmap(agRA33to35,arRA33to35,aiRA33to35,azRA33to35,0.5,4.5,10,0.001,0.041,10)
entire_likelihood_heatmap(agRA34to45,arRA34to45,aiRA34to45,azRA34to45,0.5,4.5,10,0.001,0.041,10)
entire_likelihood_heatmap(agRA35to55,arRA35to55,aiRA35to55,azRA35to55,0.5,4.5,10,0.001,0.041,10)
entire_likelihood_heatmap(agRA25to33,arRA25to33,aiRA25to33,azRA25to33,0.5,4.5,10,0.001,0.041,10)
entire_likelihood_heatmap(agRA35to34,arRA35to34,aiRA35to34,azRA35to34,0.5,4.5,10,0.001,0.041,10)
entire_likelihood_heatmap(agRA45to35,arRA45to35,aiRA45to35,azRA45to35,0.5,4.5,10,0.001,0.041,10)'''

LH_dict=create_entire_likelihood_stats_dic()
print len(LH_dict['g'])
def multi_hist_logtau_sigma(data1, data2, data3, data4,ranges,label_type,axis_type,axis,A=8,loc=2):
	f, ax = plt.subplots(2, 2,sharex=True,sharey=True)
	ax[0,0].hist(data1, range=ranges, bins=A, histtype='stepfilled',color='purple', label='g band best-fit '+label_type)
	ax[0,0].axis(axis,fontsize=30,fontweight='bold')
	ax[0,0].legend(loc=loc)
	
	ax[0,1].hist(data2, range=ranges, bins=A, histtype='stepfilled',color='blue', label='r band best-fit '+label_type)
	ax[0,1].axis(axis,fontsize=30,fontweight='bold')
	ax[0,1].legend(loc=loc)
	
	ax[1,0].hist(data3, range=ranges, bins=A, histtype='stepfilled',color='green', label='i band best-fit '+label_type)
	ax[1,0].axis(axis,fontsize=30,fontweight='bold')
	ax[1,0].legend(loc=loc)
	
	ax[1,1].hist(data4, range=ranges, bins=A, histtype='stepfilled',color='orange', label='z band best-fit '+label_type)
	ax[1,0].axis(axis,fontsize=30,fontweight='bold')
	ax[1,1].legend(loc=loc)

	# Fine-tune figure; hide x ticks for top plots and y ticks for right plots
	plt.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
	plt.setp([a.get_yticklabels() for a in ax[:, 1]], visible=False)
	ax[1,0].set_xlabel('$Maximum   Likelihood $ '+axis_type,fontsize=30, x=1.05)
	ax[0,0].set_ylabel('$Number$',fontsize=30)
	ax[1,0].set_ylabel('$Number$',fontsize=30)
	plt.show()

#multi_hist_logtau_sigma(LH_dict['g']['LOG_TAU_LH'],LH_dict['r']['LOG_TAU_LH'],LH_dict['i']['LOG_TAU_LH'],LH_dict['z']['LOG_TAU_LH'],(1, 4),'log tau','$log$ '+r'$\tau$'+' $(days)$',[0.5, 4.5,0,37],12)
#multi_hist_logtau_sigma(LH_dict['g']['SIGMA_LH'],LH_dict['r']['SIGMA_LH'],LH_dict['i']['SIGMA_LH'],LH_dict['z']['SIGMA_LH'],(0.001, 0.021),'sigma','$\sigma$'+' $(mag/day^{1/2})$',[0,0.022,0,30],10)

#vary_mean_std=create_entiresample_mean_std_lenght_table(RA33to35,vRA33to35,RA34to45,vRA34to45,RA35to55,vRA35to55,RA25to33,vRA25to33,RA35to34,vRA35to34,RA45to35,vRA45to35)

# creates four subplots sharing both x/y axes of mean magnitudes of the light curves for the four bands
def multi_hist_mean(data1, data2, data3, data4, A=9,loc=2):
	f, ax = plt.subplots(2, 2,sharex=True,sharey=True)
	ax[0,0].hist(data1, bins=A,range=(16.5,22.5),histtype='bar',color='purple', label='g mag mean')
	ax[0,0].axis([16,23,0,35],fontsize=30,fontweight='bold')
	ax[0,0].legend(loc=loc)
	
	ax[0,1].hist(data2, bins=A,range=(16.5,22.5),histtype='bar',color='blue', label='r mag mean')
	ax[0,1].axis([16,23,0,35],fontsize=30,fontweight='bold')
	ax[0,1].legend(loc=loc)
	
	ax[1,0].hist(data3, bins=A,range=(16.5,22.5),histtype='bar',color='green', label='i mag mean')
	ax[1,0].axis([16,23,0,35],fontsize=30,fontweight='bold')
	ax[1,0].legend(loc=loc)
	
	ax[1,1].hist(data4, bins=A,range=(16.5,22.5),histtype='bar',color='orange', label='z mag mean')
	ax[1,1].axis([16,23,0,35],fontsize=30,fontweight='bold')
	ax[1,1].legend(loc=loc)

	# Fine-tune figure; hide x ticks for top plots and y ticks for right plots
	plt.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
	plt.setp([a.get_yticklabels() for a in ax[:, 1]], visible=False)
	ax[1,0].set_xlabel('Magnitude Means',fontsize=16, x=1.05)
	ax[0,0].set_ylabel('Number',fontsize=16)
	ax[1,0].set_ylabel('Number',fontsize=16)
	plt.show()

#multi_hist_mean(vary_mean_std['MEAN_g'],vary_mean_std['MEAN_r'],vary_mean_std['MEAN_i'],vary_mean_std['MEAN_z'])

#creates four subplots sharing both x/y axes of magnitude standard deviation of the light curves for the four bands
def multi_hist_std(data1, data2, data3, data4, A=6,loc=2):
	f, ax = plt.subplots(2, 2,sharex=True,sharey=True)
	ax[0,0].hist(data1, bins=A,range=(0,0.3),histtype='bar',color='purple', label='g mag standard deviation')
	ax[0,0].axvline(data1.mean(), color='k', linestyle='--', linewidth=2, label='g mag mean standard deviation: %.3f'% np.mean(data1))
	ax[0,0].axis([-0.01,0.31,0,45],fontsize=30,fontweight='bold')
	ax[0,0].legend(loc=loc)
	
	ax[0,1].hist(data2, bins=A,range=(0,0.3), histtype='bar',color='blue', label='r mag standard deviation')
	ax[0,1].axvline(data2.mean(), color='k', linestyle='--', linewidth=2, label='r mag mean standard deviation: %.3f'% np.mean(data2))
	ax[0,1].axis([-0.01,0.31,0,45],fontsize=30,fontweight='bold')
	ax[0,1].legend(loc=loc)
	
	ax[1,0].hist(data3, bins=A,range=(0,0.3), histtype='bar',color='green', label='i mag standard deviation')
	ax[1,0].axvline(data3.mean(), color='k', linestyle='--', linewidth=2, label='i mag mean standard deviation: %.3f'% np.mean(data3))
	ax[1,0].axis([-0.01,0.31,0,45],fontsize=30,fontweight='bold')
	ax[1,0].legend(loc=loc)
	
	ax[1,1].hist(data4, bins=A,range=(0,0.3), histtype='bar',color='orange', label='z mag standard deviation')
	ax[1,1].axvline(data4.mean(), color='k', linestyle='--', linewidth=2, label='z mag mean standard deviation: %.3f'% np.mean(data4))
	ax[1,0].axis([-0.01,0.31,0,45],fontsize=30,fontweight='bold')
	ax[1,1].legend(loc=loc)

	# Fine-tune figure; hide x ticks for top plots and y ticks for right plots
	plt.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
	plt.setp([a.get_yticklabels() for a in ax[:, 1]], visible=False)
	ax[1,0].set_xlabel("Magnitude Standard Deviations",fontsize=16, x=1.05)
	ax[0,0].set_ylabel('Number',fontsize=16)
	ax[1,0].set_ylabel('Number',fontsize=16)
	plt.show()

#multi_hist_std(vary_mean_std['STD_g'],vary_mean_std['STD_r'],vary_mean_std['STD_i'],vary_mean_std['STD_z'])


def multi_tau_sigma(data1, data2, data3, data4, A=9,loc=2):
	f, ax = plt.subplots(2, 2,sharex=True,sharey=True)
	ax[0,0].hist(data1, bins=A,range=(16.5,22.5),histtype='bar',color='purple', label='g mag mean')
	ax[0,0].axis([16,23,0,35],fontsize=30,fontweight='bold')
	ax[0,0].legend(loc=loc)
	
	ax[0,1].hist(data2, bins=A,range=(16.5,22.5),histtype='bar',color='blue', label='r mag mean')
	ax[0,1].axis([16,23,0,35],fontsize=30,fontweight='bold')
	ax[0,1].legend(loc=loc)
	
	ax[1,0].hist(data3, bins=A,range=(16.5,22.5),histtype='bar',color='green', label='i mag mean')
	ax[1,0].axis([16,23,0,35],fontsize=30,fontweight='bold')
	ax[1,0].legend(loc=loc)
	
	ax[1,1].hist(data4, bins=A,range=(16.5,22.5),histtype='bar',color='orange', label='z mag mean')
	ax[1,1].axis([16,23,0,35],fontsize=30,fontweight='bold')
	ax[1,1].legend(loc=loc)

	# Fine-tune figure; hide x ticks for top plots and y ticks for right plots
	plt.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
	plt.setp([a.get_yticklabels() for a in ax[:, 1]], visible=False)
	ax[1,0].set_xlabel('Magnitude Means',fontsize=16, x=1.05)
	ax[0,0].set_ylabel('Number',fontsize=16)
	ax[1,0].set_ylabel('Number',fontsize=16)
	plt.show()

#multi_hist_mean(vary_mean_std['MEAN_g'],vary_mean_std['MEAN_r'],vary_mean_std['MEAN_i'],vary_mean_std['MEAN_z'])

def tau_sig_error_plot(gLH_sigma,gsig_cmin,gsig_cmax,gLH_tau,gtau_cmin,gtau_cmax,rLH_sigma,rsig_cmin,rsig_cmax,rLH_tau,rtau_cmin,rtau_cmax,
iLH_sigma,isig_cmin,isig_cmax,iLH_tau,itau_cmin,itau_cmax, zLH_sigma,zsig_cmin,zsig_cmax,zLH_tau,ztau_cmin,ztau_cmax,loc=1):
	plt.subplot(2,2,1)
	plt.errorbar(gLH_sigma, gLH_tau, xerr = [gLH_sigma - gsig_cmin, gsig_cmax - gLH_sigma], yerr=[gLH_tau - gtau_cmin, gtau_cmax - gLH_tau], marker='o',ms=10, c='purple', mec='purple', mfc='purple', ecolor='purple', label='g band', linestyle = 'None')
	plt.axis([0,0.025,0.5,4.5],fontsize=30,fontweight='bold')
	plt.ylabel('$Maximum$'+' '+'$Likelihood$'+' '+'$ log$ '+r'$\tau$'+' $(days)$',fontsize=27,fontweight='bold')
	plt.xticks((0,0.025),'',fontsize=15)
	plt.gca().yaxis.set_label_coords(-0.07,0.5)
	plt.legend(loc=loc)
	
	plt.subplot(2,2,2)
	plt.errorbar(rLH_sigma, rLH_tau, xerr = [rLH_sigma - rsig_cmin, rsig_cmax - rLH_sigma], yerr=[rLH_tau - rtau_cmin, rtau_cmax - rLH_tau], marker='o',ms=10, c='b', mec='b', mfc='b', ecolor='b', label='r band', linestyle = 'None')
	plt.xticks((0,0.025),'',fontsize=15)
	plt.yticks((0.5,0.45),'',fontsize=15)
	plt.axis([0.,0.025,0.5,4.5],fontsize=30,fontweight='bold')
	plt.legend(loc=loc)

	plt.subplot(2,2,3)
	plt.errorbar(iLH_sigma, iLH_tau, xerr = [iLH_sigma - isig_cmin, isig_cmax - iLH_sigma], yerr=[iLH_tau - itau_cmin, itau_cmax - iLH_tau], marker='o',ms=10, c='g', mec='g', mfc='g', ecolor='g', label='i band', linestyle = 'None')
	plt.axis([0.,0.025,0.5,4.5],fontsize=30,fontweight='bold')
	plt.xlabel('$Maximum$'+' '+'$Likelihood$'+' '+'$\sigma$'+' $(mag/day^{1/2})$',fontsize=27,fontweight='bold')
	plt.ylabel('$Maximum$'+' '+'$Likelihood$'+' '+'$log$ '+r'$\tau$'+' $(days)$',fontsize=27,fontweight='bold')
	plt.gca().xaxis.set_label_coords(1,-.07)
	plt.gca().yaxis.set_label_coords(-0.07,0.5)
	plt.legend(loc=loc)
	
	plt.subplot(2,2,4)
	plt.errorbar(zLH_sigma, zLH_tau, xerr = [zLH_sigma - zsig_cmin, zsig_cmax - zLH_sigma], yerr=[zLH_tau - ztau_cmin, ztau_cmax - zLH_tau], marker='o',ms=10, c='orange', mec='orange', mfc='orange', ecolor='orange', label='z band ', linestyle = 'None')
	plt.yticks((0.5,0.45),'',fontsize=15)
	plt.axis([0.,0.025,0.5,4.5],fontsize=30,fontweight='bold')
	plt.legend(loc=loc)
	plt.show()

#tau_sig_error_plot(LH_dict['g']['SIGMA_LH'],LH_dict['g']['SIGMA_Cmin'],LH_dict['g']['SIGMA_Cmax'], LH_dict['g']['LOG_TAU_LH'], LH_dict['g']['LOG_TAU_Cmin'],LH_dict['g']['LOG_TAU_Cmax'], LH_dict['r']['SIGMA_LH'],LH_dict['r']['SIGMA_Cmin'],LH_dict['r']['SIGMA_Cmax'], LH_dict['r']['LOG_TAU_LH'], LH_dict['r']['LOG_TAU_Cmin'],LH_dict['r']['LOG_TAU_Cmax'],LH_dict['i']['SIGMA_LH'],LH_dict['i']['SIGMA_Cmin'],LH_dict['i']['SIGMA_Cmax'], LH_dict['i']['LOG_TAU_LH'], LH_dict['i']['LOG_TAU_Cmin'],LH_dict['i']['LOG_TAU_Cmax'],LH_dict['z']['SIGMA_LH'],LH_dict['z']['SIGMA_Cmin'],LH_dict['z']['SIGMA_Cmax'],LH_dict['z']['LOG_TAU_LH'], LH_dict['z']['LOG_TAU_Cmin'],LH_dict['z']['LOG_TAU_Cmax'],)


'''
txt = open("variables_RA333to3335.txt", "w")
txt.write("616\n644\n1025\n1285\n2084\n2233\n2870\n2922\n3377\n3680\n3856\n4101\n4134\n4874\n4900\n5040\n5058\n5439\n5629\n6628\n7060")
txt.close()

txt = open("variables_RA333to3335_test.txt", "w")
txt.write("616\n644")
txt.close()

txt = open("variables_RA334to3345.txt", "w")
txt.write("557\n584\n593\n787\n1397\n1749\n2134\n4085\n4184\n4225\n4736\n5446\n5676\n6674\n6721\n7098\n7207\n7398\n7462\n7669\n7828\n8058\n8235\n8239")
txt.close()

txt = open("variables_RA334to3345_test.txt", "w")
txt.write("557\n584")
txt.close()

txt = open("variables_RA335to3355.txt", "w")
txt.write("2511")
txt.close()

txt = open("variables_RA335to3355_test.txt", "w")
txt.write("2511")
txt.close()

txt = open("variables_RA3325to333.txt", "w")
txt.write("240\n1863\n1932\n2658\n2887\n3392")
txt.close()

txt = open("variables_RA3325to333_test.txt", "w")
txt.write("240\n1863")
txt.close()

txt = open("variables_RA3335to334.txt", "w")
txt.write("669\n715\n747\n1137\n1525\n1840\n2094\n2363\n2375\n2462\n3002\n3069\n3385\n3389\n3902\n4216\n4342\n4516\n4569\n4990\n4995\n5091\n5111\n5589\n5782\n5795\n6273\n6750\n6817\n7059\n7686\n8203\n8258")
txt.close()

txt = open("variables_RA3335to334_test.txt", "w")
txt.write("669\n715")
txt.close()

txt = open("variables_RA3345to335.txt", "w")
txt.write("331\n2601\n2861\n2911\n3031\n3178\n3343\n4052\n4115\n4378\n4464\n4550\n4750\n5004\n5018\n5300\n6226\n7363\n7460")
txt.close()

txt = open("variables_RA3345to335_test.txt", "w")
txt.write("331\n2601")
txt.close()'''

