#! /usr/bin/python

import os

files=os.listdir(os.getcwd())

files.sort()

calc_mode='NULL'
trans_mode='NULL'

for d in files:

	fname="%s/output.doz"%d

	if os.path.isfile(fname):

		f=open(fname,'r')

		text=f.read()

		tList=text.split()

		i=1
		while i < len(tList):

			if tList[i]=="trans_mode":
				trans_mode=tList[i+2]
			if tList[i]=="calc_mode":
				calc_mode=tList[i+2]
			if (calc_mode!='NULL') and (trans_mode!='NULL'):
				break
			i+=1

		var_param='NULL'

		var_liq='NULL'
                var_glas='NULL'
		
                lamb_liq='NULL'
		lamb_glas='NULL'

		eigen_liq='NULL'
		eigen_glas='NULL'

		eigen_last='NULL'

		i=1

		if (trans_mode=='cont') or (trans_mode=='loca'):

			while i < len(tList):
				if tList[i]=="eigenvalue":
					if float(tList[i+2])>1.0:
                        			eigen_glas=tList[i+2]
						eigen_last='glas'
                        		elif float(tList[i+2])<1.0:
						eigen_liq=tList[i+2]
						eigen_last='liq'
				if tList[i]=="lambda":
                        		if eigen_last=='glas':
						lamb_glas=tList[i+2]
                        		elif eigen_last=='liq':
						lamb_liq=tList[i+2]
                        		eigen_last='NULL'
				if tList[i]=="var_param":
					var_param=tList[i+2]
					var_liq="%s_liq"%var_param
					var_glas="%s_glas"%var_param
				if tList[i]==var_liq:
                        		liq=tList[i+2]
				if tList[i]==var_glas:
                        		glas=tList[i+2]
				if tList[i]=="density":
                        		density=tList[i+2]
				if tList[i]=="delta":
					delta=tList[i+2]

				i+=1

			if var_param=='dens':
				print liq,glas,delta,delta,lamb_liq,lamb_glas,eigen_liq,eigen_glas
			elif var_param=='delt':
				print density,density,liq,glas,lamb_liq,lamb_glas,eigen_liq,eigen_glas
			else:
				print density,density,delta,delta,lamb_liq,lamb_glas,eigen_liq,eigen_glas

		elif trans_mode=='disc':

			while i < len(tList):
				if tList[i]=="eigenvalue":
                        		eigen_glas=tList[i+2]
                        	if tList[i]=="eigen_infl":
					eigen_liq=tList[i+2]
					eigen_last='liq'
				if tList[i]=="lambda":
					lamb_glas=tList[i+2]
				if tList[i]=="lamb_infl":
					lamb_liq=tList[i+2]
				if tList[i]=="var_param":
					var_param=tList[i+2]
					var_liq="%s_liq"%var_param
					var_glas="%s_glas"%var_param
				if tList[i]==var_liq:
                        		liq=tList[i+2]
				if tList[i]==var_glas:
                        		glas=tList[i+2]
				if tList[i]=="density":
                        		density=tList[i+2]
				if tList[i]=="delta":
					delta=tList[i+2]

				i+=1

			if var_param=='dens':
				print liq,glas,delta,delta,lamb_liq,lamb_glas,eigen_liq,eigen_glas
			elif var_param=='delt':
				print density,density,liq,glas,lamb_liq,lamb_glas,eigen_liq,eigen_glas
			else:
				print density,density,delta,delta,lamb_liq,lamb_glas,eigen_liq,eigen_glas


