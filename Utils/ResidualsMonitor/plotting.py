#!/usr/bin/env python
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np

# Set case name
case_name = "Default"

# Get residuals
residuals_file = ".residuals"

def plotResidualsFile(fileName):
	residuals_fid = open( fileName, 'r' )
	plt1 = plt.figure(1)
	counter = 0
	iter = []
	time = []
	continuity = []
	x_momentum = []
	y_momentum = []
	energy = [] 
	for line in residuals_fid:
		if counter > 1:
			numbers = line.split()
			iter.append(float(numbers[0]))
			time.append(float(numbers[1]))
			continuity.append(float(numbers[2]))
			x_momentum.append(float(numbers[3]))
			y_momentum.append(float(numbers[4]))
			energy.append(float(numbers[5]))
		counter = counter + 1

	ax = plt.subplot(2,1,1)
	ax.clear()
	cont = ax.semilogy(time,continuity,'-r',label='continuity')
	xmom = ax.semilogy(time,x_momentum,'-b',label='x_momentum')
	ymom = ax.semilogy(time,y_momentum,'-k',label='y_momentum')
	ener = ax.semilogy(time,energy,'-m',label='energy')
	plt.ylabel('residuals')
	plt.grid()
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles,labels)
	#plt.show()
	plt1.canvas.manager.window.attributes('-topmost',0)
	pdfname = fileName + '.pdf'
	#plt1.savefig(pdfname)

def plotProbeFile(fileName):
	probe_fid = open( fileName, 'r' )
	plt2=plt.figure(2)
	counter = 0
	iter = []
	time = []
	value = []
	for line in probe_fid:
		if counter > 6:
			numbers = line.split()
			iter.append(float(numbers[0]))
			time.append(float(numbers[1]))
			value.append(float(numbers[2]))
		counter = counter + 1

	ax = plt.subplot(1,1,1)
	ax.plot(time,value,'-k')
	plt.xlabel('time [-]')
	plt.ylabel('variable')
	#plt.show()
	pdfname = fileName + '.pdf'
	plt2.savefig(pdfname)

def plotLiftDragFile(fileName):
	lift_fid = open( fileName, 'r' )
	plt2=plt.figure(1)
	counter = 0
	iter = []
	time = []
	value = []
	for line in lift_fid:
		if counter > 7:
			numbers = line.split()
			iter.append(float(numbers[0]))
			time.append(float(numbers[1]))
			value.append(float(numbers[2]))
		counter = counter + 1

	ax = plt.subplot(2,1,2)
	ax.clear()
	ax.plot(time,value,'-k')
	plt.xlabel('time [-]')
	plt.ylabel('$c_f $ wall')
	#plt.show()
	pdfname = fileName + '.pdf'
	#plt2.savefig(pdfname)
	plt.grid()
	plt2.canvas.manager.window.attributes('-topmost',0)


plt.ion()
while True:
	plotResidualsFile(residuals_file)
	plt.pause(5)
		
