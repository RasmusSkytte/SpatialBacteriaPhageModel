#!/usr/bin/env python3

import os
from subprocess import call
from math import floor
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LogFormatterSciNotation
import shutil

def readLog( ) :
	with open('log.txt','r') as f_log :

		return_dict = {}

		# Load first line
		for line in f_log :
			line = line.split('=')
			try :
				return_dict[line[0].strip()] = float(line[1].strip())
			except:
				pass


	return return_dict

def reshape( array ) :

	# Size of space
	N = np.size(array,1)

	# Size of time
	K = int((round(np.size(array)/pow(N,3))))

	# Allocate array
	M = np.full((N,N,N,K),np.nan)

	# Do the reshaping
	for k in range(K) :
		for z in range(N) :
			for y in range(N) :
				for x in range(N) :
					M[y][x][z][k] = array[y+N*z+N*N*k][x]
	return M


def Visualize3D() :

	cwd = os.getcwd()
	path = cwd.replace('analysis','cpp/data/3D_Example_Full_Model/')
	os.chdir(path)
	exit = False

	# Check folder has the correct files
	if not os.path.isfile('CellDensity.txt') :
		exit = True
	if not os.path.isfile('PhageDensity.txt') :
		exit = True
	if not os.path.isfile('InfectedDensity.txt') :
		exit = True

	# Check all requirements are met
	if exit:
		os.chdir(cwd)
		return

	# Load the meta data
	params = readLog()

	# Recast into int
	if 'nGrid' in params:
		params['nGridXY'] = int(params['nGrid'])
		params['nGridZ']  = int(params['nGrid'])
		params['H'] = params['L']
	else:
		params['nGridXY'] = int(params['nGridXY'])
		params['nGridZ'] = int(params['nGridZ'])

	params['nSamp'] = int(params['nSamp'])
	params['T_end'] = int(params['T_end'])

	# Check model is spatial
	if params['nGridXY'] == 1:
		os.chdir(cwd)
		return

	# Check if model has been visualzed
	T = np.linspace( 0, params['T_end'], 1+params['T_end']*params['nSamp'] )
	exit = True
	for i in range(T.size) :

		outputPath = path + 'Visualize3D_py/%03i.png' % (i+1)
		if not os.path.isfile(outputPath) :
			exit = False

	if exit :
		return

	Bf = open('CellDensity.txt','r')
	Pf = open('PhageDensity.txt','r')
	If = open('InfectedDensity.txt','r')

	# Create output folders
	if not os.path.isdir('Visualize3D_py') :
		os.mkdir('Visualize3D_py')

	if not os.path.isdir('Visualize3D_py/Cells') :
		os.mkdir('Visualize3D_py/Cells')

	if not os.path.isdir('Visualize3D_py/Cells_density') :
		os.mkdir('Visualize3D_py/Cells_density')

	if not os.path.isdir('Visualize3D_py/Phages') :
		os.mkdir('Visualize3D_py/Phages')

	if not os.path.isdir('Visualize3D_py/Infected') :
		os.mkdir('Visualize3D_py/Infected')

	# Create figure
	fig = plt.figure(figsize=(3, 2.5), dpi = 300)
	ax = fig.add_subplot(111, projection='3d')
	ax.view_init(30, -42.5)
	fig.subplots_adjust(left=0.06, right=0.8, bottom=0.15, top=1)
	font = {'fontname':'arial'}

	xyticks = [-params['L']/2, -params['L']/4, 0, params['L']/4, params['L']/2]
	zticks = [-params['H']/2, -params['H']/4, 0, params['H']/4, params['H']/2]
	xlabels = ['$%.1f\cdot 10^{%.0f}$      ' % (t/pow(10,np.floor(np.log10(abs(t)))), np.floor(np.log10(abs(t)))) for t in xyticks]
	xlabels[(len(xyticks)-1)//2] = '$0.0$'
	ylabels = ['         $%.1f\cdot 10^{%.0f}$' % (t/pow(10,np.floor(np.log10(abs(t)))), np.floor(np.log10(abs(t)))) for t in xyticks]
	ylabels[(len(xyticks)-1)//2] = '$0.0$'
	zlabels = ['       $%.1f\cdot 10^{%.0f}$' % (t/pow(10,np.floor(np.log10(abs(t)))), np.floor(np.log10(abs(t)))) for t in zticks]
	zlabels[(len(zticks)-1)//2] = '$0.0$'

	# Prepare x,y,z coordinates for scatter plot
	xy = np.linspace(-params['L']/2,params['L']/2,params['nGridXY'])
	z = np.linspace(-params['H']/2,params['H']/2,params['nGridZ'])
	xs, ys, zs = np.meshgrid(xy, xy, z, sparse=False, indexing='ij')
	xs = -xs.reshape((pow(params['nGridXY'],2)*params['nGridZ'],))
	ys = -ys.reshape((pow(params['nGridXY'],2)*params['nGridZ'],))
	zs = zs.reshape((pow(params['nGridXY'],2)*params['nGridZ'],))

	# Remove header
	Bf.readline()
	Pf.readline()
	If.readline()

	i = 0

	# Read until EOF
	for l in range(params['nSamp']*params['T_end']+1):

		# Allocate
		B = np.zeros((params['nGridXY'],params['nGridXY'],params['nGridZ']))
		P = np.zeros((params['nGridXY'],params['nGridXY'],params['nGridZ']))
		I = np.zeros((params['nGridXY'],params['nGridXY'],params['nGridZ']))
		C = np.zeros((params['nGridXY'],params['nGridXY'],params['nGridZ']))

		# Read lines corresponding to this time point
		exit = False
		for k in range(params['nGridZ']):
			if exit : break
			for j in range(params['nGridXY']):
				# Writing ugly code because python is shit
				line = Bf.readline()
				if not line:
					exit = True
					break
				B[:,j,k] = np.fromstring(line,sep='\t')
				P[:,j,k] = np.fromstring(Pf.readline(),sep='\t')
				I[:,j,k] = np.fromstring(If.readline(),sep='\t')
				C[:,j,k] = B[:,j,k] + I[:,j,k]

		# Plot data if it is not already
		outputPath = path + 'Visualize3D_py/Phages/%03i.png' % i

		if not os.path.isfile(outputPath) :

			# Get the populations to plot
			Bs  = B.reshape((pow(params['nGridXY'],2)*params['nGridZ'],))
			Is  = I.reshape((pow(params['nGridXY'],2)*params['nGridZ'],))
			Cs  = C.reshape((pow(params['nGridXY'],2)*params['nGridZ'],))
			Ps  = P.reshape((pow(params['nGridXY'],2)*params['nGridZ'],))


			ind = np.where(Bs < 1)
			xsB = np.delete(xs,ind)
			ysB = np.delete(ys,ind)
			zsB = np.delete(zs,ind)
			Bs  = np.delete(Bs,ind)

			ind = np.where(Is < 1)
			xsI = np.delete(xs,ind)
			ysI = np.delete(ys,ind)
			zsI = np.delete(zs,ind)
			Is  = np.delete(Is,ind)

			ind = np.where(Cs < 1)
			xsC = np.delete(xs,ind)
			ysC = np.delete(ys,ind)
			zsC = np.delete(zs,ind)
			Cs  = np.delete(Cs,ind)

			ind = np.where(Ps < 1)
			xsP = np.delete(xs,ind)
			ysP = np.delete(ys,ind)
			zsP = np.delete(zs,ind)
			Ps  = np.delete(Ps,ind)

			ax.scatter(xsB, ysB, zsB, s=pow(np.log10(Bs),2), c='b', alpha=0.6, edgecolors='None', linewidths=0)
			ax.set_xlim(-params['L']/2,params['L']/2)
			ax.set_ylim(-params['L']/2,params['L']/2)
			ax.set_zlim(-params['H']/2,params['H']/2)

			outputPath = path + 'Visualize3D_py/Cells/%03i.png' % i

			# plt.tight_layout()
			ax.set_xticks(xyticks)
			ax.set_yticks(xyticks)
			ax.set_zticks(zticks)
			ax.set_xticklabels(xlabels,fontsize=8, verticalalignment='top', **font)
			ax.set_yticklabels(ylabels,fontsize=8, verticalalignment='top', **font)
			ax.set_zticklabels(zlabels,fontsize=8, verticalalignment='top', **font)
			ax.set_xlabel('x (µm)',fontsize=12,labelpad=10, **font)
			ax.set_ylabel('y (µm)',fontsize=12,labelpad=8, **font)
			ax.set_zlabel('z (µm)',fontsize=12,labelpad=12, **font)
			plt.savefig(outputPath, dpi=300)
			plt.cla()

			Cs = np.round(np.log10(Cs)) + 1
			sortKey = np.argsort(Cs)

			Cs = Cs[sortKey]
			xsC = xsC[sortKey]
			ysC = ysC[sortKey]
			zsC = zsC[sortKey]

			sizes = np.unique(Cs)
			for size in sizes:
				ind = np.where(Cs == size)
				ax.scatter(xsC[ind], ysC[ind], zsC[ind], s=1, c='b', alpha=min(1.0,size/6), edgecolors='None', linewidths=0)

			ax.set_xlim(-params['L']/2,params['L']/2)
			ax.set_ylim(-params['L']/2,params['L']/2)
			ax.set_zlim(-params['H']/2,params['H']/2)

			outputPath = path + 'Visualize3D_py/Cells_density/%03i.png' % i

			# plt.tight_layout()
			ax.set_xticks(xyticks)
			ax.set_yticks(xyticks)
			ax.set_zticks(zticks)
			ax.set_xticklabels(xlabels,fontsize=8, verticalalignment='top', **font)
			ax.set_yticklabels(ylabels,fontsize=8, verticalalignment='top', **font)
			ax.set_zticklabels(zlabels,fontsize=8, verticalalignment='top', **font)
			ax.set_xlabel('x (µm)',fontsize=12,labelpad=10, **font)
			ax.set_ylabel('y (µm)',fontsize=12,labelpad=8, **font)
			ax.set_zlabel('z (µm)',fontsize=12,labelpad=12, **font)
			plt.savefig(outputPath, dpi=300)
			plt.cla()


			ax.scatter(xsI, ysI, zsI, s=pow(np.log10(Is),2), c='g', alpha=0.6, edgecolors='None', linewidths=0)
			ax.set_xlim(-params['L']/2,params['L']/2)
			ax.set_ylim(-params['L']/2,params['L']/2)
			ax.set_zlim(-params['H']/2,params['H']/2)

			outputPath = path + 'Visualize3D_py/Infected/%03i.png' % i

			# plt.tight_layout()
			ax.set_xticks(xyticks)
			ax.set_yticks(xyticks)
			ax.set_zticks(zticks)
			ax.set_xticklabels(xlabels,fontsize=8, verticalalignment='top', **font)
			ax.set_yticklabels(ylabels,fontsize=8, verticalalignment='top', **font)
			ax.set_zticklabels(zlabels,fontsize=8, verticalalignment='top', **font)
			ax.set_xlabel('x (µm)',fontsize=12,labelpad=10, **font)
			ax.set_ylabel('y (µm)',fontsize=12,labelpad=8, **font)
			ax.set_zlabel('z (µm)',fontsize=12,labelpad=12, **font)
			plt.savefig(outputPath, dpi=300)
			plt.cla()


			Ps = np.round(np.log10(Ps)) + 1
			sortKey = np.argsort(Ps)

			Ps = Ps[sortKey]
			xsP = xsP[sortKey]
			ysP = ysP[sortKey]
			zsP = zsP[sortKey]

			sizes = np.unique(Ps)
			for size in sizes:
				ind = np.where(Ps == size)
				ax.scatter(xsP[ind], ysP[ind], zsP[ind], s=1, c='r', alpha=min(1.0,size/6), edgecolors='None', linewidths=0)

			ax.set_xlim(-params['L']/2,params['L']/2)
			ax.set_ylim(-params['L']/2,params['L']/2)
			ax.set_zlim(-params['H']/2,params['H']/2)

			outputPath = path + 'Visualize3D_py/Phages/%03i.png' % i

			# plt.tight_layout()
			ax.set_xticks(xyticks)
			ax.set_yticks(xyticks)
			ax.set_zticks(zticks)
			ax.set_xticklabels(xlabels,fontsize=8, verticalalignment='top', **font)
			ax.set_yticklabels(ylabels,fontsize=8, verticalalignment='top', **font)
			ax.set_zticklabels(zlabels,fontsize=8, verticalalignment='top', **font)
			ax.set_xlabel('x (µm)',fontsize=12,labelpad=10, **font)
			ax.set_ylabel('y (µm)',fontsize=12,labelpad=8, **font)
			ax.set_zlabel('z (µm)',fontsize=12,labelpad=12, **font)
			plt.savefig(outputPath, dpi=300)
			plt.cla()

			if i == 7:

				destPath = cwd + '/Fig_1'
				if not os.path.isdir(destPath) :
					os.mkdir(destPath)

				sourcePath = path + 'Visualize3D_py/Phages/%03i.png' % i
				shutil.copy2(sourcePath,destPath + '/Snapshot_Phages.png')

				sourcePath = path + 'Visualize3D_py/Cells_density/%03i.png' % i
				shutil.copy2(sourcePath,destPath + '/Snapshot_Cells.png')

			# Increase counter
			i = i + 1
		else :
			print("Path already exists: " + outputPath)
			# Increase counter
			i = i + 1
	plt.close()
	os.chdir(cwd)


# Start the program
Visualize3D()
