#!/bin/python

import sys
import math
from collections import namedtuple, defaultdict
from itertools import product
import drawSvg as draw # pip install drawSvg


Run = namedtuple('Run', ['tp', 'fn', 'fp', 'tn', 'time', 'sensitivity', 'precision', 'tool', 'rank', 'data'])

def read_csv(infile):
	for line in open(infile):
		line = line.strip()
		if not line: continue
		data, rank, tool, real, tp, fp, tn, fn = line.split(',')
		tp, fn, fp, tn = int(tp), int(fn), int(fp), int(tn)
		m, s = real.split('m')
		time = (int(m) * 60 + float(s[:-1])) / 60
		sensitivity = tp / (tp + fn)
		precision = tp / (tp + fp)
		yield Run(tp, fn, fp, tn, time, sensitivity, precision, tool, rank, data)

def show(imgfile, rank, infile='combined-metabenchmark.tsv'):
	colours = { 'kaiju': '#4e79a7'
	          , 'clark': '#f28e2b'
	          , 'kraken': '#e15759'
	          , 'umgap high precision': '#76b7b2'
	          , 'umgap high sensitivity': '#59a14f'
	          , 'umgap max precision': '#edc948'
	          , 'umgap max sensitivity': '#b07aa1'
	          , 'umgap tryptic precision': '#ff9da7'
	          , 'umgap tryptic sensitivity': '#9c755f'
	          , 'kraken2': '#bab0ac'
	          }
	scale = 1000
	bottom_cutoff = 0.65
	d = draw.Drawing(scale, scale - bottom_cutoff * scale, viewport_fill='white')
	d.append(draw.Rectangle(0, 0, scale, scale - bottom_cutoff * scale, fill='white'))
	for run in read_csv(infile):
		if run.rank != rank: continue
		# p = draw.Path(stroke_width=2, stroke=colours[run.tool])
		# p.M(scale * run.sensitivity, scale * run.precision - bottom_cutoff * scale)
		# d.append(p)
		d.append(draw.Circle(scale * run.sensitivity, scale * run.precision - bottom_cutoff * scale, 4, fill_opacity=0, stroke_width=2, stroke=colours[run.tool]))
	d.saveSvg(imgfile)

tools = [ 'umgap tryptic precision'
        , 'kraken'
        , 'umgap max precision'
        , 'kraken2'
        , 'umgap high precision'
        , 'kaiju'
        , 'umgap tryptic sensitivity'
        , 'umgap high sensitivity'
        , 'umgap max sensitivity'
        , 'clark'
        ]

def averages(rank, infile='combined-metabenchmark.tsv'):
	ppv = defaultdict(list) # positive predictive values (precision)
	tpr = defaultdict(list) # true positive rates (sensitivity)
	tnr = defaultdict(list) # true negative rates (specificity)
	npv = defaultdict(list) # negative predictive values
	mcc = defaultdict(list) # matthew's correlation coefficient
	times = defaultdict(list)
	for run in read_csv(infile):
		if run.rank != rank: continue
		tp, fp, tn, fn = run.tp, run.fp, run.tn, run.fn
		ppv[run.tool].append(run.precision)
		tpr[run.tool].append(run.sensitivity)
		tnr[run.tool].append(tn / (tn + fp))
		npv[run.tool].append(tn / (tn + fn))
		mcc[run.tool].append((tp * tn - fp * fn) / math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
		times[run.tool].append(run.time)
	print('ppv', 'tpr', 'tnr', 'npv', 'mcc', 'time', 'tool', sep='\t')
	for tool in tools:
		if ppv[tool]:
			print(
				round(100 * sum(ppv[tool])/len(ppv[tool]), 2),
				round(100 * sum(tpr[tool])/len(tpr[tool]), 2),
				round(100 * sum(tnr[tool])/len(tnr[tool]), 2),
				round(100 * sum(npv[tool])/len(npv[tool]), 2),
				round(100 * sum(mcc[tool])/len(mcc[tool]), 2),
				round(sum(times[tool])/len(times[tool]), 2),
				tool,
				sep='\t'
			)

def averages(rank, data=None, infile='combined-metabenchmark.tsv'):
	ppv = defaultdict(list) # positive predictive values (precision)
	tpr = defaultdict(list) # true positive rates (sensitivity)
	tnr = defaultdict(list) # true negative rates (specificity)
	npv = defaultdict(list) # negative predictive values
	mcc = defaultdict(list) # matthew's correlation coefficient
	times = defaultdict(list)
	for run in read_csv(infile):
		if data and run.data != data: continue
		if run.rank != rank: continue
		tp, fp, tn, fn = run.tp, run.fp, run.tn, run.fn
		ppv[run.tool].append(run.precision)
		tpr[run.tool].append(run.sensitivity)
		tnr[run.tool].append(tn / (tn + fp))
		npv[run.tool].append(tn / (tn + fn))
		mcc[run.tool].append((tp * tn - fp * fn) / math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
		times[run.tool].append(run.time)
	print('ppv', 'tpr', 'tnr', 'npv', 'mcc', 'time', 'tool', sep='\t')
	for tool in tools:
		if ppv[tool]:
			print(
				f'{round(100 * sum(ppv[tool])/len(ppv[tool]), 2)}%',
				f'{round(100 * sum(tpr[tool])/len(tpr[tool]), 2)}%',
				f'{round(100 * sum(tnr[tool])/len(tnr[tool]), 2)}%',
				f'{round(100 * sum(npv[tool])/len(npv[tool]), 2)}%',
				f'{round(100 * sum(mcc[tool])/len(mcc[tool]), 2)}%',
				f'{round(sum(times[tool])/len(times[tool]), 2)}m',
				tool,
				sep='\t'
			)

