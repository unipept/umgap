#!/bin/python

import sys
import math
from collections import namedtuple, defaultdict
from itertools import product
from functools import partial
import drawSvg as draw # pip install drawSvg


Run = namedtuple('Run', ['tp', 'fn', 'fp', 'tn', 'time', 'sensitivity', 'precision', 'tool'])

def chunks(infile):
	chunk = {}
	with open(infile) as stream:
		for line in stream:
			line = line.strip()
			if not line: pass
			elif 'FN' in line: chunk['fn'] = int(line.split()[0])
			elif 'FP' in line: chunk['fp'] = int(line.split()[0])
			elif 'total' in line: chunk['total'] = int(line.split()[0])
			elif 'TP' in line: chunk['tp'] = int(line.split()[0])
			elif 'real' in line: chunk['real'] = float(line.split()[-1][2:-1])
			elif 'user' in line: pass
			elif 'sys' in line:
				yield chunk
				chunk = {}
			else: chunk['name'] = line

def parse_results(infile):
	for chunk in chunks(infile):
		name = chunk['name']
		tp = chunk.get('tp', 0)
		#tn = chunk.get('tn', 0)
		fn = chunk.get('fn', 0)
		fp = chunk.get('fp', 0)
		real = chunk['real']
	
		sensitivity = tp / (tp + fn)
		precision = tp / (tp + fp)
		#specificity = tn / (tn + fp)
		#npv = tn / (tn + fn)
		#mcc = (tp*tn - fp*fn) / sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
	
		yield Run(tp, fn, fp, 0, real, sensitivity, precision, name)

def show(imgfile, colour_function, infile='HiSeq.benchmark.output', b=0.0, t=1.0, l=0.0, r=1.0):
	width = 1000
	height = 700
	d = draw.Drawing(width, height, viewport_fill='white')
	d.append(draw.Rectangle(0, 0, width, height, fill='white'))
	for run in parse_results(infile):
		colour = colour_function(run.tool)
		if not colour: continue
		d.append(draw.Circle(width * (run.sensitivity - l) / (r - l), height * (run.precision - b) / (t - b), 4, fill_opacity=0, stroke_width=2, stroke=colour))
	d.saveSvg(imgfile)

def plot(colour_function, name=None, **kwargs):
	show((name or colour_function.__name__) + '.svg', colour_function, **kwargs)

colours = [ blue := '#4e79a7'
          , purple := '#b07aa1'
          , red := '#e15759'
          , cyan := '#76b7b2'
          , green := '#59a14f'
          , yellow := '#edc948'
          , orange := '#f28e2b'
          , pink := '#ff9da7'
          , brown := '#9c755f'
          , gray := '#bab0ac'
          ]

def digestor(name):
	return blue  if 'tryptic' in name else orange

def tryptic_translator(name):
	if 'tryptic' not in name: return None
	if 'ft6' in name: return blue
	if 'fgspp' in name: return orange
	if 'fgs' in name: return cyan
	return None

def tryptic_length(name):
	if 'tryptic' not in name: return None
	return colours[int(name.split()[-2]) - 5]

def tryptic_freq(name):
	if 'tryptic' not in name: return None
	return colours[int(name.split()[-4]) - 1]

def kmer_seedornot(name):
	if 'kmer' in name: return blue
	if 'seedextend' in name and 'scored' not in name: return orange
	return None

def seedextend_translator(name):
	if 'seedextend' not in name or 'scored' in name: return None
	if 'ft6' in name: return blue
	if 'fgspp' in name: return orange
	if 'fgs' in name: return cyan
	return None

def seedextend_translator_seedsize(name):
	if 'seedextend' not in name or 'scored' in name: return None
	if 'ft6' in name: return colours[int(name.split()[-1]) - 2]
	if 'fgs ' in name: return colours[int(name.split()[-1]) + 1]
	return None

def seedextend_settranslator_seedsize(translator, name):
	if 'seedextend' not in name or 'scored' in name: return None
	if translator not in name: return None
	return colours[int(name.split()[-1]) + 2]

def seedextend_settranslator_freq(translator, name):
	if 'seedextend' not in name or 'scored' in name: return None
	if translator not in name: return None
	return colours[int(name.split()[-4]) - 1]

def seedextend_settranslator_profiler(translator, name):
	if 'seedextend' not in name or 'scored' in name: return None
	if translator not in name: return None
	if 'LCA*' in name or 'lca*' in name: return blue
	if 'hybrid -f 0.75' in name: return orange
	if 'hybrid -f 0.5' in name: return red
	if 'hybrid -f 0.25' in name: return cyan
	if 'mrtl' in name: return green
	return None

def main():
	plot(digestor, b=0.5, t=1.0, l=0.0, r=1.0)
	plot(tryptic_translator, b=0.5, t=1.0, l=0.0, r=0.6)
	plot(tryptic_length, b=0.5, t=1.0, l=0.0, r=0.6)
	plot(tryptic_freq, b=0.5, t=1.0, l=0.0, r=0.6)
	plot(kmer_seedornot, b=0.5, t=1.0, l=0.0, r=1.0)
	plot(seedextend_translator, b=0.9, t=1.0, l=0.4, r=1.0)
	plot(seedextend_translator_seedsize, b=0.9, t=1.0, l=0.4, r=1.0)
	plot(partial(seedextend_settranslator_seedsize, 'ft6'), name='seedextend_6ft_seedsize', b=0.9, t=1.0, l=0.4, r=1.0)
	plot(partial(seedextend_settranslator_seedsize, 'fgs '), name='seedextend_fgs_seedsize', b=0.9, t=1.0, l=0.4, r=1.0)
	plot(partial(seedextend_settranslator_freq, 'ft6'), name='seedextend_6ft_freq', b=0.9, t=1.0, l=0.4, r=1.0)
	plot(partial(seedextend_settranslator_freq, 'fgs '), name='seedextend_fgs_freq', b=0.9, t=1.0, l=0.4, r=1.0)
	plot(partial(seedextend_settranslator_profiler, 'ft6'), name='seedextend_6ft_profiler', b=0.9, t=1.0, l=0.4, r=1.0)
	plot(partial(seedextend_settranslator_profiler, 'fgs '), name='seedextend_fgs_profiler', b=0.9, t=1.0, l=0.4, r=1.0)
