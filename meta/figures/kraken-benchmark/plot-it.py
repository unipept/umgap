#!/bin/python

import sys
from math import sqrt
from collections import namedtuple, defaultdict
from itertools import product
from functools import partial
import drawSvg as draw # pip install drawSvg

class Run:
	__slots__ = [ "exectime_s", "translator", "aggregator_p", "minfreq",
		"fragmentor", "maxgap_mintryp", "minseed_maxtryp", "identified_reads",
		"fp", "tn", "fn", "tp", "sensitivity", "precision", "specificity",
		"npv", "mcc" ]

	def __init__(self, realtot, d, dt, p, a, mf, f, a1, a2, tot, fp, tn, fn, tp):
		self.exectime_s = int(dt)
		self.translator = p
		self.aggregator_p = int(a)
		self.minfreq = int(mf)
		self.fragmentor = f
		self.maxgap_mintryp = int(a1) if a1 else None
		self.minseed_maxtryp = int(a2) if a2 else None
		self.identified_reads = int(tot)
		fp = int(fp)
		tn = int(tn)
		fn = int(fn) + max(0, realtot - int(tot)) # assuming no true negatives
		tp = int(tp)
		self.fp = fp
		self.tn = tn
		self.fn = fn
		self.tp = tp
		self.sensitivity = tp / max(0.0001, tp + fn)
		self.precision = tp / max(0.0001, tp + fp)
		self.specificity = tn / max(0.0001, tn + fp)
		self.npv = tn / max(0.0001, tn + fn)
		self.mcc = (tp*tn - fp*fn) / max(0.0001, sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn)))

	def __repr__(self):
		return " ".join(f"{k}={getattr(self, k)}" for k in Run.__slots__)

def parse_results(infile, reads):
	with open(infile) as f:
		for line in f:
			yield Run(reads, *line.split("\t"))

def show(imgfile, colour_function, infile='2021-06-02.HiSeq.output', reads=10000, b=0.0, t=1.0, l=0.0, r=1.0):
	width = 1000
	height = 700
	d = draw.Drawing(width, height, viewport_fill='white')
	d.append(draw.Rectangle(0, 0, width, height, fill='white'))
	for run in parse_results(infile, reads):
		colour = colour_function(run)
		if not colour: continue
		x = width * (run.sensitivity - l) / (r - l)
		y = height * (run.precision - b) / (t - b)
		if x >= 0.001:
			if y >= 0.001:
				d.append(draw.Circle(x, y, 4, fill_opacity=0, stroke_width=2, stroke=colour))
			else:
				d.append(draw.Lines(x - 4, 8, x, 0, x + 4, 8, fill_opacity=0, stroke_width=1, stroke=colour))
		else:
			if y >= 0.001:
				d.append(draw.Lines(8, y - 4, 0, y, 8, y + 4, fill_opacity=0, stroke_width=1, stroke=colour))
			else:
				d.append(draw.Lines(8, 4, 0, 0, 4, 8, fill_opacity=0, stroke_width=1, stroke=colour))
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

def fragmentor(run):
	return blue if run.fragmentor == 'tryp' else orange

def tryptic_translator(run):
	if run.fragmentor != 'tryp': return None
	return dict(ft6=blue, fgsrs=orange, fgs=cyan).get(run.translator)

def tryptic_length(run):
	if run.fragmentor != 'tryp': return None
	return colours[run.maxgap_mintryp - 5]

def tryptic_freq(run):
	if run.fragmentor != 'tryp': return None
	return colours[run.minfreq - 1]

def kmer_seedornot(run):
	return dict(kmer=blue, seed=orange).get(run.fragmentor)

def seedextend_translator(run):
	if run.fragmentor != 'seed': return None
	return dict(ft6=blue, fgsrs=orange, fgs=cyan).get(run.translator)

def seedextend_translator_seedsize(run):
	if run.fragmentor != 'seed': return None
	if run.translator == 'ft6': return colours[run.minseed_maxtryp - 2]
	if run.translator == 'fgsrs': return colours[run.minseed_maxtryp + 1]
	return None

def seedextend_settranslator_seedsize(translator, run):
	if run.fragmentor != 'seed': return None
	if run.translator != translator: return None
	return colours[run.minseed_maxtryp + 2]

def seedextend_settranslator_freq(translator, run):
	if run.fragmentor != 'seed': return None
	if run.translator != translator: return None
	return colours[run.minfreq - 1]

def seedextend_settranslator_profiler(translator, run):
	if run.fragmentor != 'seed': return None
	if run.translator != translator: return None
	return {
		100: blue, # LCA*
		75: orange, # hybrid -f 0.75
		50: red, # hybrid -f 0.5
		25: cyan, # hybrid -f 0.25
		0: green # mrtl
	}.get(run.aggregator_p)

def main():
	plot(fragmentor, b=0.3, t=1.0, l=0.0, r=1.0) # 0.3 and 0.0
	plot(tryptic_translator, b=0.3, t=1.0, l=0.0, r=0.3)
	plot(tryptic_length, b=0.3, t=1.0, l=0.0, r=0.3)
	plot(tryptic_freq, b=0.3, t=1.0, l=0.0, r=0.3)
	plot(kmer_seedornot, b=0.5, t=1.0, l=0.0, r=1.0)
	plot(seedextend_translator, b=0.9, t=1.0, l=0.4, r=1.0)
	plot(seedextend_translator_seedsize, b=0.9, t=1.0, l=0.4, r=1.0)
	plot(partial(seedextend_settranslator_seedsize, 'ft6'), name='seedextend_6ft_seedsize', b=0.9, t=1.0, l=0.4, r=1.0)
	plot(partial(seedextend_settranslator_seedsize, 'fgsrs'), name='seedextend_fgsrs_seedsize', b=0.9, t=1.0, l=0.4, r=1.0)
	plot(partial(seedextend_settranslator_freq, 'ft6'), name='seedextend_6ft_freq', b=0.9, t=1.0, l=0.4, r=1.0)
	plot(partial(seedextend_settranslator_freq, 'fgsrs'), name='seedextend_fgsrs_freq', b=0.9, t=1.0, l=0.4, r=1.0)
	plot(partial(seedextend_settranslator_profiler, 'ft6'), name='seedextend_6ft_profiler', b=0.9, t=1.0, l=0.4, r=1.0)
	plot(partial(seedextend_settranslator_profiler, 'fgsrs'), name='seedextend_fgsrs_profiler', b=0.9, t=1.0, l=0.4, r=1.0)
