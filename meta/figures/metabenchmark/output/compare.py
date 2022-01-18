#!/usr/bin/env python3

import sys

from collections import Counter


def main(abundances, *outputs):

	abundances = Counter(dict(map(int, l.strip().split()) for l in open(abundances)))
	total = sum(abundances.values())

	print('total', total)

	print('tool', 'TP', 'FP', 'TN', 'FN', sep=',')

	for output in outputs:
		with open(output) as f:
			predictions = Counter(int(l.strip()) for l in f)

		tp = 0
		fp = 0
		for organism in predictions.keys() - {0, 1}:
			if abundances[organism] >= predictions[organism]:
				tp += predictions[organism]
			else:
				tp += abundances[organism]
				fp += predictions[organism] - abundances[organism]

		expected_negatives = abundances[0] + abundances[1]
		received_negatives = total - tp - fp
		tn = min(expected_negatives, received_negatives)
		fn = received_negatives - tn

		print(output, tp, fp, tn, fn, tp+fp+tn+fn, sep=',')


if __name__ == '__main__':
	main(*sys.argv[1:])
