
import sys

def main(lineages_file):
    with open(lineages_file) as f:
        lineages = { line.split(',')[0]: line for line in f }
    for line in sys.stdin:
        if line.startswith('>'):
            print(line, end='')
        else:
            print(lineages[line.split(',')[0]], end='')

if __name__ == "__main__":
    main(*sys.argv[1:])

