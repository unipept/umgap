import requests
import sys


url = 'http://www.uniprot.org/mapping/'


for line in sys.stdin:
    without_version=line.split(".", 1)[0]

    params = {
        'from': 'EMBL_ID',
        'to': 'ACC',
        'format': 'list',
        'query': without_version
    }

    r = requests.get(url, params=params)
    print(r.text, end="")
