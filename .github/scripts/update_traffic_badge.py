import os
import requests
import json

token = os.environ['GH_TOKEN']
repo = os.environ['REPO']

headers = {
    'Authorization': f'token {token}',
    'Accept': 'application/vnd.github.v3+json'
}

url = f'https://api.github.com/repos/{repo}/traffic/clones'

response = requests.get(url, headers=headers)

if response.status_code == 200:
    data = response.json()
    count = data.get("count", 0)
    uniques = data.get("uniques", 0)

    badge = {
        "schemaVersion": 1,
        "label": "clones (14d)",
        "message": f"{count} total / {uniques} unique",
        "color": "blue"
    }

    with open("traffic-badge.json", "w") as f:
        json.dump(badge, f)

else:
    print("Failed to fetch traffic data:", response.status_code)
    print(response.text)
    exit(1)
