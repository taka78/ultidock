import os
import requests
import json

token = os.environ.get('GH_TOKEN', '').strip()
if not token:
    raise SystemExit('Set GH_TRAFFIC_TOKEN to a token with repository traffic access.')
repo = os.environ['REPO']

headers = {
    'Authorization': f'Bearer {token}',
    'Accept': 'application/vnd.github+json'
}

url = f'https://api.github.com/repos/{repo}/traffic/clones'

response = requests.get(url, headers=headers, timeout=30)

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
