import re
from datetime import datetime
from collections import defaultdict

lines =[]
while True:
    try:
        line = input().strip()
        if line:
            lines.append(line)
    except EOFError:
        break

rules = {}

rulesLine = lines[0]

for part in rulesLine.split(','):
    if part.strip():
        k, v = part.strip().split('=')
        rules[k] = int(v)



desLine = lines[1]
fieldNames = [f.strip() for f in desLine.split(',')]



patternsD = {
    'Host': r'(\S+)',
    'Client IP': r'(\S+)',
    'Id': r'(\S+)',
    'Date': r'\[([^\]]+)\]',
    'Request': r'"([^"]*)"',
    'HTTP Status': r'(\S+)',
    'User Agent': r'"([^"]*)"',
    'Session Cookie': r'"([^"]*)"',
}



patternStr = r'^' + r'\s+'.join(patternsD[f] for f in fieldNames) + r'$'
re_pat = re.compile(patternStr)

log_data = defaultdict(lambda: defaultdict(lambda: {
    'ips': set(),
    'agents': set(),
    'sessions': set(),
    'pdf_count': 0,
    'pdf_ids': []
}))

for line in lines[2:]:
    match =re_pat.match(line)
    if not match:
        continue
    values = match.groups()
    log= {fieldNames[i]: values[i] for i in range(len(fieldNames))}
    
    status = log['HTTP Status']
    if status != '200':
        continue
    
    id_ = log['Id']
    if id_ == '-':
        continue
    date_str= log['Date']
    try:
        dt = datetime.strptime(date_str, '%d/%b/%Y:%H:%M:%S')
        day =dt.date()
    except ValueError:
        continue
    
    ip= log['Client IP']
    agent =log['User Agent']
    session= log['Session Cookie']
    request = log['Request']
    
    day_data= log_data[id_][day]
    day_data['ips'].add(ip)
    day_data['agents'].add(agent)
    day_data['sessions'].add(session)
    reqParts= request.split()
    if len(reqParts) == 3 and reqParts[0] == 'GET' and reqParts[2] == 'HTTP/1.1':
        url = reqParts[1]
        if url.startswith('/document/') and url.endswith('.pdf'):
            doc_str = url[10:-4]
            if doc_str.isdigit():
                day_data['pdf_count'] += 1
                day_data['pdf_ids'].append(int(doc_str))

violations = []
for id_ in sorted(log_data):
    for day in log_data[id_]:
        d = log_data[id_][day]
        
        if 'agent' in rules:
            cnt = len(d['agents'])
            if cnt >= rules['agent']:
                violations.append((id_, 'agent', cnt))
        
        if 'ip' in rules:
            cnt= len(d['ips'])
            if cnt >= rules['ip']:
                violations.append((id_, 'ip', cnt))
        
        if 'pdf' in rules:
            cnt= d['pdf_count']
            if cnt >= rules['pdf']:
                violations.append((id_, 'pdf', cnt))
        
        if 'session' in rules:
            cnt= len(d['sessions'])
            if cnt >= rules['session']:
                violations.append((id_, 'session', cnt))
        
        if 'crawl' in rules:
            idsList =d['pdf_ids']
            if len(idsList) == 0:
                max_streak = 0
            else:
                max_streak= 1
                current= 1
                for i in range(1, len(idsList)):
                    if idsList[i] == idsList[i-1] + 1:
                        current += 1
                        max_streak = max(max_streak, current)
                    else:
                        current = 1
            if max_streak >= rules['crawl']:
                violations.append((id_, 'crawl', max_streak))
violations.sort(key=lambda x: (x[0], x[1]))

if not violations:
    print('N/A')
else:
    for v in violations:
        print(f"{v[0]} {v[1]}={v[2]}")
# Authors: Morteza Farrokhnejad, Ali Farrokhnejad