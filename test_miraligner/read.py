from collections import defaultdict

pre = defaultdict(str)
with open('long.fa') as in_handle:
    for line in in_handle:
        if line.startswith(">"):
            name = line.strip()[1:]
        else:
            pre[name] += line.strip()

print pre
