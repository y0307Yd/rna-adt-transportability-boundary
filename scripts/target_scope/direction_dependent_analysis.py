import csv, sys, io, statistics, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

base = r'D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap'

def load_csv(path):
    with open(path, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

bidir = load_csv(base + r'\p3_target_scope_20260902\p3_transport_expanded43_v1\P3_EXPANDED43_BIDIRECTIONAL_EVIDENCE_MATRIX.csv')
cls = load_csv(base + r'\p3_target_scope_20260902\p3_transport_expanded43_v2_classified\P3_43TARGET_CLASSIFIED_MAIN_RESULT.csv')
xw = load_csv(base + r'\technical_biological_decomposition_data_v1_20260903\GSE245108_tnc_parsed\four_dataset_crosswalk_kotliarov_gse245108.csv')

# index
cls_by_target = {r['target']: r for r in cls}
xw_by_target = {}
for r in xw:
    t = r['target']
    if t not in xw_by_target:
        xw_by_target[t] = r  # first occurrence (dedup HLA-DR)

# integrated classification counts
from collections import Counter, defaultdict
ic = Counter()
groups = defaultdict(list)
for r in bidir:
    g = r['integrated_classification']
    ic[g] += 1
    groups[g].append(r['target'])

print('=== integrated_classification counts ===')
for k, v in ic.items():
    print(f'{k}: {v}')

print()
print('=== direction_dependent targets ===')
print(groups.get('direction_dependent', []))
print()
print('=== bidirectionally_stable ===')
print(groups.get('bidirectionally_stable', []))
print()
print('=== weak_or_non_transportable ===')
print(groups.get('weak_or_non_transportable', []))

# feature comparison: direction_dependent vs bidirectionally_stable vs weak
def feat(target, col, table, cast=float):
    if target in table:
        try:
            v = table[target].get(col)
            if v in (None, ''):
                return None
            return cast(v)
        except (ValueError, TypeError):
            return None
    return None

def summarize(targets, label):
    rows = {
        'n': len(targets),
        'spearman': [], 'donor_median': [], 'kot_ICC': [], 'frac_donor_kot': [],
    }
    for t in targets:
        s = feat(t, 'spearman', cls_by_target)
        dm = feat(t, 'donor_median', cls_by_target)
        icc = feat(t, 'kot_ICC', xw_by_target)
        fd = feat(t, 'frac_donor_kot', xw_by_target)
        if s is not None: rows['spearman'].append(s)
        if dm is not None: rows['donor_median'].append(dm)
        if icc is not None: rows['kot_ICC'].append(icc)
        if fd is not None: rows['frac_donor_kot'].append(fd)
    print(f'\n=== {label} (n={rows["n"]}) ===')
    for k, v in rows.items():
        if k == 'n':
            continue
        if v:
            print(f'  {k}: median={statistics.median(v):.3f}  n={len(v)}')
        else:
            print(f'  {k}: (no data)')
    return rows

dd = summarize(groups.get('direction_dependent', []), 'direction_dependent')
bs = summarize(groups.get('bidirectionally_stable', []), 'bidirectionally_stable')
wk = summarize(groups.get('weak_or_non_transportable', []), 'weak_or_non_transportable')

# Mann-Whitney U for direction_dependent vs bidirectionally_stable on key features
from math import sqrt
def mwu(a, b):
    a = [x for x in a if x is not None]
    b = [x for x in b if x is not None]
    if not a or not b:
        return None, None
    from scipy.stats import mannwhitneyu
    u, p = mannwhitneyu(a, b, alternative='two-sided')
    return u, p

print('\n=== Mann-Whitney: direction_dependent vs bidirectionally_stable ===')
for featname in ['spearman', 'donor_median', 'kot_ICC', 'frac_donor_kot']:
    a = [feat(t, featname, cls_by_target) for t in groups.get('direction_dependent', [])]
    b = [feat(t, featname, cls_by_target) for t in groups.get('bidirectionally_stable', [])]
    a = [x for x in a if x is not None]
    b = [x for x in b if x is not None]
    if a and b:
        u, p = mwu(a, b)
        print(f'  {featname}: DD median={statistics.median(a):.3f} (n={len(a)}) vs BS median={statistics.median(b):.3f} (n={len(b)})  U p={p:.4f}')
    else:
        print(f'  {featname}: insufficient data')

# also compare repeatability from crosswalk for direction-dependent targets
print('\n=== direction-dependent targets: kot_ICC (repeatability) detail ===')
for t in sorted(groups.get('direction_dependent', [])):
    icc = feat(t, 'kot_ICC', xw_by_target)
    fd = feat(t, 'frac_donor_kot', xw_by_target)
    s = feat(t, 'spearman', cls_by_target)
    print(f'  {t}: kot_ICC={icc if icc is not None else "NA"}, frac_donor={fd if fd is not None else "NA"}, spearman={s if s is not None else "NA"}')
