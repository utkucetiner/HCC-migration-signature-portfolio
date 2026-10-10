# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
import gzip,re,json
from pathlib import Path
from collections import defaultdict
hits=defaultdict(list)
with gzip.open('work/data/gencode.v38.annotation.gtf.gz','rt') as f:
 for line in f:
  if line.startswith('#'):continue
  a=line.rstrip().split('\t')
  if a[2]!='gene' or a[0] not in ['chr'+str(i) for i in range(1,23)]:continue
  gene=re.search(r'gene_name "([^"]+)"',a[8]).group(1);hits[gene].append((a[0],int(a[3]),int(a[4])))
genes=Path('work/matrix_gene_ids.txt').read_text().splitlines();alias={'LINC00152':'CYTOR','SGK223':'PRAG1'};rows=[];ambiguous=[]
for gene in genes:
 h=hits.get(gene,hits.get(alias.get(gene,''),[]))
 if len(h)==1:rows.append((gene,*h[0]))
 elif len(h)>1:ambiguous.append(gene)
rows.sort(key=lambda r:(int(r[1][3:]),r[2],r[0]))
Path('work/cnv/gene_order.tsv').write_text(''.join('\t'.join(map(str,r))+'\n' for r in rows));Path('work/cnv/gene_order_audit.json').write_text(json.dumps({'source':'GENCODE v38 GRCh38','autosomal_mapped':len(rows),'ambiguous_excluded':ambiguous,'total_genes':len(genes),'aliases':alias},indent=2));print(len(rows),len(ambiguous))
