import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import gzip,json,hashlib
B=Path(os.path.join(os.environ['PD_PROJECT_ROOT'], 'batch15'))
plan=json.loads((B/'protocol/analysis_plan_v1.json').read_text())
O=B/'raw/lpi_regions_v1';O.mkdir(exist_ok=False)
def scan(reg):
 src=B.parent/'batch12/raw/gwas'/f"{reg['study']}.tsv.gz"; dst=O/f"{reg['id']}.tsv.gz"
 before=src.stat();total=kept=0;lead=[]
 with gzip.open(src,'rt') as f,gzip.open(dst,'wt') as out:
  header=f.readline();columns=header.rstrip().split('\t');assert columns[:2]==['chromosome','base_pair_location'];out.write(header)
  for line in f:
   total+=1
   if total<=20000:
    fields=line.rstrip().split('\t');assert len(fields)==len(columns) and float(fields[5])>0
   if line.startswith(reg['chromosome']+'\t'):
    pos=int(line.split('\t',2)[1])
    if reg['start']<=pos<=reg['end']:
     assert len(line.rstrip().split('\t'))==len(columns)
     out.write(line);kept+=1
     if pos==reg['lead_position']:lead.append(line.rstrip())
 after=src.stat();assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)
 assert total>11000000 and kept>0 and lead
 result={'region':reg['id'],'source':str(src),'source_rows_scanned':total,'region_rows':kept,'lead_rows':lead,'source_stat_unchanged':True,'output_sha256':hashlib.sha256(dst.read_bytes()).hexdigest()}
 (O/f"{reg['id']}.manifest.json").write_text(json.dumps(result,indent=2));return result
with ProcessPoolExecutor(max_workers=2) as pool:results=list(pool.map(scan,plan['regions']))
(O/'COMPLETE.json').write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2))
