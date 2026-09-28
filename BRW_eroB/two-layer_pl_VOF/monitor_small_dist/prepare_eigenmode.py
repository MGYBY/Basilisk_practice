#!/usr/bin/env python3
"""Validate immutable full-NS export provenance and install its Eulerian mode."""
import argparse,configparser,hashlib,json,math,re
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('source',type=Path)
p.add_argument('--base-state',type=Path,default=Path('base_state'))
a=p.parse_args()
m=json.loads((a.source/'eigenmode.json').read_text())
if m.get('schema')!='eulerian_eigenmode_v1':
    raise ValueError('Unrecognized eigenmode schema')
if 'input_snapshot' not in m or 'profile_sha256' not in m:
    raise ValueError('Export lacks immutable input/profile provenance; regenerate it with the supplied FullNS/export_eigenmode.py')
current_path=a.base_state/'used_case_parameters.ini'
if not current_path.exists():
    raise ValueError('First generate the target base state so used_case_parameters.ini can be checked')
c=configparser.ConfigParser(inline_comment_prefixes=('#',';'))
c.read(current_path)
for section in ('physical_problem','lower_rheology','upper_rheology'):
    expected=m['input_snapshot'].get(section)
    if not isinstance(expected,dict) or not expected:
        raise ValueError('Missing immutable input snapshot section: '+section)
    for key,value in expected.items():
        actual=c.get(section,key,fallback=None)
        if actual is None:raise ValueError(f'Missing current parameter {section}.{key}')
        try:
            matches=math.isclose(float(value),float(actual),rel_tol=1e-12,abs_tol=1e-14)
        except (ValueError,TypeError):
            matches=str(value).strip().lower()==str(actual).strip().lower()
        if not matches:
            raise ValueError(f'Eigenmode/current physical mismatch: {section}.{key}: {value} vs {actual}')
derived=m['input_snapshot'].get('derived_materials')
if derived:
    header=(a.base_state/'generated_case.h').read_text()
    for name,macro in [('lambda_lower','CASE_LAMBDA_LOWER'),('lambda_upper','CASE_LAMBDA_UPPER'),('air_eta','CASE_AIR_ETA')]:
        match=re.search(r'^#define\s+'+macro+r'\s+([^\s]+)',header,re.M)
        if not match or not math.isclose(float(match.group(1)),float(derived[name]),rel_tol=1e-11,abs_tol=1e-14):
            raise ValueError('Eigenmode derived material coefficient differs from generated case: '+name)
dm=c['domain_and_mesh']
length=float(dm['lx']) if 'lx' in dm else float(dm['lx_slope_scaled'])/c.getfloat('physical_problem','slope_tan')
if not math.isclose(length,float(m['actual_top']),rel_tol=1e-12):
    raise ValueError('Eigenmode outer top/domain mismatch')
mode=float(m['k_star'])*length/(2*np.pi)
for key in ('perturbation_mode','tracked_mode'):
    if not math.isclose(c.getint('initial_condition',key),mode,rel_tol=0,abs_tol=1e-10):
        raise ValueError('Eigenmode wavenumber must match both perturbation_mode and tracked_mode')
vals=[m[k] for k in ('k_star','surface_amplitude','eta_interface_real','eta_interface_imag','eta_surface_real','eta_surface_imag')]
if not np.all(np.isfinite(vals)) or not 0<vals[1]<=.01:
    raise ValueError('Invalid eigenmode metadata/amplitude')
etaI=complex(vals[2],vals[3]);etaS=complex(vals[4],vals[5]);depth=c.getfloat('physical_problem','depth_ratio')
if vals[1]*abs(etaI)>=1 or vals[1]*abs(etaS-etaI)>=depth or vals[1]*abs(etaS)>=length-1-depth:
    raise ValueError('The displaced phase thickness can become nonpositive')
# Validate all source data before changing any installed files. The mutable
# metadata.config path is provenance text only; it never authorizes compatibility.
phase_data={}
for phase,(lo,hi) in zip(('lower','upper','air'),((0,1),(1,1+depth),(1+depth,length))):
    filename=f'mode_{phase}.tsv';source=a.source/filename
    recorded=m['profile_sha256'].get(filename)
    if not recorded or hashlib.sha256(source.read_bytes()).hexdigest()!=recorded:
        raise ValueError('Eigenmode profile hash mismatch: '+filename)
    data=np.loadtxt(source)
    if (data.ndim!=2 or data.shape[1]!=11 or data.shape[0]<2 or not np.all(np.isfinite(data))
            or not np.all(np.diff(data[:,0])>0) or abs(data[0,0]-lo)>1e-10 or abs(data[-1,0]-hi)>1e-10):
        raise ValueError('Invalid eigenmode phase table/domain: '+phase)
    phase_data[phase]=data
a.base_state.mkdir(parents=True,exist_ok=True)
for phase,data in phase_data.items():
    target=a.base_state/f'eigenmode_{phase}.dat';staged=target.with_suffix('.dat.tmp')
    with staged.open('w') as f:
        f.write(str(len(data))+'\n');np.savetxt(f,data,fmt='%.17g')
    staged.replace(target)
(a.base_state/'eigenmode_metadata.dat').write_text(' '.join(format(v,'.17g') for v in vals)+'\n')
(a.base_state/'eigenmode_source.json').write_text(json.dumps(m,indent=2)+'\n')
print('Installed verified eigenmode:',a.source,'->',a.base_state)
