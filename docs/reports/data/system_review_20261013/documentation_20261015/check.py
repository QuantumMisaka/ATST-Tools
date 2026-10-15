from pathlib import Path
import os, shutil, subprocess, sys, json
from ase.io import read
root=Path('/home/james/.codex/worktrees/3fbf/atst-tools')
work=Path('/tmp/atst-doc-onboarding-20261015/tutorial')
for case in ['09_lightweight_cli','06_relax_H2-Au']:
 shutil.copytree(root/'examples'/case, work/'examples'/case, dirs_exist_ok=True)
env=dict(os.environ, PYTHONPATH=os.environ.get('ATST_DOC_SOURCE', str(root/'src')))
print('PACKAGE_SOURCE:', env['PYTHONPATH'])
cli=[sys.executable,'-m','atst_tools.scripts.cli']
commands=[(work,['--version']),(work,['run','--list-types']),(work,['config','validate','examples/06_relax_H2-Au/config.yaml','--print-normalized']),(work,['run','--dry-run','examples/06_relax_H2-Au/config.yaml']),(work/'examples/09_lightweight_cli',['neb','make','inputs/init.xyz','inputs/final.xyz','3','-o','first_chain.traj','--method','linear']),(work/'examples/09_lightweight_cli',['neb','summary','inputs/neb_result.extxyz','--n-max','1','--tail','5']),(work/'examples/09_lightweight_cli',['relax','post','inputs/relax_result.extxyz','--output-format','traj','--output','restart.traj'])]
for cwd,args in commands:
 result=subprocess.run(cli+args,cwd=cwd,env=env,text=True,capture_output=True)
 print('COMMAND:', 'atst '+' '.join(args), '\nCWD:',cwd,'\nEXIT:', result.returncode, '\n',result.stdout,result.stderr)
 assert result.returncode==0
chain=read(work/'examples/09_lightweight_cli/first_chain.traj',index=':')
assert len(chain)==5
assert (work/'examples/09_lightweight_cli/restart.traj').is_file()
print('ONBOARDING_OK: 7 commands; 5-frame chain and restart artifact; no backend calculation')
