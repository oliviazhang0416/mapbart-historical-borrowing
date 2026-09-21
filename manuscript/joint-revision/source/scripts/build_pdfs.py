"""Compile the three standalone documents without retaining temporary build files."""
from pathlib import Path
import shutil, subprocess, tempfile

source=Path(__file__).resolve().parents[1]
manuscript=source.parent
texbin=Path('/Library/TeX/texbin')
with tempfile.TemporaryDirectory(prefix='build-',dir=source) as temporary:
    build=Path(temporary)
    relative_build=build.relative_to(manuscript)
    for name in ['main_tables_figures','main','appendix','main']:
        commands=[[str(texbin/'pdflatex'),'-interaction=nonstopmode','-halt-on-error',
                   f'-output-directory={relative_build}',name+'.tex']]
        if name!='main_tables_figures':
            commands.append([str(texbin/'bibtex'),str(relative_build/name)])
        commands+= [commands[0]]*3
        for i,cmd in enumerate(commands):
            result=subprocess.run(cmd,cwd=manuscript,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            if result.returncode:
                print(result.stdout[-5000:])
                raise SystemExit(f'Compilation failed: {name}, pass {i+1}')
        log=(build/(name+'.log')).read_text(errors='replace')
        warnings=[line for line in log.splitlines() if 'Overfull' in line or 'undefined' in line.lower()]
        if warnings:
            print(name+': '+ '\n'.join(warnings))
        if name in ['main_tables_figures','appendix','main']:
            refname = {'main_tables_figures':'tables_refs.aux', 'appendix':'appendix_refs.aux', 'main':'main_refs.aux'}[name]
            labels = [line for line in (build/(name+'.aux')).read_text().splitlines() if line.startswith(r'\newlabel')]
            (source/'generated'/refname).write_text('\n'.join(labels)+'\n')
        shutil.copy2(build/(name+'.pdf'),manuscript/(name+'.pdf'))
        print(f'Compiled {name}.pdf')
