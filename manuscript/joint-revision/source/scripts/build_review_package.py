"""Build the reconciled manuscript and retain logs and recorder provenance."""
from pathlib import Path
import shutil
import subprocess

root = Path(__file__).resolve().parents[2]
build = root / "build"
build.mkdir(exist_ok=True)
def tex_program(name):
    found = shutil.which(name)
    if found:
        return found
    mac_path = Path("/Library/TeX/texbin") / name
    if mac_path.is_file():
        return str(mac_path)
    raise SystemExit(f"Required TeX program not found: {name}; install TeX Live and add its bin directory to PATH")

pdflatex = tex_program("pdflatex")
bibtex = tex_program("bibtex")
documents = [
    ("main_tables_figures", "tables_refs"),
    ("main_critical_final", "main_refs"),
    ("appendix", "appendix_refs"),
    ("main_critical_final", "main_refs"),
    ("appendix", "appendix_refs"),
    ("main_critical_final", "main_refs"),
]
for cycle, (name, refs) in enumerate(documents, 1):
    tex = [pdflatex, "-recorder", "-interaction=nonstopmode",
           "-halt-on-error", "-output-directory=build", name + ".tex"]
    commands = [tex]
    if name != "main_tables_figures":
        commands.append([bibtex, "build/" + name])
    commands += [tex] * 3
    for pass_no, command in enumerate(commands, 1):
        result = subprocess.run(command, cwd=root, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (build / f"{cycle:02d}-{name}-{pass_no}.txt").write_text(result.stdout)
        if result.returncode:
            print(result.stdout[-6000:])
            raise SystemExit(f"Failed: {name}, cycle {cycle}, pass {pass_no}")
    labels = [line for line in (build / (name + ".aux")).read_text().splitlines()
              if line.startswith(r"\newlabel")]
    (root / "source/generated" / (refs + ".aux")).write_text("\n".join(labels) + "\n")
    shutil.copy2(build / (name + ".pdf"), root / (name + ".pdf"))
    if name == "main_critical_final":
        shutil.copy2(build / (name + ".pdf"), root / "main_draft.pdf")
        shutil.copy2(build / (name + ".pdf"), root / "main.pdf")
    print(f"Compiled {name}.pdf (cycle {cycle})", flush=True)

issues = []
for name in ("main_critical_final", "main_tables_figures", "appendix"):
    for line in (build / (name + ".log")).read_text(errors="replace").splitlines():
        if any(word in line.lower() for word in
               ("undefined", "overfull", "multiply defined", "duplicate ignored", "fatal")):
            issues.append(f"{name}: {line}")
(build / "FINAL-WARNINGS.txt").write_text("\n".join(issues) + ("\n" if issues else ""))
print("Warnings: " + str(len(issues)))
for line in issues:
    print(line)
