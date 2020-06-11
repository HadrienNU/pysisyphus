import io
import subprocess
import tempfile

from jinja2 import Template

from pysisyphus.config import get_cmd


TPL = """tolerance {{ tolerance }}
output {{ output_fn }}
filetype pdb

# Solute
structure {{ solute_fn }}
 number {{ solute_num }}
 {% if solute_num == 1 %}
     fixed 0. 0. 0. 0. 0. 0.
     centerofmass
 {% endif %}
end structure

# Solvent
structure {{ solvent_fn }}
 number {{ solvent_num }}
 inside sphere 0. 0. 0. {{ sphere_radius }}
end structure
"""


def make_input(output_fn, solute_fn, solute_num, solvent_fn, solvent_num,
               sphere_radius, tolerance=2.0):
    assert solute_num == 1, "Handling of > solutes is not yet implemented"

    tpl = Template(TPL)
    rendered = tpl.render(
                tolerance=tolerance,
                output_fn=output_fn,
                solute_fn=solute_fn,
                solute_num=solute_num,
                solvent_fn=solvent_fn,
                solvent_num=solvent_num,
                sphere_radius=sphere_radius,
    )
    return rendered


def call_packmol(inp):
    packmol_cmd = get_cmd("packmol")

    with tempfile.NamedTemporaryFile(mode="w", dir=".") as tmp:
        tmp.write(inp)
        cmd = f"{packmol_cmd} < {tmp.name}"
        proc = subprocess.run(cmd, shell=True,
                              stdout=subprocess.PIPE, #stderr=subprocess.PIPE,
                              text=True)
    out = proc.stdout
    # err = proc.stderr

    return proc
