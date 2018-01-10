#!/usr/bin/env python3

import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

def run_ext(out_fn, args, path, env):
    #print("out_fn", out_fn)
    #print("args", args)
    #print("path", path)
    #print("env", env)
    with open(path / out_fn, "w") as handle:
        result = subprocess.Popen(args, cwd=path, stdout=handle, env=env)
        result.wait()
    #return result


class Calculator:
    logger = logging.getLogger("calculator")

    def __init__(self, charge=0, mult=1, name="calculator", client=None):
        self.charge = int(charge)
        self.mult = int(mult)
        self.name = name
        self.client = client

        self.counter = 0
        self._energy = None
        self._forces = None
        self._hessian = None

        self.inp_fn = "calc.inp"
        self.out_fn = "calc.out"

    def log(self, message):
        self.logger.debug(f"{self.name}_{self.counter:03d}, " + message)

    def get_energy(self, atoms, coords):
        raise Exception("Not implemented!")

    def get_hessian(self, atoms, coords):
        raise Exception("Not implemented!")

    def make_fn(self, ext):
        return f"{self.name}.{self.counter:03d}.{ext}"

    def prepare(self, inp, path=None):
        if not path:
            prefix = f"{self.name}_{self.counter:03d}_"
            path = Path(tempfile.mkdtemp(prefix=prefix))
        inp_path = path / self.inp_fn
        with open(inp_path, "w") as handle:
            handle.write(inp)

        return path

    def run(self, inp, calc, add_args=None, env=None):
        path = self.prepare(inp)
        self.log(f"running in {path}")
        args = [self.base_cmd, self.inp_fn]
        if add_args:
            args.extend(add_args)
        if not env:
            env = os.environ.copy()
        # Use dask for parallel computing
        if self.client:
        #if False:
            #print("using dask")
            #print(self.client)
            #import pdb; pdb.set_trace()
            task = self.client.submit(run_ext, self.out_fn, args, path, env)
            #print(self.client)
            self.client.gather(task)
            #run_ext(self.out_fn, args, path, env)
        else:
            run_ext(self.out_fn, args, path, env)
        try:
            results = self.parser_funcs[calc](path)
            self.keep(path)
        except Exception as err:
            print(err)
            print()
            print("Crashed input:")
            print(inp)
            backup_dir = Path(os.getcwd()) / f"crashed_{self.name}"
            if backup_dir.exists():
                shutil.rmtree(backup_dir)
            shutil.copytree(path, backup_dir)
            sys.exit()
        finally:
            #self.clean(path)
            self.counter += 1

        return results

    def keep(self, path, exts=()):
        kept_fns = dict()
        for ext in exts:
            pattern = f"*.{ext}"
            globbed = list(path.glob(pattern))
            assert(len(globbed) <= 1), (f"Expected at most one {pattern} in {path}."
                                        f" Found {len(globbed)} instead!"
            )
            if len(globbed) == 0:
                continue
            old_fn = globbed[0]
            new_fn = os.path.abspath(self.make_fn(ext))
            shutil.copy(old_fn, new_fn)
            kept_fns[ext] = new_fn
        return kept_fns

    def clean(self, path):
        shutil.rmtree(path)
        self.log(f"cleaned {path}")


if __name__ == "__main__":
    calc = Calculator()
    calc.base_cmd = "sleep"
    calc.run("dummy_inp", "dummy_calc_type")