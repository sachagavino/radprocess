#!/usr/bin/env python3
"""
Static API check of a script that drives the radprocess Pipeline.

Reads the *current* radprocess source (Pipeline.py and utils/config.py) with the
``ast`` module -- nothing is imported, so neither radprocess's dependencies nor
POLARIS / RADMC-3D are needed -- and reports, for the script(s) given:

  * calls to Pipeline methods with keyword arguments the method does not accept,
    or with too many positional arguments;
  * calls to Pipeline methods that do not exist;
  * reads / assignments of ``configparams.<group>.<field>`` where the group or the
    field does not exist (the config objects are strict: this raises at run time);
  * names that were removed from the API, with their replacement.

Usage
-----
    python tools/check_pipeline_api.py radprocess_pipeline.py [more.py ...]
    python tools/check_pipeline_api.py --src /path/to/radprocess script.py
    python tools/check_pipeline_api.py --list          # print the current API

Exit status is 1 if any error was found, so it can run in CI before a long job.

Limits: the check is static. It recognises Pipeline method calls by method name
(``anything.run_polaris_opacity(...)``) and config access by attribute chain
(``....dust.mixtures``, ``cfg.imaging.npix``); it does not validate *values*
(types, units, the content of ``dust.mixtures``) -- see the Parameters reference
in the documentation for those.
"""

import argparse
import ast
import importlib.util
import sys
from pathlib import Path


# Names removed from the public API, with what replaces them.
REMOVED = {
    # run_polaris_opacity / render_images arguments
    "dust_components": "use configparams.dust.mixtures (or the dust_mixtures= argument); "
                       "see the 'Dust mixtures' section of the Parameters reference",
    "dust_size_min": "grain sizes now come from each component's 'amin' in "
                     "configparams.dust.mixtures (or configparams.dust.bin_edges)",
    "dust_size_max": "grain sizes now come from each component's 'amax' in "
                     "configparams.dust.mixtures (or configparams.dust.bin_edges)",
    "dust_size_powerlaw": "the size-distribution exponent is now each component's "
                          "'index' in configparams.dust.mixtures",
    "midplane_zoom": "removed; use fov_au= to zoom",
    "fov_m": "removed; use fov_au= (field of view in au, full width)",
    # config fields
    "dust": "amrsource.dust does not exist; use amrsource.dustratios "
            "(or amrsource.fluiddens for multi-fluid outputs)",
}


# ---------------------------------------------------------------------------
#  Read the current API from the source
# ---------------------------------------------------------------------------

def find_src(explicit=None):
    """Directory of the radprocess package (the one containing pipeline/)."""
    if explicit:
        p = Path(explicit)
        return p / "radprocess" if (p / "radprocess" / "pipeline").is_dir() else p
    here = Path(__file__).resolve().parent.parent / "radprocess"
    if (here / "pipeline" / "Pipeline.py").exists():
        return here
    spec = importlib.util.find_spec("radprocess")  # does not execute the package
    if spec and spec.submodule_search_locations:
        return Path(list(spec.submodule_search_locations)[0])
    sys.exit("Cannot locate the radprocess source; pass --src.")


def pipeline_methods(src):
    """{method_name: (positional_names, keyword_names, has_varargs, has_varkw)}"""
    tree = ast.parse((src / "pipeline" / "Pipeline.py").read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "Pipeline")
    methods = {}
    for node in cls.body:
        if not isinstance(node, ast.FunctionDef) or node.name.startswith("_"):
            continue
        if any(isinstance(d, ast.Name) and d.id == "property" for d in node.decorator_list):
            methods[node.name] = None                                   # property
            continue
        a = node.args
        pos = [x.arg for x in a.posonlyargs + a.args][1:]              # drop self
        kw = pos + [x.arg for x in a.kwonlyargs]
        methods[node.name] = (pos, kw, a.vararg is not None, a.kwarg is not None)
    return methods


def config_schema(src):
    """{group_name: set(field_names)} for ConfigParams, plus top-level fields."""
    tree = ast.parse((src / "utils" / "config.py").read_text())
    classes, class_methods = {}, {}
    for node in tree.body:
        if isinstance(node, ast.ClassDef):
            classes[node.name] = [
                (s.target.id, ast.unparse(s.annotation))
                for s in node.body
                if isinstance(s, ast.AnnAssign) and isinstance(s.target, ast.Name)
            ]
            class_methods[node.name] = {
                s.name for s in node.body
                if isinstance(s, ast.FunctionDef) and not s.name.startswith("_")
            }
    groups, top = {}, set(class_methods["ConfigParams"])     # e.g. validate()
    for name, ann in classes["ConfigParams"]:
        top.add(name)
        if ann in classes:
            groups[name] = {f for f, _ in classes[ann]} | class_methods[ann]
    return groups, top


# ---------------------------------------------------------------------------
#  Check a script
# ---------------------------------------------------------------------------

def attr_chain(node):
    """['pipe', 'configparams', 'dust', 'mixtures'] for pipe.configparams.dust.mixtures"""
    parts = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        node = node.value
    if isinstance(node, ast.Name):
        parts.append(node.id)
    return parts[::-1]


def check_file(path, methods, groups, top):
    tree = ast.parse(Path(path).read_text(), filename=str(path))
    issues = []

    def report(node, msg):
        issues.append((node.lineno, msg))

    # Names bound to pipe.configparams (e.g. cfg = pipe.configparams)
    cfg_aliases = {"configparams"}
    for node in ast.walk(tree):
        if (isinstance(node, ast.Assign) and isinstance(node.value, ast.Attribute)
                and node.value.attr == "configparams"):
            for t in node.targets:
                if isinstance(t, ast.Name):
                    cfg_aliases.add(t.id)

    seen_attr = set()
    for node in ast.walk(tree):
        # ---- Pipeline method calls ----
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            name = node.func.attr
            chain = attr_chain(node.func.value)
            if name in methods and methods[name] is not None and "configparams" not in chain:
                pos, kw, varargs, varkw = methods[name]
                if len(node.args) > len(pos) and not varargs:
                    report(node, f"{name}() takes {len(pos)} positional argument(s), "
                                 f"{len(node.args)} given")
                for k in node.keywords:
                    if k.arg is None or k.arg in kw or varkw:
                        continue
                    hint = REMOVED.get(k.arg)
                    msg = f"{name}() has no argument '{k.arg}'"
                    if hint:
                        msg += f" -> {hint}"
                    report(node, msg)

        # ---- configparams.<group>.<field> ----
        if isinstance(node, ast.Attribute) and id(node) not in seen_attr:
            chain = attr_chain(node)
            for i, part in enumerate(chain[:-1]):
                if part in cfg_aliases:
                    rest = chain[i + 1:]
                    # mark inner attribute nodes so the same chain is reported once
                    inner = node.value
                    while isinstance(inner, ast.Attribute):
                        seen_attr.add(id(inner))
                        inner = inner.value
                    if rest[0] not in top:
                        report(node, f"configparams has no group '{rest[0]}' "
                                     f"(valid: {', '.join(sorted(top))})")
                    elif rest[0] in groups and len(rest) > 1 and rest[1] not in groups[rest[0]]:
                        msg = f"configparams.{rest[0]} has no field '{rest[1]}'"
                        hint = REMOVED.get(rest[1])
                        if hint:
                            msg += f" -> {hint}"
                        else:
                            msg += f" (valid: {', '.join(sorted(groups[rest[0]]))})"
                        report(node, msg)
                    break

    return sorted(set(issues))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("scripts", nargs="*", help="Python files to check")
    ap.add_argument("--src", help="radprocess repository or package directory")
    ap.add_argument("--list", action="store_true", help="print the current API and exit")
    args = ap.parse_args(argv)

    src = find_src(args.src)
    methods = pipeline_methods(src)
    groups, top = config_schema(src)

    if args.list:
        print(f"radprocess source: {src}\n\nPipeline methods:")
        for m, sig in sorted(methods.items()):
            if sig is not None:
                print(f"  {m}({', '.join(sig[1])})")
        print("\nconfigparams:")
        for g in sorted(top):
            print(f"  {g}: {', '.join(sorted(groups[g]))}" if g in groups else f"  {g}")
        return 0

    n_err = 0
    for script in args.scripts:
        issues = check_file(script, methods, groups, top)
        n_err += len(issues)
        status = "OK" if not issues else f"{len(issues)} problem(s)"
        print(f"{script}: {status}")
        for line, msg in issues:
            print(f"  line {line}: {msg}")
    return 1 if n_err else 0


if __name__ == "__main__":
    sys.exit(main())
