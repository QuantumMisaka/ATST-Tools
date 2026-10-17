"""Summarize syntactic first-party imports without importing scientific code."""
import argparse
import ast
import collections
import hashlib
import json
import pathlib
import subprocess

root = pathlib.Path.cwd()
package = root / 'src/atst_tools'
files = sorted(p for p in package.rglob('*.py') if 'external' not in p.relative_to(package).parts and '__pycache__' not in p.parts)
modules = {}
for path in files:
    parts = path.relative_to(root / 'src').with_suffix('').parts
    name = '.'.join(parts[:-1] if parts[-1] == '__init__' else parts)
    modules[name] = path
edges = []
dynamic = []
for source, path in modules.items():
    base = source.split('.') if path.name == '__init__.py' else source.split('.')[:-1]
    class Visitor(ast.NodeVisitor):
        depth = 0
        def visit_FunctionDef(self, node):
            self.depth += 1
            self.generic_visit(node)
            self.depth -= 1
        visit_AsyncFunctionDef = visit_FunctionDef
        def record(self, target, node):
            if target in modules and target != source:
                edges.append({'source': source, 'target': target, 'line': node.lineno, 'inside_function': bool(self.depth)})
        def visit_Import(self, node):
            for alias in node.names:
                self.record(alias.name, node)
        def visit_ImportFrom(self, node):
            head = '.'.join(base[:len(base) - node.level + 1]) if node.level else ''
            target = '.'.join(p for p in (head, node.module or '') if p)
            self.record(target, node)
            for alias in node.names:
                self.record(target + '.' + alias.name, node)
        def visit_Call(self, node):
            name = ast.unparse(node.func)
            if name.endswith('import_module') or name == '__import__':
                dynamic.append({'source': source, 'line': node.lineno, 'expression': ast.unparse(node)})
            self.generic_visit(node)
    Visitor().visit(ast.parse(path.read_text()))

def group(name):
    return name.split('.')[1] if '.' in name else '(root)'

def scc(selected):
    graph = {name: set() for name in modules}
    for edge in selected:
        graph[edge['source']].add(edge['target'])
    index, stack, active, indices, low, output = 0, [], set(), {}, {}, []
    def visit(v):
        nonlocal index
        indices[v] = low[v] = index
        index += 1
        stack.append(v)
        active.add(v)
        for w in sorted(graph[v]):
            if w not in indices:
                visit(w)
                low[v] = min(low[v], low[w])
            elif w in active:
                low[v] = min(low[v], indices[w])
        if low[v] == indices[v]:
            component = []
            while True:
                w = stack.pop()
                active.remove(w)
                component.append(w)
                if w == v:
                    break
            if len(component) > 1:
                output.append(sorted(component))
    for v in sorted(graph):
        if v not in indices:
            visit(v)
    return output

summary = {
    'head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'method': 'AST explicit imports; excludes vendored external and dynamic imports; non-function means syntactic module/class/conditional body, not guaranteed runtime execution; SCC is syntactic, not proof of an import failure.',
    'files': [{'path': str(p.relative_to(root)), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in files],
    'module_count': len(modules),
    'module_counts_by_directory': dict(sorted(collections.Counter(group(m) for m in modules).items())),
    'cross_directory_edges': {g: sorted({group(e['target']) for e in edges if group(e['source']) == g and group(e['target']) != g}) for g in sorted({group(m) for m in modules})},
    'scc_all_explicit_imports': scc(edges),
    'scc_excluding_function_bodies': scc([e for e in edges if not e['inside_function']]),
    'dynamic_imports': dynamic,
    'edges': edges,
}
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=pathlib.Path, help='Optional JSON evidence path')
arguments = parser.parse_args()
if arguments.output is not None:
    arguments.output.write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps({k: v for k, v in summary.items() if k not in {'files', 'edges'}}, indent=2))
