"""Check ATST-Tools documentation invariants."""

from __future__ import annotations

import argparse
import re
import sys
from html.parser import HTMLParser
from pathlib import Path


ACTIVE_MARKDOWN_ROOTS = (
    Path("README.md"),
    Path("src/atst_tools/README.md"),
    Path("examples/README.md"),
    Path("AGENTS.md"),
    Path("docs"),
)
REQUIRED_ENTRY_DOCUMENTS = (
    Path("README.md"),
    Path("docs/index.md"),
    Path("docs/developer/DOCUMENTATION.md"),
    Path("src/atst_tools/README.md"),
    Path("examples/README.md"),
)
CONFLICT_MARKER = re.compile(r"^(<<<<<<<|=======|>>>>>>>)", re.MULTILINE)
MARKDOWN_LINK = re.compile(r"(?<!!)\[[^\]]+\]\(([^)]+)\)")


def active_markdown_files(root: Path) -> list[Path]:
    """Return active Markdown files that participate in link checks."""
    files: list[Path] = []
    for entry in ACTIVE_MARKDOWN_ROOTS:
        path = root / entry
        if path.is_file() and path.suffix == ".md":
            files.append(path)
        elif path.is_dir():
            for markdown in path.rglob("*.md"):
                relative = markdown.relative_to(root)
                if len(relative.parts) >= 2 and relative.parts[:2] == ("docs", "archive"):
                    continue
                files.append(markdown)
    return sorted(set(files))


def _strip_anchor(target: str) -> str:
    target = target.strip()
    if target.startswith("<") and target.endswith(">"):
        target = target[1:-1]
    return target.split("#", 1)[0]


def check_conflict_markers(root: Path) -> list[str]:
    """Check active Markdown files for unresolved conflict markers."""
    issues: list[str] = []
    for path in active_markdown_files(root):
        text = path.read_text(encoding="utf-8")
        if CONFLICT_MARKER.search(text):
            issues.append(f"{path.relative_to(root)} contains merge conflict marker")
    return issues


def check_markdown_links(root: Path) -> list[str]:
    """Check local relative Markdown links in active Markdown files."""
    issues: list[str] = []
    for path in active_markdown_files(root):
        text = path.read_text(encoding="utf-8")
        for match in MARKDOWN_LINK.finditer(text):
            target = _strip_anchor(match.group(1))
            if not target or target.startswith(("http://", "https://", "mailto:")):
                continue
            if "://" in target:
                continue
            resolved = (path.parent / target).resolve()
            try:
                resolved.relative_to(root)
            except ValueError:
                issues.append(f"{path.relative_to(root)} links outside repository: {target}")
                continue
            if not resolved.exists():
                issues.append(f"{path.relative_to(root)} has missing link target: {target}")
    return issues


def check_entry_documents(root: Path) -> list[str]:
    """Require the stable navigation entry points to exist."""
    return [
        f"required entry document is missing: {path.as_posix()}"
        for path in REQUIRED_ENTRY_DOCUMENTS
        if not (root / path).is_file()
    ]


def check_html_reports(root: Path) -> list[str]:
    """Parse active HTML reports with Python's standard HTML parser."""
    issues: list[str] = []
    for path in sorted((root / "docs" / "reports").glob("*.html")):
        parser = HTMLParser()
        try:
            parser.feed(path.read_text(encoding="utf-8"))
            parser.close()
        except Exception as exc:  # pragma: no cover - HTMLParser rarely raises.
            issues.append(f"{path.relative_to(root)} failed HTMLParser parse: {exc}")
    return issues


def check_repository(root: Path) -> list[str]:
    """Return documentation governance issues for a repository root."""
    root = root.resolve()
    issues: list[str] = []
    issues.extend(check_conflict_markers(root))
    issues.extend(check_markdown_links(root))
    issues.extend(check_entry_documents(root))
    issues.extend(check_html_reports(root))
    return issues


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="Repository root to check")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    args = parse_args(argv)
    issues = check_repository(args.root)
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}", file=sys.stderr)
        return 1
    print("documentation governance checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
