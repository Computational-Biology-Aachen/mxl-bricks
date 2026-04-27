# ruff: noqa: INP001
"""CST-based transforms for normalising mxlbricks enzyme-function signatures."""

from __future__ import annotations

import re
from pathlib import Path  # noqa: TC003

import libcst as cst
import typer

app = typer.Typer()

_FLOAT = cst.Annotation(annotation=cst.Name("float"))
_COMMA = cst.Comma(whitespace_after=cst.SimpleWhitespace(""))
_MIN_CALL_ARGS = 2


def _with_trailing_comma_param(p: cst.Param) -> cst.Param:
    return (
        p
        if not isinstance(p.comma, cst.MaybeSentinel)
        else p.with_changes(comma=_COMMA)
    )


def _with_trailing_comma_arg(a: cst.Arg) -> cst.Arg:
    return (
        a
        if not isinstance(a.comma, cst.MaybeSentinel)
        else a.with_changes(comma=_COMMA)
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _to_snake(name: str) -> str:
    s = re.sub(r"([A-Z]+)([A-Z][a-z])", r"\1_\2", name)
    s = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", s)
    return s.lower()


def _normalize_fn_name(name: str) -> str:
    snake = _to_snake(name.lstrip("_"))
    return f"_{snake}"


def _returns_model(node: cst.FunctionDef) -> bool:
    if node.returns is None:
        return False
    ann = node.returns.annotation
    return isinstance(ann, cst.Name) and ann.value == "Model"


def _apply(path: Path, transformer: cst.CSTTransformer) -> bool:
    source = path.read_text()
    new_code = cst.parse_module(source).visit(transformer).code
    if new_code == source:
        return False
    path.write_text(new_code)
    return True


# ---------------------------------------------------------------------------
# Transformers
# ---------------------------------------------------------------------------


class AddFloatAnnotations(cst.CSTTransformer):
    """CST transformer that adds `float` type annotations to unannotated params and return types."""

    def leave_Param(  # noqa: N802
        self,
        original_node: cst.Param,  # noqa: ARG002
        updated_node: cst.Param,
    ) -> cst.Param:
        """Add float annotation to unannotated non-star parameters."""
        if updated_node.annotation is not None or updated_node.star in ("*", "**"):
            return updated_node
        return updated_node.with_changes(annotation=_FLOAT)

    def leave_FunctionDef(  # noqa: N802
        self,
        original_node: cst.FunctionDef,  # noqa: ARG002
        updated_node: cst.FunctionDef,
    ) -> cst.FunctionDef:
        """Add float return annotation to functions missing one."""
        if updated_node.returns is not None:
            return updated_node
        return updated_node.with_changes(returns=_FLOAT)


class AddTrailingCommas(cst.CSTTransformer):
    """CST transformer that adds trailing commas to function parameters and call arguments."""

    def leave_Parameters(  # noqa: N802
        self,
        original_node: cst.Parameters,  # noqa: ARG002
        updated_node: cst.Parameters,
    ) -> cst.Parameters:
        """Add trailing comma to the last parameter in a function signature."""
        p = updated_node
        if p.star_kwarg is not None:
            return p.with_changes(star_kwarg=_with_trailing_comma_param(p.star_kwarg))
        if p.kwonly_params:
            last = _with_trailing_comma_param(p.kwonly_params[-1])
            return p.with_changes(kwonly_params=(*p.kwonly_params[:-1], last))
        if isinstance(p.star_arg, cst.Param):
            return p.with_changes(star_arg=_with_trailing_comma_param(p.star_arg))
        if p.params:
            last = _with_trailing_comma_param(p.params[-1])
            return p.with_changes(params=(*p.params[:-1], last))
        if p.posonly_params:
            last = _with_trailing_comma_param(p.posonly_params[-1])
            return p.with_changes(posonly_params=(*p.posonly_params[:-1], last))
        return p

    def leave_Call(  # noqa: N802
        self,
        original_node: cst.Call,  # noqa: ARG002
        updated_node: cst.Call,
    ) -> cst.Call:
        """Add trailing comma to the last argument in a function call."""
        if len(updated_node.args) < _MIN_CALL_ARGS:
            return updated_node
        last = _with_trailing_comma_arg(updated_node.args[-1])
        return updated_node.with_changes(args=(*updated_node.args[:-1], last))


class _CollectNonModelFns(cst.CSTVisitor):
    def __init__(self) -> None:
        self.renames: dict[str, str] = {}

    def visit_FunctionDef(self, node: cst.FunctionDef) -> None:  # noqa: N802
        if _returns_model(node):
            return
        old = node.name.value
        new = _normalize_fn_name(old)
        if old != new:
            self.renames[old] = new


class _ApplyRenames(cst.CSTTransformer):
    def __init__(self, renames: dict[str, str]) -> None:
        self.renames = renames

    def leave_Name(  # noqa: N802
        self,
        original_node: cst.Name,  # noqa: ARG002
        updated_node: cst.Name,
    ) -> cst.Name:
        new = self.renames.get(updated_node.value)
        return updated_node.with_changes(value=new) if new else updated_node

    def leave_Attribute(  # noqa: N802
        self, original_node: cst.Attribute, updated_node: cst.Attribute
    ) -> cst.Attribute:
        # leave_Name would otherwise rename the attr part of dotted access
        return updated_node.with_changes(attr=original_node.attr)


class _CollectPrivateFnParamRenames(cst.CSTVisitor):
    def __init__(self) -> None:
        self.renames: dict[str, str] = {}

    def visit_FunctionDef(self, node: cst.FunctionDef) -> None:  # noqa: N802
        if not node.name.value.startswith("_"):
            return
        all_params = [
            *node.params.params,
            *node.params.posonly_params,
            *node.params.kwonly_params,
        ]
        if isinstance(node.params.star_arg, cst.Param):
            all_params.append(node.params.star_arg)
        if node.params.star_kwarg is not None:
            all_params.append(node.params.star_kwarg)
        for param in all_params:
            old = param.name.value
            new = _to_snake(old)
            if old != new:
                self.renames[old] = new


class _ApplyParamRenames(cst.CSTTransformer):
    """Rename params and their usages, scoped to _-prefixed functions only."""

    def __init__(self, renames: dict[str, str]) -> None:
        self.renames = renames
        self._in_private: list[bool] = []

    def visit_FunctionDef(self, node: cst.FunctionDef) -> bool:  # noqa: N802
        self._in_private.append(node.name.value.startswith("_"))
        return True

    def leave_FunctionDef(  # noqa: N802
        self,
        original_node: cst.FunctionDef,  # noqa: ARG002
        updated_node: cst.FunctionDef,
    ) -> cst.FunctionDef:
        self._in_private.pop()
        return updated_node

    def leave_Param(  # noqa: N802
        self,
        original_node: cst.Param,  # noqa: ARG002
        updated_node: cst.Param,
    ) -> cst.Param:
        if not self._in_private or not self._in_private[-1]:
            return updated_node
        new = self.renames.get(updated_node.name.value)
        return updated_node.with_changes(name=cst.Name(new)) if new else updated_node

    def leave_Name(  # noqa: N802
        self,
        original_node: cst.Name,  # noqa: ARG002
        updated_node: cst.Name,
    ) -> cst.Name:
        if not self._in_private or not self._in_private[-1]:
            return updated_node
        new = self.renames.get(updated_node.value)
        return updated_node.with_changes(value=new) if new else updated_node

    def leave_Attribute(  # noqa: N802
        self,
        original_node: cst.Attribute,
        updated_node: cst.Attribute,
    ) -> cst.Attribute:
        return updated_node.with_changes(attr=original_node.attr)


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------


@app.command()
def annotate(
    files: list[Path] = typer.Argument(..., help="Python files to annotate"),  # noqa: B008
) -> None:
    """Add float annotations to all unannotated function params and return types."""
    for path in files:
        changed = _apply(path, AddFloatAnnotations())
        typer.echo(f"{'updated' if changed else 'no changes'} {path}")


@app.command()
def normalize_names(
    files: list[Path] = typer.Argument(..., help="Python files to normalise"),  # noqa: B008
) -> None:
    """Prefix non-Model-returning fns with _ and convert names to snake_case."""
    for path in files:
        source = path.read_text()
        tree = cst.parse_module(source)

        collector = _CollectNonModelFns()
        tree.visit(collector)

        if not collector.renames:
            typer.echo(f"no changes {path}")
            continue

        for old, new in collector.renames.items():
            typer.echo(f"  {old} → {new}")

        new_code = tree.visit(_ApplyRenames(collector.renames)).code
        if new_code != source:
            path.write_text(new_code)
            typer.echo(f"updated {path}")


@app.command()
def normalize_fn_pars(
    files: list[Path] = typer.Argument(..., help="Python files to normalise"),  # noqa: B008
) -> None:
    """Convert camelCase/PascalCase param names to snake_case in _-prefixed functions."""
    for path in files:
        source = path.read_text()
        tree = cst.parse_module(source)

        collector = _CollectPrivateFnParamRenames()
        tree.visit(collector)

        if not collector.renames:
            typer.echo(f"no changes {path}")
            continue

        for old, new in collector.renames.items():
            typer.echo(f"  {old} → {new}")

        new_code = tree.visit(_ApplyParamRenames(collector.renames)).code
        if new_code != source:
            path.write_text(new_code)
            typer.echo(f"updated {path}")


@app.command()
def trailing_commas(
    files: list[Path] = typer.Argument(..., help="Python files to fix"),  # noqa: B008
) -> None:
    """Add trailing comma after last param/arg in all function defs and call sites."""
    for path in files:
        changed = _apply(path, AddTrailingCommas())
        typer.echo(f"{'updated' if changed else 'no changes'} {path}")


if __name__ == "__main__":
    app()
