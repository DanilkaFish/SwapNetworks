from dataclasses import dataclass, field
from typing import Any, List, Optional, Sequence, Union
import pathlib
import re

_LATEX_SPECIALS = {
    '&': r'\&', '%': r'\%', '$': r'\$', '#': r'\#', '_': r'\_',
    '{': r'\{', '}': r'\}', '~': r'\textasciitilde{}',
    '^': r'\textasciicircum{}', '\\': r'\textbackslash{}',
}


def escape_latex(s: Any) -> str:
    if s is None:
        return ""
    s = str(s)
    if s.startswith("$") and s.endswith("$") and len(s) > 1:
        return s
    return "".join(s)

    # return "".join(_LATEX_SPECIALS.get(ch, ch) for ch in s)


def centered_m(width: str) -> str:
    return rf">{{\centering\arraybackslash}}m{{{width}}}"


_unit_re = re.compile(r"^\d+(\.\d+)?(em|cm|pt|mm)$")


def _normalize_col_spec(spec: Union[str, Sequence[Union[str, float, int]]], ncols: int) -> List[str]:

    spec_list = list(spec)
    if len(spec_list) != ncols:
        raise ValueError(f"Column-spec sequence length {len(spec_list)} must equal number of columns {ncols}.")

    normalized: List[str] = []
    for item in spec_list:
        if isinstance(item, (int, float)):
            normalized.append(centered_m(f"{item}em"))
            continue

        if isinstance(item, str):
            s = item.strip()
            if _unit_re.match(s):  # looks like "5em" or "2.5cm"
                normalized.append(centered_m(s))
            else:
                normalized.append(s)
            continue
        raise TypeError("Each column spec must be a str, int, or float.")

    return normalized


@dataclass
class Table:
    caption: str
    column_names: List[str]
    align: List[float] = "c"
    bordered: bool = False
    label: Optional[str] = None
    placement: str = "h"
    use_booktabs: bool = True
    use_adjustbox: bool = True
    rows: List[Sequence[Any]] = field(default_factory=list)
    float_format: Optional[str] = None

    def __post_init__(self):
        self.ncols = len(self.column_names)
        self.align_list = _normalize_col_spec(self.align, self.ncols)

    def add_row(self, *values: Any):
        if len(values) != self.ncols:
            raise ValueError("Number of values must equal number of columns.")
        self.rows.append(values)

    def _format_cell(self, val: Any) -> str:
        if val is None:
            return ""
        if self.float_format and isinstance(val, float):
            try:
                return self.float_format.format(val)
            except Exception:
                pass
        return escape_latex(val)

    def _column_spec_latex(self) -> str:
        col_spec = "".join(self.align_list)
        if self.bordered:
            return "|" + "|".join(self.align_list) + "|"
        return col_spec

    def header(self) -> List[str]:
        parts = [f"\\begin{{table}}[{self.placement}]", "\\centering"]
        if self.caption:
            parts.append(f"\\caption{{{escape_latex(self.caption)}}}")
        if self.label:
            parts.append(f"\\label{{{escape_latex(self.label)}}}")
        if self.use_adjustbox:
            parts.append("\\begin{adjustbox}{max width=\\linewidth}")
        parts.append(f"\\begin{{tabular}}{{{self._column_spec_latex()}}}")
        if self.use_booktabs:
            parts.append("\\toprule")
        parts.append(" & ".join(f"\\textbf{{{escape_latex(c)}}}" for c in self.column_names) + " \\\\")
        if self.use_booktabs:
            parts.append("\\midrule")
        return parts

    def body(self) -> List[str]:
        out = []
        for r in self.rows:
            cells = [self._format_cell(c) for c in r]
            out.append(" & ".join(cells) + " \\\\")
        return out

    def footer(self) -> List[str]:
        parts = []
        if self.use_booktabs:
            parts.append("\\bottomrule")
        parts.append("\\end{tabular}")
        if self.use_adjustbox:
            parts.append("\\end{adjustbox}")
        parts.append("\\end{table}")
        return parts

    def generate_latex(self) -> str:
        return "\n".join(self.header() + self.body() + self.footer())

    def save(self, path: Union[str, pathlib.Path]) -> None:
        p = pathlib.Path(path)
        p.write_text(self.generate_latex(), encoding="utf8")
