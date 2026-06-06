from typing import List, Iterable, Optional

class AxisScale:
    @staticmethod
    def logy(): return "semilogyaxis"
    @staticmethod
    def logx(): return "semilogxaxis"
    @staticmethod
    def loglog(): return "loglogaxis"
    @staticmethod
    def usual(): return "axis"

class Marks:
    @staticmethod
    def pentagon(): return 'pentagon*'
    @staticmethod
    def square(): return 'square*'
    @staticmethod
    def diamond(): return 'diamond*'
    @staticmethod
    def triangle(): return 'triangle*'

class Line:
    def __init__(self, X, Y,
                 mark: Optional[str] = None,
                 color: Optional[str] = None,
                 legend: Optional[str] = None,
                 custom_options: Optional[str] = ""):
        self.X = X
        self.Y = Y
        self.mark = mark
        self.color = color
        self.legend = legend
        self.custom_options = custom_options or ""

    def get_tikz_options(self):
        options = []
        if self.color:
            options.append(f"color={self.color}")
        if self.mark:
            options.append(f"mark={self.mark}")
        if self.custom_options:
            options.append(self.custom_options)
        return ", ".join(options)

    def get_coordinates(self):
        return " ".join(f"({x}, {y})" for x, y in zip(self.X, self.Y))

    def get_legend(self):
        return f"\\addlegendentry{{{self.legend}}}" if self.legend else ""

    def add_custom_option(self, option: str):
        if self.custom_options:
            self.custom_options += f", {option}"
        else:
            self.custom_options = option

class Axis:
    def __init__(self, lines: Iterable[Line],
                 xlabel: str = "",
                 ylabel: str = "",
                 axistype: str = AxisScale.usual(),
                 xshift: float = 0,
                 yshift: float = 0,
                 title: str = "Title"):
        self.lines = list(lines)
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.xshift = xshift
        self.yshift = yshift
        self.axis = axistype
        self.title = title

    def position_cm(self, xshift, yshift):
        self.xshift = xshift
        self.yshift = yshift

    def add_title(self, title): self.title = title
    def add_line(self, line: Line): self.lines.append(line)
    def setxlabel(self, s: str): self.xlabel = s
    def setylabel(self, s: str): self.ylabel = s
    def setaxistype(self, s): self.axis = s

    def generate_axis_tikz(self) -> str:
        header = rf"""
\begin{{{self.axis}}}[
    title={{{self.title}}},
    xlabel={{{self.xlabel}}},
    ylabel={{{self.ylabel}}},
    xshift={self.xshift}cm,
    yshift={self.yshift}cm,
    legend pos=north west,
    ymajorgrids=true,
    grid style=dashed,
    legend style={{nodes={{scale=0.5, transform shape}}}}
]
"""
        body = ""
        for line in self.lines:
            options = line.get_tikz_options()
            coords = line.get_coordinates()
            legend = line.get_legend()
            body += f"\\addplot[{options}]\n\tcoordinates {{ {coords} }};\n"
            if legend:
                body += f"\t{legend}\n"

        footer = f"\\end{{{self.axis}}}\n"
        return header + body + footer

class GraphData:
    def __init__(self, axes: List[Axis]):
        self.axes = axes

    def add_axis(self, axis: Axis):
        self.axes.append(axis)

    def generate_tikz(self) -> str:
        return "\n".join(axis.generate_axis_tikz() for axis in self.axes)